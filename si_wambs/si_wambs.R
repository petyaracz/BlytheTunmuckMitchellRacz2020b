# This is the model fitting procedure for the Blythe et al paper.
# First, we load in data and fix factor type and factor levels.
# Then, we fit models on the comprehension data from the nominal experiment, the entire nominal experiment, and the entire kintax experiment, finally, on the participants of the nom exp who are in the kintax exp.
# For each model, we consider two variants, one with a potential interaction of interest, one without it. 
# We go through the wambs checklist with both variants. More or less.
# Once they LGTM, we compare the two variants using waic and loo. We reloo when needed.
# We keep the model with a better overall elpd. If there is no robust difference in elpd we keep the model with less stuff in it.

## --------- header --------- ##

setwd('~/Github/BlytheTunmuckMitchellRacz2020/')

options(mc.cores=parallel::detectCores())
library(tidyverse)
library(broom)
library(brms)
library(ggthemes)
library(bayesplot)
library(tidybayes)
library(patchwork)

kinap = read_tsv('data/nominal_experiment_data.tsv')

participants = kinap %>% 
  distinct(
    Identifier
  ) %>% 
  mutate(
    count = 1:n(),
    `Participant ID` = paste0('Participant ', count)
  )

kinap = left_join(kinap, participants, by = 'Identifier')

# we need to make sure that various ordered factors are ordered the way we want them.

kinap = kinap %>% 
  filter(
    keep == T,
    !is.na(correct)
  ) %>% 
  mutate(
    `Age band` = case_when(
      `Age Band` == 'A 5_7;9' ~ 'A (5 - 7;9)',
      `Age Band` == 'B 7;10_10;6' ~ 'B (7;10 - 10;6)',
      `Age Band` == 'C 10;7_13;3' ~ 'C (10;7 - 13;3)',
      `Age Band` == 'D 13;4_16' ~ 'D (13;4 - 16)'
    ),
    Task2 = case_when(
      Task == 1 ~ 'comprehension',
      Task == 2 ~ 'production',
      Task == 3 ~ 'altercentric_production'
    ),
    `Age band` = ordered(`Age band`),
    Gender = as.factor(Gender),
    Kin.distance = ordered(`Kin Distance`),
    Genealogical.distance = ordered(`Genealogical distance`),
    Half.sibling.merger = ordered(`Half-sibing merger`),
    Same.sex.sibling.merger = ordered(`Same-sex Sib Merger`),
    Step.kin.merger = as.factor(`Step-kin Merger`),
    Task_name = factor(Task2, levels = c('comprehension', 'production', 'altercentric_production')),
    Question_type = case_when(
      Task_name == 'comprehension' & Task_name2 == 'apposite' ~ 'apposite',
      Task_name == 'comprehension' & Task_name2 == 'inapposite' ~ 'inapposite'
    )
  )

kinap$Task_name3 = ifelse(
  kinap$Task_name == 'altercentricity', 'altercentric production', as.character(kinap$Task_name)
)
kinap$Task_name3 = factor(kinap$Task_name3, levels = c('comprehension', 'production', 'altercentric production'))

# we report on comprehension questions separately

k1 = kinap %>% 
  filter(
    Task_name == 'comprehension',
    !is.na(correct_comp)
  ) %>% 
  droplevels %>% 
  mutate(
    Apposite = Task_name2 == 'apposite',
    Veracity = ifelse(Apposite, 'Apposite question', 'Inapposite question')
  )

kintax = read_csv('../data/Witjpi_Kintax_Data_proc_2019_anonymised.txt')

kintax = kintax %>% 
  mutate(
    `Age band` = case_when(
      Age.bands == 'A 5_7;9' ~ 'A (5 - 7;9)',
      Age.bands == 'B 7;10_10;6' ~ 'B (7;10 - 10;6)',
      Age.bands == 'C 10;7_13;3' ~ 'C (10;7 - 13;3)',
      Age.bands == 'D 13;4_16' ~ 'D (13;4 - 16)',
      Age.bands == 'E 16;1_38' ~ 'E (16;1 - 38)',
    ),
    Condition = factor(Condition, levels = c('male v female', 'dual v paucal', 'cross- v parallel-cousin', 'sibling v non-sibling')),
    Condition2 = ordered(Condition),
    Age.bands = ordered(Age.bands)
  ) %>% 
  filter(
    !is.na(Age.bands) # "Christina Alye Cumaiyi"  "Dave Nulurn Mullumbuk"   "Lucy Thangkirra Tcherna": Dave Nurlurn is too young, Lucy and Christina are not from the same genealogy (Joe)
  )

kinap.av = kinap %>% 
  group_by(Identifier, Age.Band) %>% 
  summarise(kinap.av = mean(na.omit(correct)))

kintax2 = kintax %>% 
  inner_join(kinap.av, by = 'Identifier')

## --------- functions --------- ##

comparePosteriors = function(model1, model2, id1, id2){
  posts = posterior_samples(model1) %>% 
    rownames_to_column() %>% 
    pivot_longer(- rowname, names_to = 'predictor', values_to = 'draw') %>% 
    filter(str_detect(predictor, 'b_')) %>% 
    mutate(model = id1)
  posts2 = posterior_samples(model2) %>% 
    rownames_to_column() %>% 
    pivot_longer(- rowname, names_to = 'predictor', values_to = 'draw') %>% 
    filter(str_detect(predictor, 'b_')) %>% 
    mutate(model = id2) %>% 
    bind_rows(posts)
  
  p1 = posts2 %>% 
    ggplot(aes(x = draw, y = predictor, fill = model)) +
    geom_halfeyeh(alpha = .5) +
    scale_fill_brewer(palette = 'Dark2')
  
return(p1)
}

wambsCheck = function(my_fit, my_data){
  
  my_formula = my_fit$formula
  my_data_r = sample_n(my_data, nrow(my_data))
  
  fit_unif = brm(formula = my_formula, data = my_data, family = bernoulli, save_all_pars = T)
  fit_long = brm(formula = my_formula, data = my_data, family = bernoulli, prior = gelman_hill_binom_priors, save_all_pars = T, iter = 5000)
  fit_rand = brm(formula = my_formula, data = my_data_r, family = bernoulli, prior = gelman_hill_binom_priors, save_all_pars = T)
  
  my_text = paste0(
    paste0('Rhat below 1.1 for fit: ', all(rhat(fit) < 1.1), '\n'),  
    paste0('Rhat below 1.1 with noninf prior: ', all(rhat(fit_unif) < 1.1), '\n'), 
    paste0('Rhat below 1.1 with 5000 iters: ', all(rhat(fit_long) < 1.1), '\n'), 
    paste0('Rhat below 1.1 with rand ord: ', all(rhat(fit_rand) < 1.1)) 
  )
  
  p1 = comparePosteriors(my_fit, fit_unif, 'weak prior', 'noninf prior') + ggtitle('non-informative prior')
  p2 = comparePosteriors(my_fit, fit_long, '2000 iter', '5000 iter') + ggtitle('more iterations')
  p3 = comparePosteriors(my_fit, fit_rand, 'orig', 'reshuffled') + ggtitle('reshuffled data')
  
  p4 = wrap_elements((grid::textGrob(my_text))) + p1 + p2 + p3
  
  return(p4)
}

## --------- weakly informative priors --------- ##

gelman_hill_binom_priors = c(
  prior(student_t(1, 0, 2.5), class = "Intercept"),
  prior(student_t(1, 0, 2.5), class = "b")
)


## --------- fit 1 --------- ##

fit1 = brm(correct_comp ~ Age.Band * Apposite + Kin.distance + (1|Identifier) + (1|Referent_ID), data=k1, family=bernoulli, prior = gelman_hill_binom_priors, save_all_pars = T)

p1 = wambsCheck(fit1, k1)
ggsave('si_wambs/fit1_wambs.pdf', width = 16, height = 8)

fit1b = brm(correct_comp ~ Age.Band + Apposite + Kin.distance + (1|Identifier) + (1|Referent_ID), data=k1, family=bernoulli, prior = gelman_hill_binom_priors, save_all_pars = T, control = list(max_treedepth = 15))

p2 = wambsCheck(fit1b, k1)
ggsave('si_wambs/fit1b_wambs.pdf', width = 16, height = 8)

crit1 = add_criterion(fit1, c('waic', 'loo'), reloo = T)
crit1b = add_criterion(fit1b, c('waic', 'loo'), reloo = T)

loo_compare(crit1, crit1b)
summary(fit1)
save(fit1b, file = 'si_models/fit1_prior.rda')

## --------- fit 2 --------- ##

fit2 = brm(correct ~ Age.Band * Task2 + Genealogical.distance + Same.sex.sibling.merger + Step.kin.merger + Half.sibling.merger + (1|Identifier) + (1|Referent_ID), data=kinap, family=bernoulli, prior = gelman_hill_binom_priors, save_all_pars = T)

p3 = wambsCheck(fit2, kinap)
ggsave(plot = p3, filename = 'si_wambs/fit2_wambs.pdf', width = 16, height = 8)

fit2b = brm(correct ~ Age.Band + Task2 + Genealogical.distance + Same.sex.sibling.merger + Step.kin.merger + Half.sibling.merger + (1|Identifier) + (1|Referent_ID), data=kinap, family=bernoulli, prior = gelman_hill_binom_priors, save_all_pars = T)

p4 = wambsCheck(fit2b, kinap)
ggsave(plot = p4, filename = 'si_wambs/fit2b_wambs.pdf', width = 16, height = 8)

crit2 = add_criterion(fit2, c('waic', 'loo'), reloo = T)
crit2b = add_criterion(fit2b, c('waic', 'loo'), reloo = T)

loo_compare(crit2, crit2b)
summary(fit2b)
save(fit2b, file = 'si_models/fit2_prior.rda')


## --------- fit 3 --------- ##

fit3 = brm(correct ~ Age.bands + Condition + (1|Identifier), data = kintax, family=bernoulli, prior = gelman_hill_binom_priors, save_all_pars = T)

p5 = wambsCheck(fit3, kintax)
ggsave(plot = p5, filename = 'si_wambs/fit3_wambs.pdf', width = 16, height = 8)

fit3b = brm(correct ~ Age.bands * Condition + (1|Identifier), data = kintax, family=bernoulli, prior = gelman_hill_binom_priors, save_all_pars = T)

p6 = wambsCheck(fit3b, kintax)
ggsave(plot = p6, filename = 'si_wambs/fit3b_wambs.pdf', width = 16, height = 8)

crit3 = add_criterion(fit3, c('waic', 'loo'), reloo = T)
crit3b = add_criterion(fit3b, c('waic', 'loo'), reloo = T)

loo_compare(crit3, crit3b)
summary(fit3)
save(fit3, file = 'si_models/fit3_prior.rda')

## --------- fit 4 --------- ##

fit4 = brm(correct ~ Age.bands + Condition + (1|Identifier), data = kintax2, family=bernoulli, prior = gelman_hill_binom_priors, save_all_pars = T)
fit4b = brm(correct ~ kinap.av + Condition + (1|Identifier), data = kintax2, family=bernoulli, prior = gelman_hill_binom_priors, save_all_pars = T)

p7 = wambsCheck(fit4, kintax2)

p8 = wambsCheck(fit4b, kintax2)

ggsave(plot = p6, filename = 'si_wambs/fit3b_wambs.pdf', width = 16, height = 8)

crit4 = add_criterion(fit4, c('waic', 'loo'), reloo = T)
crit4b = add_criterion(fit4b, c('waic', 'loo'), reloo = T)

loo_compare(crit4, crit4b)
summary(fit4)
save(fit4, file = 'si_models/fit4_prior.rda')
save(fit4b, file = 'si_models/fit4b_prior.rda')
