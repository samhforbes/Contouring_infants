######################################################################
#Analysis of visual closure

#
#Note: Runs analysis on data created with the new
#     script. Final analysis
#####################################################################
#Pre-Processing
#####################################################################
#Load packages
library(tidyr)
library(readr)
library(dplyr)
library(rstan)
library(brms)
library(tidybayes)
library(sjstats)
library(ggsci)
library(ggplot2); theme_set(theme_classic(base_size = 30))
library(eyetrackingR)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = T)

#load file
visclo <- read_csv('data/visual_closure.csv') #obscured
tinytiming <- read_csv('data/tdata.csv') #control experiment
CDI <- read_csv('data/CDI_words.csv') #CDI

cdi2 <- CDI %>% 
  rename(ID = TestID) %>% 
  gather(key = Target, value = Score,
         butterfly, cat, dog, duck, apple, tree)

visclo$Condition <- 'Closure'
tinytiming$Condition <- 'Control'

clodata <- inner_join(visclo, cdi2)

contdata <- inner_join(tinytiming, cdi2)
#2 participants never did CDI

clodata2 <- clodata %>% 
  filter(Score > 0) %>% 
  mutate(Trial = paste(Condition, Trial, sep = ''))

contdata2 <- contdata %>% 
  filter(Score > 0) %>% 
  mutate(Trial = paste(Condition, Trial, sep = ''))
#################################################################
#ERT

CLdata <- make_eyetrackingr_data(clodata2, 
                                 participant_column = "ID",
                                 trial_column = "Trial",
                                 time_column = "Timestamp",
                                 trackloss_column = "Other",
                                 aoi_columns = c('TargetLook','DistLook'),
                                 treat_non_aoi_looks_as_missing = TRUE)

CNdata <- make_eyetrackingr_data(contdata2, 
                                 participant_column = "ID",
                                 trial_column = "Trial",
                                 time_column = "Timestamp",
                                 trackloss_column = "Other",
                                 aoi_columns = c('TargetLook','DistLook'),
                                 treat_non_aoi_looks_as_missing = TRUE)

#########################################################################
#Analyse Trackloss CL
trackloss <- trackloss_analysis(data = CLdata)
trackloss_subjects <- unique(trackloss[,c('ID','TracklossForParticipant')])

#remove files with trackloss
CLdata_clean <- clean_by_trackloss(data = CLdata,
                                            trial_prop_thresh = 0.75)

#(removed 11 trials)
CLtrackloss_clean <- trackloss_analysis(data = CLdata_clean)
CLtrackloss_clean_subjects <- unique(CLtrackloss_clean[, c('ID','TracklossForParticipant')])

#find out trackloss
CLtrackloss_mean <- mean(1 - CLtrackloss_clean_subjects$TracklossForParticipant)
CLtrackloss_sd <- sd(1- CLtrackloss_clean_subjects$TracklossForParticipant)

#See trials contributed by each participant
CLfinal_summary <- describe_data(CLdata_clean, 'TargetLook', 'ID')
CLmean_num_trials <- mean(CLfinal_summary$NumTrials)
CLsd_sum_trials <- sd(CLfinal_summary$NumTrials)

#see summary stats
CLsummaries <- CLdata_clean %>% 
  group_by(ID, Age, Gender) %>% 
  summarise(nTrial = length(Trial)) %>% 
  ungroup() %>% 
  group_by(Age, Gender) %>% 
  summarise(nTrial = mean(nTrial),
            N = length(ID)) %>% 
  ungroup()
#########################################################################
#Analyse Trackloss CN
trackloss <- trackloss_analysis(data = CNdata)
trackloss_subjects <- unique(trackloss[,c('ID','TracklossForParticipant')])

#remove files with trackloss
CNdata_clean <- clean_by_trackloss(data = CNdata,
                                   trial_prop_thresh = 0.75)

#(removed 64 trials, removed 0 participants)
CNtrackloss_clean <- trackloss_analysis(data = CNdata_clean)
CNtrackloss_clean_subjects <- unique(CNtrackloss_clean[, c('ID','TracklossForParticipant')])

#find out trackloss
CNtrackloss_mean <- mean(1 - CNtrackloss_clean_subjects$TracklossForParticipant)
CNtrackloss_sd <- sd(1- CNtrackloss_clean_subjects$TracklossForParticipant)

#See trials contributed by each participant
CNfinal_summary <- describe_data(CNdata_clean, 'TargetLook', 'ID')
CNmean_num_trials <- mean(CNfinal_summary$NumTrials)
CNsd_sum_trials <- sd(CNfinal_summary$NumTrials)

#see summary stats
CNsummaries <- CNdata_clean %>% 
  group_by(ID, Age, Gender) %>% 
  summarise(nTrial = length(Trial)) %>% 
  ungroup() %>% 
  group_by(Age, Gender) %>% 
  summarise(nTrial = mean(nTrial),
            N = length(ID)) %>% 
  ungroup()
#########################################################################
# Setup plot

CL_clean <- subset_by_window(CLdata_clean,
                                   window_start_time = 2000,
                                   window_end_time = 7000,
                                   rezero = T,
                                   remove = T)

CN_clean <- subset_by_window(CNdata_clean,
                             window_start_time = 2000,
                             window_end_time = 5000,
                             rezero = T,
                             remove = T)

clean <- rbind(CL_clean, CN_clean)

timecourse <- make_time_sequence_data(clean,
                                      time_bin_size = 100,
                                      aois = 'TargetLook',
                                      predictor_columns = c('Age', 'Gender', 'Condition'),
                                      summarize_by = 'ID')

timecourse %>%
  mutate(Experiment = ifelse(Condition == 'Closure', 'Experiment 1', 'Experiment 2')) %>% 
  ggplot(aes(x = Time, y = Prop,
             colour = Experiment, linetype = Experiment, shape = Experiment)) +
  stat_summary(geom = 'pointrange', fun.data = 'mean_se') +
  geom_hline(yintercept = 0.5, colour = 'red', linetype = 2) +
  theme(legend.position = c(0.85, 0.85)) +
  xlab('Time (ms)') +
  ylab('Proportion looks to target') +
  #scale_colour_manual(values = c('navyblue', 'turquoise')) +
  scale_color_manual(values = c('forestgreen', 'purple')) +
  ggsave('plots/Overall.pdf', width = 11.69, height = 8.27)
#########################################################################
# analyse Exp1
CL_clean2 <- subset_by_window(CLdata_clean,
                                    window_start_time = 2000,
                                    window_end_time = 5000,
                                    rezero = T,
                                    remove = T) # All 5 not needed

CLtimecourse <- make_time_sequence_data(CL_clean2,
                                       time_bin_size = 100,
                                       aois = 'TargetLook',
                                       predictor_columns = c('Age', 'Gender'),
                                       summarize_by = 'ID')

CLtimecourse <- CLtimecourse %>% 
  filter(SamplesTotal > 0) %>% 
  mutate(Age = as.character(Age)) %>% 
  mutate(Age_C = ifelse(Age == '16', 0.5, -0.5))

#let's take a look
# modtest <- glmmTMB::glmmTMB(cbind(SamplesInAOI, SamplesTotal-SamplesInAOI) ~
#                               (ot1 + ot2 + ot3) *
#                               Age_C +
#                               (ot1+ot2+ot3|ID),
#                             family = binomial,
#                             data = CLtimecourse)
# 
# summary(modtest)
# CLtimecourse$predict <- predict(modtest, type = 'response')
# 
# CLtimecourse %>%
#   ggplot(aes(x = Time, y = Prop,
#              colour = Age, linetype = Age, shape = Age)) +
#   stat_summary(geom = 'pointrange', fun.data = 'mean_se') +
#   stat_summary(aes(y = predict), geom = 'line', fun.y = 'mean') +
#   geom_hline(yintercept = 0.5, colour = 'red', linetype = 2) +
#   theme(legend.position = c(0.85, 0.85)) +
#   xlab('Time (ms)') +
#   ylab('Prop looks to target') +
#   scale_colour_manual(values = c('navyblue', 'turquoise'))
# 
# acf(resid(modtest))

prior <- set_prior('student_t(3, 0, 2.5)', class = 'b') #mean 0, sd of 10
#1 is close to a cauchy recommended by Gelman 2008 
#could try nu of 3 perhaps

clo_brm_short <- brm(SamplesInAOI|trials(SamplesTotal) ~
                       (ot1 + ot2 + ot3) *
                       Age_C +
                       (1 + ot1 + ot2 + ot3|ID),
                     family = binomial('logit'),
                     prior = prior,
                     data = CLtimecourse,
                     warmup = 1000,
                     iter = 2000,
                     control = list(adapt_delta = 0.8, max_treedepth = 30),
                     chains = 4,
                     cores = 4,
                     thin = 1)
save(clo_brm_short, file = 'models/clo_brm_rejigged')

preds <- data.frame(fitted(clo_brm_short, type = 'response'))
summary(clo_brm_short)

CLdatamodel <- cbind(CLtimecourse, preds)

CLdatamodel <- CLdatamodel %>% 
  mutate(Estimate = Estimate/SamplesTotal,
         Est.Error = Est.Error/SamplesTotal,
         Q2.5 = Q2.5/SamplesTotal,
         Q97.5 = Q97.5/SamplesTotal)

se <- function(x) sqrt(var(x)/length(x))

CLshortmodel <- CLdatamodel %>% 
  group_by(Age, Time) %>% 
  summarize(Estimate = mean(Estimate),
            Q2.5 = mean(Q2.5),
            Q97.5 = mean(Q97.5),
            S.E. = se(Prop),
            Prop = mean(Prop))%>% 
  mutate(Condition = 'Closure')

CLshortmodel %>% 
  ggplot(aes(x = Time, y = Prop,
             colour = Age, linetype = Age, shape = Age, fill = Age)) +
  geom_pointrange(aes(ymin = Prop - S.E., ymax = Prop + S.E.), alpha = 0.3) +
  stat_summary(aes(y = Estimate), geom = 'line', fun.y = 'mean', size = 1.5) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), colour = NA, alpha = 0.2) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = 'red') +
  theme(legend.position = c(0.85, 0.85)) +
  coord_cartesian(ylim = c(0.4, 0.75), xlim = c(0, 3000)) +
  xlab('Time (ms)') +
  ylab('Prop looks to target') +
  # scale_colour_manual(values = c('navyblue', 'turquoise')) +
  # scale_fill_manual(values = c('navyblue', 'turquoise')) +
  scale_color_jco() +
  scale_fill_jco() +
  ggsave('plots/BayesExp1_short.pdf', width = 11.69, height = 8.27)



#######################################################################
# Exp2
CN_clean2 <- subset_by_window(CNdata_clean,
                              window_start_time = 2000,
                              window_end_time = 5000,
                              rezero = T,
                              remove = T) 

CNtimecourse <- make_time_sequence_data(CN_clean2,
                                        time_bin_size = 100,
                                        aois = 'TargetLook',
                                        predictor_columns = c('Age', 'Gender'),
                                        summarize_by = 'ID')

CNtimecourse <- CNtimecourse %>% 
  filter(SamplesTotal > 0) %>% 
  mutate(Age = as.character(Age)) %>% 
  mutate(Age_C = ifelse(Age == '16', 0.5, -0.5))

#let's take a look
# modtest2 <- glmmTMB::glmmTMB(cbind(SamplesInAOI, SamplesTotal-SamplesInAOI) ~
#                               (ot1 + ot2 + ot3) *
#                               Age_C +
#                               (ot1+ot2+ot3|ID),
#                             family = binomial,
#                             data = CNtimecourse)
# 
# summary(modtest2)
# CNtimecourse$predict <- predict(modtest2, type = 'response')
# 
# CNtimecourse %>%
#   ggplot(aes(x = Time, y = Prop,
#              colour = Age, linetype = Age, shape = Age)) +
#   stat_summary(geom = 'pointrange', fun.data = 'mean_se') +
#   stat_summary(aes(y = predict), geom = 'line', fun.y = 'mean') +
#   geom_hline(yintercept = 0.5, colour = 'red', linetype = 2) +
#   theme(legend.position = c(0.85, 0.85)) +
#   xlab('Time (ms)') +
#   ylab('Prop looks to target') +
#   scale_colour_manual(values = c('navyblue', 'turquoise'))







prior <- set_prior('student_t(3, 0, 2.5)', class = 'b') #mean 0, sd of 10
#1 is close to a cauchy recommended by Gelman 2008 
#could try nu of 3 perhaps

con_brm_short <- brm(SamplesInAOI|trials(SamplesTotal) ~
                       (ot1 + ot2 + ot3) *
                       Age_C +
                       (1 + ot1 + ot2 + ot3|ID),
                     family = binomial('logit'),
                     prior = prior,
                     data = CNtimecourse,
                     warmup = 1000,
                     iter = 2000,
                     control = list(adapt_delta = 0.8, max_treedepth = 30),
                     chains = 4,
                     cores = 4,
                     thin = 1)
save(con_brm_short, file = 'models/con_brm_rejigged')

preds <- data.frame(fitted(con_brm_short, type = 'response'))
summary(con_brm_short)

CNdatamodel <- cbind(CNtimecourse, preds)

CNdatamodel <- CNdatamodel %>% 
  mutate(Estimate = Estimate/SamplesTotal,
         Est.Error = Est.Error/SamplesTotal,
         Q2.5 = Q2.5/SamplesTotal,
         Q97.5 = Q97.5/SamplesTotal)

se <- function(x) sqrt(var(x)/length(x))

CNshortmodel <- CNdatamodel %>% 
  group_by(Age, Time) %>% 
  summarize(Estimate = mean(Estimate),
            Q2.5 = mean(Q2.5),
            Q97.5 = mean(Q97.5),
            S.E. = se(Prop),
            Prop = mean(Prop))%>% 
  mutate(Condition = 'Closure')

CNshortmodel %>% 
  ggplot(aes(x = Time, y = Prop,
             colour = Age, linetype = Age, shape = Age, fill = Age)) +
  geom_pointrange(aes(ymin = Prop - S.E., ymax = Prop + S.E.), alpha = 0.3) +
  stat_summary(aes(y = Estimate), geom = 'line', fun.y = 'mean', size = 1.5) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), colour = NA, alpha = 0.2) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = 'red') +
  theme(legend.position = c(0.85, 0.85)) +
  coord_cartesian(ylim = c(0.4, 0.75), xlim = c(0, 3000)) +
  xlab('Time (ms)') +
  ylab('Prop looks to target') +
  # scale_colour_manual(values = c('navyblue', 'turquoise')) +
  # scale_fill_manual(values = c('navyblue', 'turquoise')) +
  scale_color_jco() +
  scale_fill_jco() +
  ggsave('plots/BayesExp2_short.pdf', width = 11.69, height = 8.27)

################################################################
# team plot and model 

names <- inner_join(CNtimecourse[1], CLtimecourse[1], by = 'ID')
names <- distinct(names)
CLtimecourse$Age <- as.character(CLtimecourse$Age)
CNtimecourse$Age <- as.character(CNtimecourse$Age)

CL_clean33 <- subset_by_window(CLdata_clean,
                               window_start_time = 2000,
                               window_end_time = 4000,
                               rezero = T,
                               remove = T)

CLtimecourse3 <- make_time_sequence_data(CL_clean33,
                                        time_bin_size = 100,
                                        aois = 'TargetLook',
                                        predictor_columns = c('Age', 'Gender'),
                                        summarize_by = 'ID')

CLtimecourse3 <- CLtimecourse3 %>% 
  filter(SamplesTotal > 0) %>% 
  mutate(Age = as.character(Age)) %>% 
  mutate(Age_C = ifelse(Age == '16', 0.5, -0.5))

CN2 <- left_join(names, CNtimecourse) %>% 
  mutate(Condition = 'Control',
         Cond_C = 0.5)
CL2 <- left_join(names, CLtimecourse3) %>% 
  mutate(Condition = 'Closure',
         Cond_C = -0.5)

combined <- bind_rows(CN2, CL2) %>% 
  filter(Time <= 2000)

# modtest3 <- glmmTMB::glmmTMB(cbind(SamplesInAOI, SamplesTotal-SamplesInAOI) ~
#                               (ot1 + ot2 + ot3) *
#                               Age_C * Cond_C +
#                               (ot1 + ot2 + ot3+ Cond_C|ID),
#                             family = binomial,
#                             data = combined)
# # 
# summary(modtest3)
# # acf(resid(modtest3))
# # 
#  combined$predict <- predict(modtest3, type = 'response')
# # 
# combined %>%
#   ggplot(aes(x = Time, y = Prop,
#              colour = Condition, linetype = Condition, shape = Condition)) +
#   stat_summary(geom = 'pointrange', fun.data = 'mean_se') +
#   stat_summary(aes(y = predict), geom = 'line', fun.y = 'mean') +
#   geom_hline(yintercept = 0.5, colour = 'red', linetype = 2) +
#   theme(legend.position = c(0.85, 0.85)) +
#   xlab('Time (ms)') +
#   ylab('Prop looks to target') +
#   scale_colour_manual(values = c('navyblue', 'turquoise')) +
#   facet_wrap(~Age)

# combined_short <- brm(SamplesInAOI|trials(SamplesTotal) ~
#                        (ot1 + ot2 + ot3) *
#                        Age_C * Cond_C +
#                        (1 + ot1 + ot2 + ot3|ID),
#                      family = binomial('logit'),
#                      prior = prior,
#                      data = combined,
#                      warmup = 1000,
#                      iter = 2000,
#                      control = list(adapt_delta = 0.8, max_treedepth = 30),
#                      chains = 4,
#                      cores = 4,
#                      thin = 1)
# save(combined_short, file = 'models/combined_short')

combined_short2 <- brm(SamplesInAOI|trials(SamplesTotal) ~
                         (ot1 + ot2 + ot3) *
                         Age_C * Cond_C +
                         (1 + ot1 + ot2 + ot3|ID:Condition),
                       family = binomial('logit'),
                       prior = prior,
                       data = combined,
                       warmup = 1000,
                       iter = 2000,
                       control = list(adapt_delta = 0.8, max_treedepth = 30),
                       chains = 4,
                       cores = 4,
                       thin = 1)
save(combined_short2, file = 'models/combined_short2')

preds <- data.frame(fitted(combined_short2, type = 'response'))

summary(combined_short2)

combmodel <- cbind(combined, preds)

combmodel <- combmodel %>% 
  mutate(Estimate = Estimate/SamplesTotal,
         Est.Error = Est.Error/SamplesTotal,
         Q2.5 = Q2.5/SamplesTotal,
         Q97.5 = Q97.5/SamplesTotal)

se <- function(x) sqrt(var(x)/length(x))

combshortmodel <- combmodel %>% 
  group_by(Age, Time, Condition) %>% 
  summarize(Estimate = mean(Estimate),
            Q2.5 = mean(Q2.5),
            Q97.5 = mean(Q97.5),
            S.E. = se(Prop),
            Prop = mean(Prop)) %>% 
  ungroup()

combshortmodel %>% 
  mutate(Experiment = ifelse(Condition == 'Closure', 'Experiment 1', 'Experiment 2')) %>% 
  ggplot(aes(x = Time, y = Prop,
             colour = Experiment, linetype = Experiment, 
             shape = Experiment, fill = Experiment)) +
  geom_pointrange(aes(ymin = Prop - S.E., ymax = Prop + S.E.), alpha = 0.3, linetype = 1) +
  stat_summary(aes(y = Estimate), geom = 'line', fun = 'mean', size = 1.5) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), colour = NA, alpha = 0.2) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = 'red') +
  theme(legend.position = c(0.2, 0.85)) +
  coord_cartesian(ylim = c(0.4, 0.75), xlim = c(0, 2000)) +
  xlab('Time (ms)') +
  ylab('Proportion looks to target') +
  scale_color_manual(values = c('forestgreen', 'purple')) +
  scale_fill_manual(values = c('forestgreen', 'purple')) +
  facet_wrap(~Age) +
  ggsave('plots/BayesComb_short.pdf', width = 11.69, height = 8.27)
################################################################
# Onset

cl <- left_join(names, CL_clean2) %>% 
  mutate(Condition = 'Closure')

cn <- left_join(names, CN_clean2) %>% 
  mutate(Condition = 'Control')

overall <- bind_rows(cl, cn)

response_window_clean3 <- make_eyetrackingr_data(overall, 
                                                 participant_column = "ID",
                                                 trial_column = "Trial",
                                                 time_column = "Timestamp",
                                                 trackloss_column = "Other",
                                                 aoi_columns = c('TargetLook','DistLook'),
                                                 treat_non_aoi_looks_as_missing = TRUE)

response_clean3 <- subset_by_window(response_window_clean3,
                                    window_start_time = 0,
                                    window_end_time = 3000,
                                    rezero = T,
                                    remove = T)

onsets <- make_onset_data(response_clean3, onset_time = 0, 
                          fixation_window_length = 100, 
                          target_aoi='TargetLook')

plot(onsets, predictor_columns = 'Experiment') +
  theme_classic(base_size = 24) +
  theme(legend.position = c(.8,.65)) +
  #facet_wrap(~ Experiment,  labeller = labeller(Experiment = xplab)) +
  #scale_color_manual(values = c('navyblue', 'turquoise')) +
  scale_color_manual(values = c('forestgreen', 'purple')) +
  scale_linetype_manual(labels = c('Target', 'Other'), values = c('solid', 'dotted')) +
  guides(linetype = guide_legend(title = 'AOI'))
onset_switches <- make_switch_data(onsets, 
                                   predictor_columns = c("Condition", 'Age'))

onsetsw <- onset_switches %>% 
  mutate(Age = as.character(Age)) %>% 
  mutate(Age_C = ifelse(Age == '16', 0.5, -0.5)) %>% 
  mutate(Cond_C = ifelse(Condition == 'Closure', 0.5, -0.5)) %>% 
  mutate(FirstAOIC = ifelse(FirstAOI == 'TargetLook', 0.5, -0.5))

plot(onsetsw, predictor_columns = c('Condition', 'Age'))

prior1 <- set_prior("normal(0,800)", class = "b")

model_switches <- brm(FirstSwitch ~ FirstAOIC * Cond_C * Age_C +
                        (1|Trial) + (1|ID),
                      data = onsetsw,
                      prior = prior1,
                      family = gaussian,
                      warmup = 1000,
                      iter = 2000,
                      chains = 4,
                      thin = 1)
save(model_switches, file = 'models/model_switches')

#gamma
model_switches2 <- brm(FirstSwitch ~ FirstAOIC * Cond_C * Age_C +
                         (1|Trial) + (1|ID),
                       data = onsetsw,
                       prior = c(prior(normal(0,800),class="Intercept"),
                                 prior(normal(0,800),class="b"),
                                 prior(gamma(0.01,0.01),class="shape")),
                       family = Gamma(link = 'log'),
                       inits = c(1,2,4,3),
                       warmup = 1000,
                       iter = 2000,
                       chains = 4,
                       thin = 1)
save(model_switches2, file = 'models/model_switches2')
summary(model_switches2)

#better priors
model_switches3 <- brm(FirstSwitch ~ FirstAOIC * Cond_C * Age_C +
                         (1|Trial) + (1|ID),
                       data = onsetsw,
                       prior = c(prior(cauchy(0,10),class="Intercept"),
                                 prior(cauchy(0,2.5),class="b"), #Priors from Gelman 2008
                                 prior(gamma(0.01,0.01),class="shape")),
                       family = Gamma(link = 'log'),
                       inits = c(1,2,4,3),
                       warmup = 1000,
                       iter = 2000,
                       chains = 4,
                       thin = 1)
save(model_switches3, file = 'models/model_switches3')

summary(model_switches3)
launch_shinystan(model_switches3)
################################################################

#################################################################
#plot switches
onset_switches2 <- make_switch_data(onsets, 
                                    predictor_columns = c("Condition", 'Age'),
                                    summarize_by = 'ID')

source('misc/flatviolin.R')

onset_switches2 %>% 
  mutate(Experiment = ifelse(Condition == 'Closure', 'Experiment 1', 'Experiment 2')) %>%
  mutate(Age = as.character(Age)) %>% 
  filter(!is.na(FirstSwitch)) %>% 
  mutate(AOI = ifelse(FirstAOI == 'TargetLook', 'Target', 'Other')) %>% 
  ggplot(aes(y = FirstSwitch, x = Age, fill = AOI, colour = AOI, shape = AOI)) +
  geom_flat_violin(alpha = 0.5, position = position_nudge(x = .2, y = 0), colour = NA) +
  geom_boxplot(width = .2,  outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitter(width = .2), size = 1, alpha = 0.8) +
  facet_wrap(~Experiment) +
  theme(legend.position = c(0.9, 0.8)) +
  xlab('Age group') +
  ylab('Mean switch time (ms)') +
  coord_flip(ylim = c(0,2000)) +
  scale_colour_manual(values = c('red', 'grey2')) +
  scale_fill_manual(values = c('red', 'grey2')) +
  # scale_colour_jco() +
  # scale_fill_jco() +
  ggsave('plots/OnsetSwitches.pdf', width = 11.69, height = 8.27)

######################################################################
#C Check split over chance
plot(CLtimecourse, predictor_column = 'Age')
plot(CNtimecourse, predictor_column = 'Age')
CL_sig19 <- make_cool_splines_data.time_sequence_data(subset(CLtimecourse, Age == 19),
                                                    comp_value = 0.5,
                                                    within_subj = T,
                                                    aoi = 'TargetLook')

CL19_split <- analyze_boot_splines(CL_sig19) 

plot(CL19_split)

CL_sig16 <- make_cool_splines_data.time_sequence_data(subset(CLtimecourse, Age == 16),
                                                      comp_value = 0.5,
                                                      within_subj = T,
                                                      aoi = 'TargetLook')

CL16_split <- analyze_boot_splines(CL_sig16) 

plot(CL16_split)


CN_sig19 <- make_cool_splines_data.time_sequence_data(subset(CNtimecourse, Age == 19),
                                                      comp_value = 0.5,
                                                      within_subj = T,
                                                      aoi = 'TargetLook')

CN19_split <- analyze_boot_splines(CN_sig19) 

plot(CN19_split)

CN_sig16 <- make_cool_splines_data.time_sequence_data(subset(CNtimecourse, Age == 16),
                                                      comp_value = 0.5,
                                                      within_subj = T,
                                                      aoi = 'TargetLook')

CN16_split <- analyze_boot_splines(CN_sig16) 

plot(CN16_split)

#########################################################
#combined posterior draws.
combinedoutput <- add_predicted_draws(newdata = combined, 
                                      model = combined_short2,
                                      scale = 'response') 

combinedoutput2 <- combinedoutput %>% 
  ungroup() %>% 
  mutate(Estimate = .prediction/SamplesTotal) %>% 
  select(ID, Age, Gender, Time, Condition, .draw, Estimate) %>% 
  spread(Condition, Estimate) %>% 
  mutate(Difference = Control - Closure)

combinedoutput3 <- combinedoutput2 %>% 
  group_by(.draw, Time, Age) %>% 
  summarise(Difference = mean(Difference, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(Age, Time) %>% 
  summarise(min_Diff = min(Difference, na.rm = T),
            max_Diff = max(Difference, na.rm = T),
            Diff = mean(Difference, na.rm = T)) %>% 
  ungroup()

combinedoutput3 %>% 
  ggplot(aes(x = Time, y = Diff, fill = Age, colour = Age)) +
  # stat(, alpha = 0.3) +
  #stat_summary(geom = 'line', fun = 'mean', alpha = 0.3) +
  geom_ribbon(aes(ymin = min_Diff, ymax = max_Diff), alpha = 0.2, colour = NA) +
  geom_line(size = 1) +
  # geom_line(alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2, colour = 'red') +
  theme(legend.position = 'none') +
  xlab('Time (ms)') +
  ylab('Difference between conditions') +
  # scale_colour_manual(values = c('navyblue', 'turquoise')) + #NOT TURQUOISE!
  # scale_fill_manual(values = c('navyblue', 'turquoise')) +
  scale_colour_jco() +
  scale_fill_jco() +
  facet_wrap(~Age)

##########################################################
#Exp 1 posterior draws.
ex1output <- add_predicted_draws(newdata = CLtimecourse, 
                                 model = clo_brm_short,
                                 scale = 'response') 

ex1output2 <- ex1output %>% 
  ungroup() %>% 
  mutate(Estimate = .prediction/SamplesTotal) %>% 
  select(ID, Age, Gender, Time, .draw, Estimate) 

ex1output3 <- ex1output2 %>% 
  group_by(.draw, Time, Age) %>% 
  summarise(Estimate = mean(Estimate, na.rm = T)) %>% 
  ungroup() %>% 
  spread(Age, Estimate) %>% 
  mutate(Difference = `19` - `16`) %>% 
  group_by(Time) %>% 
  summarise(min_Diff = min(Difference, na.rm = T),
            max_Diff = max(Difference, na.rm = T),
            Diff = mean(Difference, na.rm = T)) %>% 
  ungroup()

ex1output3 %>% 
  ggplot(aes(x = Time, y = Diff)) +
  # stat(, alpha = 0.3) +
  #stat_summary(geom = 'line', fun = 'mean', alpha = 0.3) +
  geom_ribbon(aes(ymin = min_Diff, ymax = max_Diff), alpha = 0.2, colour = NA) +
  geom_line(size = 1) +
  # geom_line(alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2, colour = 'red') +
  theme(legend.position = c(0.2, 0.85)) +
  xlab('Time (ms)') +
  ylab('Difference between age groups') +
  # scale_colour_manual(values = c('navyblue', 'turquoise')) + #NOT TURQUOISE!
  # scale_fill_manual(values = c('navyblue', 'turquoise')) +
  scale_colour_jco() +
  scale_fill_jco() 
##########################################################