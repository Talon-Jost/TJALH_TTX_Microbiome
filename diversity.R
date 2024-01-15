rm(list=ls())
library(usethis)
library(dplyr)
library(MASS)
library(ggplot2)
library(lme4)
library(nlme)
library(vegan)
library(readr)
library(plotrix)
library(reticulate)
library(effects)
library(remotes)
library(MuMIn)
library(ggtext)
library(car)
library(MuMIn)
library(minpack.lm)
#remotes::install_github('ProcessMiner/nlcor')
library(nlcor)
#remotes::install_github('femiguez/nlraa')
library(nlraa)
library(mgcv)
library(tidymv)
#remotes::install_github("gavinsimpson/gratia")
library(gratia)
library(lefser)
library(devtools)

setwd('/home/tally/Documents/Thesis_Part_2/Newt_TTX_Project/R_Stuff')
list.files()
df_diversity <- read.csv('TJALH_Metadata_Diversity.csv')
names(df_diversity)[22] = 'Whole_TTX'
df_diversity$Whole_TTX = as.numeric(df_diversity$Whole_TTX)
df_diversity$Mass_g = as.numeric(df_diversity$Mass_g)
df_diversity$Bd_Presence = as.numeric(df_diversity$Bd_Presence)
df_diversity$shannon = as.numeric(df_diversity$shannon)

bd_presence_overall_multivariate <- glmer(data = df_diversity, Bd_Presence ~ Whole_TTX * Mass_g * Sex * shannon + (1 | Site), family = binomial)
summary(bd_presence_overall_multivariate)
Anova(bd_presence_overall_multivariate, type = 3)
plot(allEffects(bd_presence_overall_multivariate), multiline = T, ci.style = 'band')

#trial output
#library(sjPlot)
#library(sjmisc)
#library(sjlabelled)
#devtools::install_github("strengejacke/strengejacke")
#sjp.glmer(bd_presence_overall_multivariate)
#sjp.aov1()

library(tibble)
library(flextable)

diversity_log_bd = log(positive$ZE.avg.Bd + 1)
observed_trial <- lm(data = positive, diversity_log_bd ~ Observed_Features)
summary(observed_trial)

positive <- subset(df_diversity, Site %in% (c('GRP', 'SLNF', 'TMC', 'WBC')))
bd_presence_overall_multivariate <- glmer(data = positive, Bd_Presence ~ Whole_TTX * Mass_g * Sex * shannon + (1 | Site), family = binomial)
summary(bd_presence_overall_multivariate)
plot(allEffects(bd_presence_overall_multivariate), multiline = T, ci.style = 'band')
plot(allEffects(bd_presence_overall_multivariate), multiline = T, ci.style = 'band', type = 'response')
table(positive$Bd_Presence, positive$Sex)
anova1 <- Anova(bd_presence_overall_multivariate, type = 3)


multivariate_2 <- glm(data = positive, Bd_Presence ~ Whole_TTX * Mass_g * Sex * shannon, family = quasibinomial)
summary(multivariate_2)
Anova(multivariate_2, type = 3)

multivariate_3 <- glm(data = positive, Bd_Presence ~ Whole_TTX * Mass_g * Sex * shannon * Site, family = quasibinomial)
summary(multivariate_3)
Anova(multivariate_3, type = 3)

tmc_absent <- subset(df_diversity, Site %in% (c('GRP', 'SLNF', 'WBC')))
tmc_multi <- glmer(data = tmc_absent, Bd_Presence ~ Whole_TTX * Mass_g * Sex * shannon + (1 | Site), family = binomial)
summary(tmc_multi)

multivariate_4 <- glm(data = tmc_absent, Bd_Presence ~ Whole_TTX * Mass_g * Sex * shannon * Site, family = quasibinomial)
summary(multivariate_4)
Anova(multivariate_4, type = 3)
Anova(multivariate_4, type = 2)

#do the same with other diversity metrics

sex_glmer <- glmer(data = df_diversity, Bd_Presence ~ Whole_TTX * Sex + (1 | Site), family = binomial)
summary(sex_glmer)
Anova(sex_glmer, type = 3)
plot(allEffects(sex_glmer))

mass_sex_anova = aov(df_diversity$Bd_Presence ~ df_diversity$Sex)
summary(mass_sex_anova)

#10323 start
no_na_diversity <- na.omit(df_diversity)
no_na_diversity$faith <- as.numeric(no_na_diversity$faith)
no_na_diversity$Observed_Features <- as.numeric(no_na_diversity$Observed_Features)
no_na_diversity$pielou_e <- as.numeric(no_na_diversity$pielou_e)
multivariate_5 <- glmer(data = no_na_diversity, Bd_Presence ~ Whole_TTX * Sex * faith + (1 | Site), family = binomial)
summary(multivariate_5)
Anova(multivariate_5, type = 3)
plot(allEffects(multivariate_5), multiline = T, ci.style = 'band')


multivariate_6 <- glmer(data = no_na_diversity, Bd_Presence ~ Whole_TTX * Sex * Observed_Features + (1 | Site), family = binomial)
summary(multivariate_6)
Anova(multivariate_6, type = 3)
plot(allEffects(multivariate_6), multiline = T, ci.style = 'band')

multivariate_7 <- glmer(data = no_na_diversity, Bd_Presence ~ Whole_TTX * Sex * Mass_g * faith + (1 | Site), family = binomial)
summary(multivariate_7)
Anova(multivariate_7, type = 3)
plot(allEffects(multivariate_7), multiline = T, ci.style = 'band', rescale.axis = F)

multivariate_10 <- glmer(data = no_na_diversity, Bd_Presence ~ Whole_TTX * Mass_g * faith + (1 | Site), family = binomial)
summary(multivariate_10)
Anova(multivariate_10, type = 3)
plot(allEffects(multivariate_10), multiline = F, ci.style = 'band')
#Males
Males <- subset(no_na_diversity, Sex == 'M')

multivariate_8 <- glmer(data = Males, Bd_Presence ~ Whole_TTX * Mass_g * faith + (1 | Site), family = binomial)
summary(multivariate_8)
Anova(multivariate_8, type = 3)
plot(allEffects(multivariate_8))

table(df_diversity$Bd_Presence, df_diversity$Sex)
boxplot(df_diversity$Whole_TTX ~ df_diversity$Sex)
ggplot(no_na_diversity, aes(x = Sex, y = Whole_TTX)) + theme_bw() + geom_jitter(width = 0.2)
ggplot(no_na_diversity, aes(x = Sex, y = Mass_g)) + geom_jitter(width = 0.2) + theme_bw()
ggplot(no_na_diversity, aes(x = Sex, y = Observed_Features)) + geom_jitter(width = 0.3) + theme_bw()
ggplot(no_na_diversity, aes(x = Whole_TTX, y = Bd_Presence, shape = Sex, color = Mass_g)) + geom_point() + theme_bw()
ggplot(no_na_diversity, aes(x = Whole_TTX, y = Mass_g, color = Sex)) + geom_point()
#Females
Females <- subset(no_na_diversity, Sex == 'F')
high_females <- subset(Females, Whole_TTX > 2)
multivariate_9 <- glmer(data = Females, Bd_Presence ~ Whole_TTX * Mass_g * faith + (1 | Site), family = binomial)
summary(multivariate_9)
Anova(multivariate_9, type = 3)
plot(allEffects(multivariate_9))

#observed________
new_df <- na.omit(positive)
new_df$ZE.avg.Bd <- as.numeric(new_df$ZE.avg.Bd)
log_bd <- log(new_df$ZE.avg.Bd + 1)
new_df$faith <- as.numeric(new_df$faith)
new_df$Observed_Features <- as.numeric(new_df$Observed_Features)
observed_lm <- lm(data = new_df, log_bd ~ Observed_Features)
summary(observed_lm)

faith_lm <- lm(data=new_df, log_bd ~ faith)
summary(faith_lm)

positive$ZE.avg.Bd <- as.numeric(positive$ZE.avg.Bd)
pos_log_bd <- log(positive$ZE.avg.Bd + 1)
positive_observed_lm <- lm(pos_log_bd ~ Observed_Features, data = positive)
summary(positive_observed_lm)

positive_faith_lm <- lm(pos_log_bd ~ faith, data = positive)
summary(positive_faith_lm)

#____LEFSE
library(indicspecies)
indval=multipatt
freqtable <- read.csv('updated_frequency_table-lvl6.csv', header = TRUE)
dim(freqtable)
nrow(freqtable)

row.names(freqtable)=freqtable[,1]
freqtable=freqtable[,-1]
#classes <- unlist(strsplit(row.names(freqtable), '\\;D_3__'))[seq(2,2*1256,2)]
#row.names(freqtable) = classes
freqtable=t(freqtable)
freqtable=freqtable[-89,]
dim(freqtable)
#colnames(freqtable)=classes
#colnames(freqtable)
site <- unlist(strsplit(row.names(freqtable), '\\.'))[seq(1,176,2)]
indval = multipatt(freqtable,site)
result = summary(indval)
print

table2 = read.csv('updated_frequency_table-lvl6.csv', header = T)
sampleid = names(table2)[2:89]
sampleid = gsub('\\.','\\-',sampleid[1:88])
mysamples=NULL
for (i in 1:88) mysamples=c(mysamples,which(metadata$SwabID==sampleid[i]))
metadata = read.csv('TJAH_TTX_DATA_CLEANED.csv', header = T)

#lvl seven
level_7 <- read.csv('r_frequancy_table_lvl_7.csv', header = T)
row.names(level_7)=level_7[,1]
level_7=level_7[,-1]
#classes <- unlist(strsplit(row.names(freqtable), '\\;D_3__'))[seq(2,2*1256,2)]
#row.names(freqtable) = classes
level_7=t(level_7)
level_7=level_7[-89,]
dim(level_7)
vl_site <- unlist(strsplit(row.names(level_7), '\\.'))[seq(1,176,2)]
indval = multipatt(level_7, vl_site)
result = summary(indval)
print
