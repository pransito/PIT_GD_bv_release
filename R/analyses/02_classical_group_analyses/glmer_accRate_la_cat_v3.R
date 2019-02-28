## PREAMBLE ===================================================================
# PIT GD behav study
# classical group-mean-differences tests/mixed ANOVAs (using lme4)
# script to describe acceptance rate
# script to fit glmer models to answer questions on acceptance rate,
# loss aversion, effects of category and group
# model comparison to get significance of overall effect of e.g. group or category
# on gain and loss
# plotting sensitivity to gain, loss and plotting loss aversion per group

# BEFORE RUNNING THIS SCRIPT
# YOU HAVE TO run the R/select_study.R script with which_study = "POSTPILOT_HCPG"

## SETTINGS [nothing to change here ==============================================
# control object for the glmer fitting
cur_control = glmerControl(check.conv.grad="ignore",
                           check.conv.singular="ignore",
                           check.conv.hess="ignore",
                           optCtrl=list(optimizer = "nloptwrap",maxfun=250))
# do the boot / permutation test or just load prepared results?
# attention: it takes over an hours to run the permuation test for loss aversion group comparison
doBootPerm = 0
# bootstrap for plotting cfint's of beta gain and beta loss (already done)
doBoot     = 0
# wd for saving the results of the bootstraps
setwd(root_wd)
setwd('02_classical_group_analyses/results/effects_under_0_la')
if (which_study == 'MRI') {
  setwd('MRI')
} else if (which_study == 'POSTPILOT_HCPG') {
  setwd('behav')
} else {
  stop('which_study is set to an unknown value!')
}
bootResWd = getwd()
# how many bootstraps/permutations?
cur_num   = 100
# how many cpus to use?
cur_cpus  = detectCores()-1
# fit glmer models?
# careful this takes on an intel core i7 with 8GB RAM about 10 minutes
doFitGlmer = 1
# put in the original fixed effect estimation
put_in_original_fe = 1

## glmer models acceptance rate ===============================================
if (doFitGlmer) {
  moda_00  = glmer(accept_reject ~ 1 + (1|subject) + (1|stim) + (1|cat),data = data_pdt,family = 'binomial')
  moda_01  = glmer(accept_reject ~ cat + (cat|subject) + (1|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  moda_02  = glmer(accept_reject ~ cat*HCPG + (cat|subject) + (1|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  moda_01b = glmer(accept_reject ~ HCPG + (1|subject) + (1|stim) + (1|cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
}

## glmer models la ============================================================
if (doFitGlmer) {
  modla_00  = glmer(accept_reject ~ 1 + (1|subject) + (1|stim) + (1|cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_01  = glmer(accept_reject ~ gain + loss  + (gain + loss |subject) + (gain + loss |stim) + (gain + loss |cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_0g  = glmer(accept_reject ~ (gain + loss )*HCPG + (gain + loss |subject) + (gain + loss |stim)  + (gain + loss |cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_c0  = glmer(accept_reject ~ (gain + loss ) + cat + (gain + loss + cat|subject) + (gain + loss |stim) + (gain + loss |cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_cg  = glmer(accept_reject ~ (gain + loss )*HCPG + cat*HCPG + (gain + loss + cat |subject) + (gain + loss |stim)  + (gain + loss |cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_cgi = glmer(accept_reject ~ (gain + loss )*cat*HCPG + ((gain + loss )*cat|subject) + (gain + loss |stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_ci  = glmer(accept_reject ~ (gain + loss )*cat + ((gain + loss )*cat|subject) + (gain + loss |stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
}

## check model fit per subject
cur_dp         = modla_cg@frame
cur_dp$pred_00 = as.numeric(as.numeric(predict(modla_00) >= 0) == cur_dp$accept_reject)
cur_dp$pred_cg = as.numeric(as.numeric(predict(modla_cg) >= 0) == cur_dp$accept_reject)
dens_df        = aggregate(cbind(pred_00,pred_cg) ~ subject + HCPG, data = cur_dp, FUN = mean)
message('The model modla_cg makes better predictions than modla_00 in ...')
message(' ')
message(mean(dens_df$pred_cg > dens_df$pred_00)*100)
message(' ')
message('... percent subjects.')

## acceptance rate and under different cue conditions #########################
if (which_study == 'POSTPILOT_HCPG') {
  # acceptance rate graph, descriptives (CIs over subjects per category)
  # acceptance rate not based on lmer model
  mod_acc        = aggregate(as.numeric(as.character(data_pdt$accept_reject)),by=list(data_pdt$subject,data_pdt$cat), FUN=mean.rmna)
  names(mod_acc) = c('subject','category','mean_acceptance')
  mod_acc$Group  = agk.recode.c(mod_acc$subject,dat_match$VPPG,dat_match$HCPG)
  mod_acc        = aggregate(mod_acc$mean_acceptance,by=list(mod_acc$Group,mod_acc$cat),FUN=agk.boot.ci,R=2000,lower=0.025,upper=0.975,cur_fun=mean)
  mod_acc        = data.frame(mod_acc[[1]], mod_acc[[2]],mod_acc[[3]])
  names(mod_acc) = c('Group','category','mean_acceptance','ci_0025','ci_0975')
  mod_acc$Group  = agk.recode.c(mod_acc$Group,'PG','GD')
  
  mRat  = ggplot(mod_acc, aes(category, mean_acceptance,fill=Group))
  mRat  = mRat + labs(x='category', y=paste('Mean of acceptance (',0.95*100,'% CI, bootstrapped)'))
  mRat  = mRat + ggtitle("Mean acceptance across categories")
  mRat  = mRat + geom_bar(position="dodge", stat="identity")
  dodge = position_dodge(width=0.9)
  mRat  = mRat + geom_bar(position=dodge, stat="identity")
  mRat  = mRat + geom_errorbar(aes(ymin = ci_0025, ymax = ci_0975), position=dodge, width=0.25) + theme_bw()
  print(mRat)
}

## sensitivity to gain, loss, LA ==============================================
laecg_effects = summary(modla_0g)
  
HC_gain = laecg_effects$coefficients['gain',1]
GD_gain = HC_gain + laecg_effects$coefficients['gain:HCPGPG',1]
HC_loss = laecg_effects$coefficients['loss',1]
GD_loss = HC_loss + laecg_effects$coefficients['loss:HCPGPG',1]  

message('Sensitivity to gain, loss for HC, GD was: ')
message('gain')
message(HC_gain)
message(GD_gain)
message('loss')
message(HC_loss)
message(GD_loss)
message('Loss aversion for HC, GD was according to fixed effects of gain and loss: ')
message(-HC_loss/HC_gain)
message(-GD_loss/GD_gain)

print (anova(modla_00,modla_01,modla_0g,modla_cg,modla_cgi))
print(summary(modla_cg))

# acceptance rate only between group
mod_accnc        = aggregate(as.numeric(as.character(data_pdt$accept_reject)),by=list(data_pdt$subject), FUN=mean.rmna)
names(mod_accnc) = c('subject','mean_acceptance')
mod_accnc$Group  = agk.recode.c(mod_accnc$subject,dat_match$VPPG,dat_match$HCPG)
mod_accnc        = aggregate(mod_accnc$mean_acceptance,by=list(mod_accnc$Group),FUN=agk.boot.ci,R=2000,lower=0.025,upper=0.975,cur_fun=mean)
mod_accnc        = data.frame(mod_accnc[[1]],mod_accnc[[2]])
names(mod_accnc) = c('Group','mean_acceptance','ci_0025','ci_0975')
mod_accnc$Group  = agk.recode.c(mod_accnc$Group,'PG','GD')
message('Acceptance rate per group descriptively across all categories.')
print(mod_accnc)

message('effect of category and group')
print(anova(moda_00,moda_01,moda_02)) # effect of category and group

message('stats without cat (simple acceptance rate difference between groups)')
print(anova(moda_00,moda_01b))

## plot loss aversion between groups ==========================================
# bootstrap the CIs
setwd(bootResWd)
library(dplyr)
if (doBoot == 1) {
  # bootstrap p-value modla_0g (permutation)
  effects_under_0_0g = agk.boot.p.mermod(mermod = modla_0g,mermod0 = modla_01,num_cpus = cur_cpus,num = cur_num,fun_extract = fixef,cur_control = cur_control,permvars = c('HCPG'),type='perm')
  save(file= 'effects_under_0_0g_perm_300.RData',list=c('effects_under_0_0g'))
  
  # bootstrap cfint modla_0g (np boot)
  nonp_unit        = 'batch'
  nonp_batch_names = c('HCPG','subject')
  cur_type         = 'non-parametric'
  boot_cfint_0g = agk.boot.cfint.mermod(mermod = modla_0g,num_cpus = cur_cpus,num = cur_num,fun_extract = get_la_fixef_pdt,cur_control = cur_control,
                                        cur_type = cur_type, nonp_unit = nonp_unit,nonp_batch_names = nonp_batch_names,
                                        agk.boot.cfint.mermod.subfun = agk.boot.cfint.mermod.subfun)
  save(file = 'boot_cfint_0g_100_npmu.RData',list=c('boot_cfint_0g'))
}

# graph fixed effects la model
# alternative graph using cfint from fixed effects bootstrap
cie    = new.env()
load('boot_cfint_0g_100_npmu.RData',envir = cie)
rcfint = cie$boot_cfint_0g[[1]]
for (cc in 2:length(cie$boot_cfint_0g)) {
  rcfint = rbind(rcfint,cie$boot_cfint_0g[[cc]])
}

# add the original
#rcfint = rbind(rcfint,get_la_fixef_pdt(modla_0g))

# df of bootstrap data
rcfint           = data.frame(rcfint)
x_la_HCgrPG        = rcfint$x_la_HCgrPG   
rcfint$x_la_HCgrPG = NULL
rcfint_HC          = rcfint[grep('_HC',names(rcfint))]
rcfint_HC$group    = 'HC'
rcfint_PG          = rcfint[grep('_PG',names(rcfint))]
rcfint_PG$group    ='PG'
names(rcfint_HC)   = gsub('_HC','',names(rcfint_HC))
names(rcfint_PG)   = gsub('_PG','',names(rcfint_PG))
rcfint             = rbind(rcfint_HC,rcfint_PG)

rcfinta         = aggregate(.~group,data=rcfint,FUN=agk.mean.quantile,lower=0.025,upper=0.975)
rcfintdf        = data.frame(rcfinta[[1]],rcfinta[[2]])
names(rcfintdf) = c('Group','mean','ci_0025','ci_0975')
rcfintdf$var    = names(rcfinta)[2] 
for (rr in 3:length(rcfinta)) {
  cur_df        = data.frame(rcfinta[[1]],rcfinta[[rr]])
  names(cur_df) = c('Group','mean','ci_0025','ci_0975')
  cur_df$var    = names(rcfinta)[rr] 
  rcfintdf      = rbind(rcfintdf,cur_df)
}
la_overall = rcfintdf

if (put_in_original_fe) {
  # put in the original fe
  obs_fixef      = get_la_fixef_pdt(modla_0g)
  cur_vars       = unique(la_overall$var)
  cur_fe_HC      = obs_fixef[c(1,3,5,7)]
  cur_fe_PG      = obs_fixef[c(2,4,6,8)]
  for (gg in 1:length(cur_vars)) {
    la_overall$mean[la_overall$Group == 'HC' & la_overall$var == cur_vars[gg]] = cur_fe_HC[gg]
    la_overall$mean[la_overall$Group == 'PG' & la_overall$var == cur_vars[gg]] = cur_fe_PG[gg]
  }
}

# actul plotting
la_overall$var = agk.recode.c(la_overall$var,la_overall$var,c('intercept','intercept','gain','gain','loss','loss','la','la'))
la_overall$var = factor(la_overall$var, levels = c('la','gain','loss','intercept'))
mRat           = ggplot(la_overall, aes(Group, mean,fill = var))
mRat           = mRat + labs(x='Group', y=paste('Mean of LA (',0.95*100,'% CI, bootstrapped)'))
mRat           = mRat + ggtitle("Fixed effects for loss aversion parameters per group")
mRat           = mRat + geom_bar(position="dodge", stat="identity")
dodge          = position_dodge(width=0.9)
mRat           = mRat + geom_bar(position=dodge, stat="identity")
mRat           = mRat + geom_errorbar(aes(ymin = ci_0025, ymax = ci_0975), position=dodge, width=0.25) + theme_bw()
print(mRat)

# cfint-based boot-based LA p-value
agk.density_p(x_la_HCgrPG,0)

# p-value for LA based in permuation
cie    = new.env()
load('effects_under_0_0g_perm_300.RData',envir = cie)
rcfint = cie$effects_under_0_0g[[1]]
for (cc in 2:length(cie$effects_under_0_0g)) {
  rcfint = rbind(rcfint,cie$effects_under_0_0g[[cc]])
}
rcfint = data.frame(rcfint)

rcfint$gain.HCPGPG = rcfint$gain + rcfint$gain.HCPGPG
rcfint$loss.HCPGPG = rcfint$loss + rcfint$loss.HCPGPG
rcfint$HCPGPG      = rcfint$X.Intercept. + rcfint$HCPGPG
rcfint$LA          = -rcfint$loss/rcfint$gain
rcfint$LA.HCPGPG   = -rcfint$loss.HCPGPG/rcfint$gain.HCPGPG
rcfint$x_LA_HC_PG  = rcfint$LA - rcfint$LA.HCPGPG

# get the p-value
obs_LA_diff = get_la_fixef_pdt(modla_0g)
message('p-value that fixed effect LA is different between groups under assumption of null-hypothesis (permuted labels):')
print(1-agk.density_p(rcfint$x_LA_HC_PG,obs_LA_diff["x_la_HCgrPG"]))
