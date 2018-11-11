# analysis of ratings (comparison of cue categires and groups)
# run select_study.R with which_study == 'POSTPILOT_HCPG' first
# the function agk.summarize.models will do a one-group analysis and compare
# all categories and it will do a group comparison (HC vs. GD)
# here: "GD" == "PG"

## LOAD PACKAGES AND DATA =====================================================
agk.load.ifnot.install('lme4')
agk.load.ifnot.install('multcomp')
agk.load.ifnot.install('nloptr')
agk.load.ifnot.install('pracma')

## PUT DATA_PDT in safe-keeping ===============================================
data_pdt_ra_bcp = data_pdt

## PARAMS =====================================================================
which_group           = 'both'
pic_value_vars        = c("arousal","dominance","valence","imageRating1s",
                          "imageRating2s","imageRating3s","imageRating4s")
pic_value_vars_labels = c('Valence', 'Arousal','Dominance', "Elicits_Craving",
                          "Representative_for_Gambles","Representative_for_Negative","Representative_for_Positive")
cur_control           = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE)

## FUNCTIONS ==================================================================
nlopt <- function(par, fn, lower, upper, control) {
  .nloptr <<- res <- nloptr(par, fn, lb = lower, ub = upper, 
                            opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
                                        maxeval = 1000, xtol_abs = 1e-6, ftol_abs = 1e-6))
  list(par = res$solution,
       fval = res$objective,
       conv = if (res$status > 0) 0 else res$status,
       message = res$message
  )
}

agk.lme.summary = function(model,type) {
  # will produce no approximation (default),
  # normal distr. approximation
  # to get p-values for lme4 models
  # model: lme4 model
  # type: 'none','norm'
  # 'satter' out cause anoying summary function, takes long
  # and does not allow REML
  # 'ken' also slow
  # https://en.wikipedia.org/wiki/Restricted_maximum_likelihood
  # value: returns the normal summary, and if type
  # is either 'norm', coefficient table with p-value
  
  if (type == 'none') {
    return(summary(model))
  } else if (type == 'norm') {
    cur_sum = summary(model)
    # extract coefficients
    coefs <- data.frame(coef(cur_sum))
    # use normal distribution to approximate p-value
    coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
    cur_sum$coefficients = coefs
    return(cur_sum)
  }
}

agk.estimate.models = function(des_var,cur_control) {
  # des_var is a string that indicates the variable
  # you would like to do your tests on
  
  # define commands
  cmd_1 = paste0('mod0  = lmer(',des_var,' ~ 1 + (1| subject) + (1|stim),data = data_pdt,REML = F)')
  cmd_2 = paste0('modc  = lmer(',des_var,' ~ (0+cat) + (0+cat| subject) ,data = data_pdt, REML = F,control=cur_control)')
  cmd_3 = paste0('modcg = lmer(',des_var,' ~ (0+cat)*HCPG + (0+cat| subject) ,data = data_pdt, REML = F,control=cur_control)')
  
  # evaluate commands
  eval(parse(text = cmd_1))
  eval(parse(text = cmd_2))
  eval(parse(text = cmd_3))
  
  # return
  return(list(mod0 = mod0,modc = modc,modcg = modcg))
} 

agk.summarize.models = function(est_mods) {
  # function that summarizes two lmer mods
  
  disp('####################')
  disp('# MODEL COMPARISON #')
  disp('####################')
  print(anova(est_mods$mod0,est_mods$modc,est_mods$modcg))
  disp('')
  disp('####################')
  disp('# MODEL MOD0       #')
  disp('####################')
  print(agk.lme.summary(est_mods$mod0,type='norm'))
  disp('')
  disp('####################')
  disp('# MODEL MODC       #')
  disp('####################')
  print(agk.lme.summary(est_mods$modc,type='norm'))
  disp('####################')
  disp('# MODEL MODCG      #')
  disp('####################')
  print(agk.lme.summary(est_mods$modcg,type='norm'))
  
  disp('')
  disp('#########################')
  disp('# ALL COMPARISONS MODC  #')
  disp('#########################')
  cons = glht(est_mods$modc, linfct = mcp(cat = "Tukey"))
  print(summary(cons,test=adjusted('none')))
}

## PREP THE DATA ==============================================================
# get a valid HCPG variable
for (hh in 1:length(data_pdt$HCPG)) {
  
  if (any(c('HC','PG') %in% data_pdt$HCPG[hh])) {
    cur_res = data_pdt$HCPG[hh]
  } else {
    cur_res = 'HC'
  }
  data_pdt$HCPG[hh] = cur_res  
}

# do it for HC or PG? PG or everything else (== HC)
if (any(c('HC','PG') %in% which_group)) {
  data_pdt = subset(data_pdt,HCPG == which_group)
} else if(which_group == 'both') {
  # do both
} else {
  stop('Do not no this value of which_group.')
}


# debug
#data_pdt = subset(data_pdt, Cohort == 'PhysioPilot')
#data_pdt = data_pdt[grep(pattern = '^5',x = data_pdt$stim),]

# CATEGORY LABELS =============================================================
# Main effect of the final experimental categories: gam, pos, neg, neu_aw
# data_pdt_finCat = data_pdt
# data_pdt_finCat$cat = factor(as.numeric(as.character(data_pdt_finCat$cat)),levels = c(1,2,3,4,5,6,7,8),
#                              labels = c('gam','neg', 'pos','neuIAPS_NAPS','gray','neuIAPSAW','negIAPS','posIAPS'))
# data_pdt = data_pdt_finCat

if(sum(is.na(data_pdt$cat))) {
  stop('There are NAs in the data_pdt$cat variable!')
}

## CLEANING RATINGS DATA ======================================================
# get rid of all subjects that just haven't rated anything at all
sub_ok = c()
all_subs = unique(data_pdt$subject)
for (ss in 1:length(all_subs)) {
  cur_dat  = data_pdt[data_pdt$subject == all_subs[ss],]
  cur_dat  = cur_dat[pic_value_vars]
  if(all(is.na(cur_dat))) {
    sub_ok[ss] = F
  } else {
    sub_ok[ss] = T
  }
}
data_pdt = data_pdt[data_pdt$subject %in% all_subs[sub_ok],]


# get rid of double stimuli within subjects
all_subs = unique(data_pdt$subject)
cur_dat  = data_pdt[data_pdt$subject == all_subs[1],]
cur_dat  = cur_dat[!duplicated(cur_dat$stim),]
dpdt_cln = cur_dat
for (ss in 2:length(all_subs)) {
  cur_dat  = data_pdt[data_pdt$subject == all_subs[ss],]
  cur_dat  = cur_dat[!duplicated(cur_dat$stim),]
  dpdt_cln = rbind(dpdt_cln,cur_dat) 
}
data_pdt = dpdt_cln

# get rid of stimuli not being rated at all
stim_ok = c()
all_stim = unique(data_pdt$stim)
for (ss in 1:length(all_stim)) {
  cur_dat  = data_pdt[data_pdt$stim == all_stim[ss],]
  cur_dat  = cur_dat[pic_value_vars]
  if(all(is.na(cur_dat))) {
    stim_ok[ss] = F
  } else {
    stim_ok[ss] = T
  }
}
data_pdt = data_pdt[data_pdt$stim %in% all_stim[stim_ok],]


## TEST RATINGS HYPOTHESES ONE GROUP ==========================================
for (dd in 1:length(pic_value_vars_labels)) {
  message(paste0('Univariate tests for the rating variable: ', pic_value_vars_labels[dd]))
  des_var = pic_value_vars[dd]
  res     = agk.estimate.models(des_var,cur_control)
  sumres  = agk.summarize.models(res) 
}

## DATA_PDT into initial state ================================================
data_pdt_ra_bcp = data_pdt