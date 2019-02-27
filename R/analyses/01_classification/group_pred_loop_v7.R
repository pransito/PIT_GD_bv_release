# PREAMBLE ====================================================================
# Version 7.0
# script to run group prediction with CV in a loop to see how stable results are
# i.e. group prediction cross validated; in a loop because
# we work with 10-fold or 5-fold cross-validation for tuning and this can be run
# many different times in different ways; needs to be sampled, like a bootstrap
# note: prediction here used interchangeably with "classification"

# through running a control model getting a null-distribution
# also you can run the nooCV version to get the model that is most likely 
# estimated by the given data; nooCV: no outer cross validation; complete
# data is used

# use classifier_setting.R to set what you want to run

# author: Alexander Genauck
# email:  alexander.genauck@charite.de
# date:   14.01.2019

# PREPARATION FOR FOREIGN RUN =================================================
setwd(paste0(root_wd,'/01_classification'))

# LIBRARIES ===================================================================
agk.load.ifnot.install("psych")
agk.load.ifnot.install("pracma")
agk.load.ifnot.install("pls")
agk.load.ifnot.install("Hmisc")
agk.load.ifnot.install("lme4")
agk.load.ifnot.install("reshape2")
agk.load.ifnot.install("R.matlab")
agk.load.ifnot.install("gtools")
agk.load.ifnot.install("plyr")
agk.load.ifnot.install("ggplot2")
agk.load.ifnot.install("rgl")
agk.load.ifnot.install("gridExtra")
agk.load.ifnot.install("boot")
agk.load.ifnot.install("simpleboot")
agk.load.ifnot.install("corrplot")
agk.load.ifnot.install("glmnet")
agk.load.ifnot.install("glmnetUtils")
agk.load.ifnot.install("foreign")
agk.load.ifnot.install("parallel")
agk.load.ifnot.install("foreach")
agk.load.ifnot.install("doSNOW")
agk.load.ifnot.install("GPArotation")
agk.load.ifnot.install("nnet")
agk.load.ifnot.install("msm")
agk.load.ifnot.install("foreign")
agk.load.ifnot.install('readxl')
agk.load.ifnot.install('ptw')
agk.load.ifnot.install('lmPerm')
agk.load.ifnot.install('pROC')
agk.load.ifnot.install('cvTools')
agk.load.ifnot.install('matlib')
agk.load.ifnot.install('robust')
agk.load.ifnot.install('e1071')
agk.load.ifnot.install('compiler')
agk.load.ifnot.install('caret')
#agk.load.ifnot.install('StatMatch')
#agk.load.ifnot.install('MatchIt')
#agk.load.ifnot.install('optmatch')
#agk.load.ifnot.install('WhatIf')
#agk.load.ifnot.install('Matching')
agk.load.ifnot.install('OptimalCutpoints')
agk.load.ifnot.install('dplyr')

# what to run =================================================================
if (!init_run) {
  if (init_done == F) {
    stop('First run the "select_study.R script to initialize all analyses')
  }
  setwd(root_wd)
  setwd('01_classification/')
  source('classifier_settings.R')
}

# PARAMETERS TO SET: General ==================================================
# set some seed to ensure reproducability
des_seed              = 7777 #6993
# run the models (for param extraction in exp)
est_models            = 1
# ridge regression binomial;
# measure for loss function; default is "deviance"; check out options
# help("cv.glmnet")
# in case of gaussian (metric), then 'mse' is used
# 'class' is possible; or 'auc': careful auc needs a certain number of test samples
# I programmed to use k = 5 for innerCV in case of auc
type_measure_binomial = "auc"
# should unit tests for balanced CVfolds be used in wioCV script?
# should always be 1
unit_test_strat_outCV = T
# should unit tests for balanced inner CVfolds be used in wioCV script?
# should always be 1 (see below!)
unit_test_strat_innCV = T
# should scaling be used for innerCV?
# meaning the innerCV training data will be scaled and the scaling params
# will be applied to test data (exp, add feat, convars)
# this is the optional first step of the algorithm (part of the innerCV test)
scale_for_CV          = T
# initial scale: should the whole data set be scaled first
# note that this is a preprocessing step which IS NOT submitted to innerCV testing
# ideally initial_scale and scale_for_CV should either be both on or both off
# both OFF is definitely the recommended version
initial_scale         = F
# standardize in glmnet? coefficients always returned on orig scale
# default is TRUE
des_stand_glmnet      = T
# stratification of inner CV training and(?) test data
# should be T, ideally, but has been not 
strat_innerCV         = T
# estimate the control model (never except below it gets switched on)
c_mod                 = F
# all alphas: then no model selection with just ridge but complete elastic net 
# for all models immediately
all_alphas            = F
# message box width
box_width             = 800
# what predictors to control for
if (which_study == 'MRI') {
  pred_to_control = c('edu_years')
} else if (which_study == 'POSTPILOT_HCPG') {
  pred_to_control = c('smoking_ftdt')
} else {
  stop ('Unknown value of "which_study".')
}
# how many times should full model fit repeated with different folds to get
# mean model? glmnet or svm
fullm_reps            = 10
# set to ML; meaning glmnet elastic net machine learning
which_ML              = 'ML'
# regress out covs (third option; valid option now)
regress_out_covs      = 0
# regress out covs (cleaning MRI data) information criterion
clean_inf_crit        = 'AIC'
# what predictors to clean for
if (which_study == 'MRI') {
  pred_to_clean = c('edu_years')
} else {
  pred_to_clean = c('smoking_ftdt')
}

# PARAMETERS TO SET: Behavior =================================================
# add nonlinear behav models? currently not available; must be 0
add_nonlinear_behav     = F
# repeat the model selection over many kinds of foldings, take modus
ms_reps                 = 10
# ridge the behavioral models after having been fit? (use the ridged versions)
# not recommended; slower computation
ridge_behav_models      = F
# in case of lmlist should residual sum of squares be pooled?
do_pool                 = F
# ridge the behavioral models after having been fit? (fit anew, otherwise load)
# not recommended; slower computation
ridge_behav_models_anew = F
# what alphas to use when ridging those behavioral models?
ridge_bm_alphas         = c(0,0.5,1)
# how many repetitions to ensure stable params?
ridge_bv_reps           = 20
# add physio parameters (pp) in the behavioral models instead of category
# this is an alternative way to model physio PIT (for p.physio paper)
mod_physio_val          = F
# standardize feat vector within subject? (only if at least two values)
behav_within_sub_z      = F
# kitchen sink model: put all variables of all behav models in one
put_all_behav_vars_in   = F

# PARAMETERS TO SET: Physio/MRI ===============================================
# plot rating params
# only for generating plots of ratings/p. physio
plot_ratings         = T
# ... plots for p. physio? not relevant in MRI study, nor in behav only; only for p.physio
plot_physio          = F
# which aggregate fun to be used for rating/physio/MRI variables
agg_fun              = mean.rmna

# PARAMETERS TO SET: General ==================================================
# in glmnet final model which alphas to be tested?
# CAREFUL: MODEL SELECTION PUT ON RIDGE (alpha=0) BY HAND; ALPHAS IGNORED!
# This is to not unselect many or all predictors in predictor sets
alphas          = c(0,0.05,0.1,0.2,0.4,0.8,0.95,1)
# desired_lambdas
# should always be NULL
# however there seems to be a bug:
# https://github.com/lmweber/glmnet-error-example/blob/master/glmnet_error_example.R
# so need to use predefined lambdas, for the case this bug occurs
des_lambdas      = exp(seq(log(0.001), log(10), length.out=100))
# use model matrix (glmnetUtils)
useModelFrame    = T
doGrouped        = F
# how many inner CV folds to be used for picking lambda
# k = 10 recommended, 'LOOCV' is leave-one-out CV
#cnfolds         = "LOOCV"
cnfolds          = 10
if (type_measure_binomial == 'auc') {
  cnfolds        = 5
}
# k of the outer CV
# if left empty then LOOCV will be used
outer_k          = c(10)

# PROCESS PREPS ===============================================================
# DO NOT make any changes here
pred_grp = T
if (strat_innerCV == F) {
  unit_test_strat_innCV = F
}

if (initial_scale == T & scale_for_CV == F) {
  stop('initial_scale is on, then should also be scale_for_CV')
}

# SVM
if(which_ML == 'SVM' & fullm_reps > 3) {
  warning('Since we are using SVM will set fullm_reps to 5')
  fullm_reps = 5
}

# for report use c_mod_report = T
# for p-value against covariate-only classification model
c_mod_report = T

# what features to add
add_cr_pp   = add_cr_pp_ma
add_cr_ra   = add_cr_ra_ma

if (which_study == 'MRI_and_POSTPILOT') {
  add_cr_pp   = 0
  add_cr_ra   = 0
}

# get the operating system
cur_os = Sys.info()
cur_os = cur_os['sysname']
if (cur_os == 'Windows') {
  curpbfun = winProgressBar
  curpbset = setWinProgressBar
} else {
  curpbfun = txtProgressBar
  curpbset = setTxtProgressBar
}

if (outer_cv_addfeaton | noout_cv_addfeaton) {
  message('Careful, the addfeat only option will overwrite results with current addfeat selection under name "phys".')
}

# get the function that wraps the looping of CV rounds, and for initialization
source('group_pred_loop_subf_v7.R')
source('group_pred_init_v7.R')

# run the init of (the CV of) group pred
cur_res = agk.group.pred.init()
agk.assign.envtoenv(cur_res,globalenv())

# get the original matching subjects group
sub_grp_matching = featmod_coefs_bcp[[1]][c('HCPG')]
if (add_cr_pp_ma == T | add_cr_ra_ma == T) {
  stopifnot(all(row.names(featmod_coefs_bcp[[1]]) == row.names(feature_clusters_bcp[[2]])))
}


# just the behavioral parameter sets ==========================================
if (outer_cv_noaddfeat) {
  agk.pred.group.CV(outer_CV = T,addfeat = F,add_cr_pp_fn = F,add_cr_ra_fn = F,des_seed, fm = fm)
}
if (noout_cv_noaddfeat) {
  agk.pred.group.CV(outer_CV = F,addfeat = F,add_cr_pp_fn = F,add_cr_ra_fn = F,des_seed, fm = fm)
}

# behavioral parameter sets plus additional features ==========================
if (outer_cv_wiaddfeat | noout_cv_wiaddfeat) {stopifnot(add_cr_pp_ma | add_cr_ra_ma)}
if (outer_cv_wiaddfeat) {
  agk.pred.group.CV(outer_CV = T,addfeat = T,add_cr_pp_ma,add_cr_ra_ma,des_seed, fm = fm)
}
if (noout_cv_wiaddfeat) {
  agk.pred.group.CV(outer_CV = F,addfeat = T,add_cr_pp_ma,add_cr_ra_ma,des_seed, fm = fm)
}

# only additional features ===================================================
if (outer_cv_addfeaton | noout_cv_addfeaton) {stopifnot(add_cr_pp_ma | add_cr_ra_ma)}
if (outer_cv_addfeaton) {
  agk.pred.group.CV(outer_CV = T,addfeat = T,add_cr_pp_ma,add_cr_ra_ma,des_seed,addfeat_only = T, fm = fm)
}
if (noout_cv_addfeaton) {
  agk.pred.group.CV(outer_CV = F,addfeat = T,add_cr_pp_ma,add_cr_ra_ma,des_seed,addfeat_only = T, fm = fm)
}

# OUTERCV, CONTROL MODEL =====================================================
# the control model; intercept only, or only the control variables
if (outer_cv_c_model) {
  agk.pred.group.CV(outer_CV = T,addfeat=F,add_cr_pp_fn = F,add_cr_ra_fn = F,des_seed,addfeat_only = F,c_mod = T, fm = fm)
}