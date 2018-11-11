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
# date:   09.11.2018

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
agk.load.ifnot.install('StatMatch')
agk.load.ifnot.install('MatchIt')
agk.load.ifnot.install('optmatch')
agk.load.ifnot.install('WhatIf')
agk.load.ifnot.install('Matching')
agk.load.ifnot.install('OptimalCutpoints')
agk.load.ifnot.install('caret')

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
des_seed              = 6993
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
if (which_study == 'MRT') {
  pred_to_control = c()
} else {
  pred_to_control = c('smoking_ftdt')
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
if (which_study == 'MRT') {
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

if (which_study == 'MRT_and_POSTPILOT') {
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
  agk.pred.group.CV(outer_CV = T,addfeat = F,add_cr_pp_fn = F,add_cr_ra_fn = F,des_seed)
}
if (noout_cv_noaddfeat) {
  agk.pred.group.CV(outer_CV = F,addfeat = F,add_cr_pp_fn = F,add_cr_ra_fn = F,des_seed)
}

# behavioral parameter sets plus additional features ==========================
if (outer_cv_wiaddfeat | noout_cv_wiaddfeat) {stopifnot(add_cr_pp_ma | add_cr_ra_ma)}
if (outer_cv_wiaddfeat) {
  agk.pred.group.CV(outer_CV = T,addfeat = T,add_cr_pp_ma,add_cr_ra_ma,des_seed)
}
if (noout_cv_wiaddfeat) {
  agk.pred.group.CV(outer_CV = F,addfeat = T,add_cr_pp_ma,add_cr_ra_ma,des_seed)
}

# only additional features ===================================================
if (outer_cv_addfeaton | noout_cv_addfeaton) {stopifnot(add_cr_pp_ma | add_cr_ra_ma)}
if (outer_cv_addfeaton) {
  agk.pred.group.CV(outer_CV = T,addfeat = T,add_cr_pp_ma,add_cr_ra_ma,des_seed,addfeat_only = T)
}
if (noout_cv_addfeaton) {
  agk.pred.group.CV(outer_CV = F,addfeat = T,add_cr_pp_ma,add_cr_ra_ma,des_seed,addfeat_only = T)
}

# OUTERCV, CONTROL MODEL =====================================================
# the control model; intercept only, or only the control variables
if (outer_cv_c_model) {
  agk.pred.group.CV(outer_CV = T,addfeat=F,add_cr_pp_ma = F,add_cr_ra_ma = F,des_seed,addfeat_only = F,c_mod = T)
}

# REPORTING: PREPARATION =====================================================
if (do_report) {
  cur_home = getwd()
  setwd(root_wd)
  setwd('01_classification/results/')
  setwd(as.character(runs))
  if (pred_grp) {
    all_res_files = dir(pattern = paste0(which_study,'_predGrp1'))
  } else {
    all_res_files = dir(pattern = paste0(which_study,'_predGrp0'))
  }
  
  if (do_report_with_added_feat) {
    all_res_files_Ha = all_res_files[grep('wiaddfeat',all_res_files)]
  } else if (do_report_no_added_feat) {
    all_res_files_Ha = all_res_files[grep('noaddfeat',all_res_files)]
  } else if (do_report_feat_only) {
    all_res_files_Ha = all_res_files[grep('onlyPhys',all_res_files)]
  } else {
    stop('Ha not defined for report.')
  }
  
  for (ii in 1:length(all_res_files)) {
    load(all_res_files[ii])
  }
  
  # over-load the correct Ha
  for (ii in 1:length(all_res_files_Ha)) {
    load(all_res_files_Ha[ii])
  }
  
  setwd(cur_home)
  if (c_mod_report & report_CV_p) {
    # using control model, renaming variables
    CVp_res_list    = CVcm_res_list
    CVp_res_list_op = CVcm_res_list
  }
  
  # get the current relevant CV_res_list
  if (exists('CV_res_list_op')) {
    CV_res_list = CV_res_list_op
  }
}



# REPORTING: MAIN CODE ========================================================
if (do_report_no_added_feat | do_report_with_added_feat | do_report_feat_only) {
  if (report_CV_p) {
    # p-value for the algorithm
    # get the accuracy null-distrib/alt hypothesis distribution
    cor    = c() # correlation
    cor_p  = c() # correlation control
    coru   = c() # correlation unpooled (mean of k-fold)
    coru_p = c() # correlation unpooled (mean of k-fold) control
    mse    = c() # mean squared error
    mse_p  = c() # mean squared error control
    accs_p = c()
    accs   = c()
    sens_p = c()
    sens   = c()
    spec_p = c()
    spec   = c()
    auc_p  = c()
    auc    = c()
    roc_pl = list()
    rocl   = list()
    osens    = c()   # optimal sensitivity (optimal cutoff post-hoc)
    ospec    = c()   # optimal specificity (optimal cutoff post-hoc)
    osens_p  = c()   # optimal sensitivity (optimal cutoff post-hoc)
    ospec_p  = c()   # optimal specificity (optimal cutoff post-hoc)
    prec     = c()   # precision rate (of the detected how many are actually detected ones)
    prec_p   = c()   # precision under 0
    bacc     = c()   # balanced accuracy
    bacc_p   = c()   # b.a. under 0
    if (length(CVp_res_list) >= 1000) {
      bootstr_length = 1000
    } else {
      bootstr_length = length(CVp_res_list)
    }

    for (ii in 1:bootstr_length) {
      accs_p[ii]   = CVp_res_list[[ii]]$acc$accuracy
      accs[ii]     = CV_res_list[[ii]]$acc$accuracy
      sens_p[ii]   = CVp_res_list[[ii]]$acc$cur_sens
      sens[ii]     = CV_res_list[[ii]]$acc$cur_sens
      spec_p[ii]   = CVp_res_list[[ii]]$acc$cur_spec
      spec[ii]     = CV_res_list[[ii]]$acc$cur_spec
      roc_pl[[ii]] = CVp_res_list[[ii]]$roc
      auc_p[ii]    = as.numeric(roc_pl[[ii]]$auc)
      rocl[[ii]]   = CV_res_list[[ii]]$roc
      auc[ii]      = as.numeric(rocl[[ii]]$auc)
      cor[[ii]]    = CV_res_list[[ii]]$cor$estimate
      cor_p[[ii]]  = CVp_res_list[[ii]]$cor$estimate
      coru[[ii]]   = CV_res_list[[ii]]$coru
      coru_p[[ii]] = CVp_res_list[[ii]]$coru
      mse[[ii]]    = CV_res_list[[ii]]$mse
      mse_p[[ii]]  = CVp_res_list[[ii]]$mse
      bacc[ii]     = mean(c(spec[ii], sens[ii]))
      bacc_p[ii]   = mean(c(spec_p[ii], sens_p[ii]))
    }
    
    # make new ROC objects
    for (ii in 1:bootstr_length) {
      # Ha
      predictor  = rocl[[ii]]$original.predictor
      response   = rocl[[ii]]$original.response
      # get the precision
      cur_cm    = confusionMatrix(as.factor(response),as.factor(ifelse(predictor>=0, 'PG','HC')))
      prec[ii]  = as.numeric(cur_cm$byClass['Precision'])
      response   = ifelse(response == 'PG',1,0)
      predictor  = agk.scale_range(predictor,0,1)
      # get the optimal sen and spe
      cur_df    = data.frame(X = rocl[[ii]]$original.predictor, status = rocl[[ii]]$original.response) 
      cur_o     = optimal.cutpoints(X = "X", status = "status", data =cur_df,methods ="MaxSpSe",tag.healthy = 'HC')
      cur_sum   = summary(cur_o)
      osens[ii] = cur_sum$MaxSpSe$Global$optimal.cutoff$Se[1]
      ospec[ii] = cur_sum$MaxSpSe$Global$optimal.cutoff$Sp[1]
      
      rocl[[ii]] = roc(response,predictor)
      auc[ii]    = auc(rocl[[ii]])
      rocl[[ii]] = smooth( rocl[[ii]],n = 500)
      
      # control model
      predictor    = roc_pl[[ii]]$original.predictor
      response     = roc_pl[[ii]]$original.response
      # get the precision
      cur_cm     = confusionMatrix(as.factor(response),as.factor(ifelse(predictor>=0, 'PG','HC')))
      prec_p[ii] = as.numeric(cur_cm$byClass['Precision'])
      response     = ifelse(response == 'PG',1,0)
      predictor    = agk.scale_range(predictor,0,1)
      # get the optimal sen and spe
      cur_df    = data.frame(X = roc_pl[[ii]]$original.predictor, status = roc_pl[[ii]]$original.response) 
      cur_o     = optimal.cutpoints(X = "X", status = "status", data =cur_df,methods ="MaxSpSe",tag.healthy = 'HC')
      cur_sum   = summary(cur_o)
      osens_p[ii] = cur_sum$MaxSpSe$Global$optimal.cutoff$Se[1]
      ospec_p[ii] = cur_sum$MaxSpSe$Global$optimal.cutoff$Sp[1]
      
      roc_pl[[ii]] = roc(response,predictor)
      auc_p[ii]    = auc(roc_pl[[ii]])
      
      if (var(roc_pl[[ii]]$original.predictor) == 0) {
        cur_obj = list()
        cur_obj$sensitivities = seq(0,1,length.out = 502)
        cur_obj$specificities = seq(0,1,length.out = 502)
        roc_pl[[ii]]          = cur_obj
      } else {
        roc_pl[[ii]] = smooth(roc_pl[[ii]],n = 500)
      }
    }
    
    # print the mean of acc, sens, spec
    acc_sens_spec       = data.frame(auc,accs,sens,spec,osens,ospec,prec,bacc)
    acc_sens_spec_p     = data.frame(auc_p,accs_p,sens_p,spec_p,osens_p,ospec_p,prec_p,bacc_p)
    med_ac_se_sp        = lapply(acc_sens_spec,FUN=mean)
    med_ac_se_sp        = lapply(acc_sens_spec,FUN=mean)
    med_ac_se_sp_np_ci  = lapply(acc_sens_spec,FUN=agk.mean.quantile,
                                 lower=0.025,upper=0.975)
    med_ac_se_sp_npp_ci = lapply(acc_sens_spec_p,FUN=agk.mean.quantile,
                                 lower=0.025,upper=0.975)
    disp('Mean (CI) AUC, accuracy, sensitivity, specificity, opt sens, opt spec, precision, balanced acc. across CV rounds.')
    print(med_ac_se_sp_np_ci )
    disp('UNDER H0: Mean (CI) AUC, accuracy, sensitivity, specificity, opt sens, opt spec, precision, balanced acc. across CV rounds.')
    print(med_ac_se_sp_npp_ci )
    
    # p-values
    # auc
    diffs_auc    = auc - auc_p
    cur_p        = agk.density_p.c(diffs_auc,0)
    disp("Probability that AUC has improved under assumption of 0-hypothesis.")
    print(cur_p)
    # acc
    diffs_accs   = accs - accs_p
    cur_p        = agk.density_p.c(diffs_accs,0)
    disp("Probability that acc has improved under assumption of 0-hypothesis.")
    print(cur_p)
    # sens
    diffs_sens   = sens - sens_p
    cur_p        = agk.density_p.c(diffs_sens,0)
    disp("Probability that sens has improved under assumption of 0-hypothesis.")
    print(cur_p)
    # spec
    diffs_spec   = spec - spec_p
    cur_p        = agk.density_p(diffs_spec,0)
    disp("Probability that spec has improved under assumption of 0-hypothesis.")
    print(cur_p)
    # osens
    diffs_osens  = osens - osens_p
    cur_p        = agk.density_p(diffs_osens,0)
    disp("Probability that optimal sens has improved under assumption of 0-hypothesis.")
    print(cur_p)
    # ospec
    diffs_ospec  = ospec - ospec_p
    cur_p        = agk.density_p(diffs_ospec,0)
    disp("Probability that optimal spec has improved under assumption of 0-hypothesis.")
    print(cur_p)
    # precision
    diffs_prec   = prec - prec_p
    cur_p        = agk.density_p(diffs_prec,0)
    disp("Probability that precision has improved under assumption of 0-hypothesis.")
    print(cur_p)
    # balanced accuracy
    diffs_bacc   = bacc - bacc_p
    cur_p        = agk.density_p(diffs_bacc,0)
    disp("Probability that balanced accuracy has improved under assumption of 0-hypothesis.")
    print(cur_p)
    # all performance measures combined
    diffs_ospec  = diffs_auc + diffs_accs + diffs_sens + diffs_spec + diffs_osens + diffs_ospec + diffs_prec + diffs_bacc
    cur_p        = agk.density_p(diffs_ospec,0)
    disp("Probability that combined performance measure improvements have occured under assumption of 0-hypothesis.")
    print(cur_p)
    
    # ROC curve (get all the coordinates)
    # Ha
    specificities  = t(as.matrix(rocl[[1]]$specificities))
    sensitivities  = t(as.matrix(rocl[[1]]$sensitivities))
    
    # control
    specificitiesl = t(as.matrix(roc_pl[[1]]$specificities))
    sensitivitiesl = t(as.matrix(roc_pl[[1]]$sensitivities))
    
    # get all the ROC curves
    for (ii in 2:bootstr_length) {
      if (length(rocl[[ii]]$specificities) != length(specificities[ii-1,])) {
        stop('length not same!')
      }
      
      if (length(roc_pl[[ii]]$specificities) != length(specificitiesl[ii-1,])) {
        stop('length not same!')
      }
      
      specificities = rbind(specificities,t(as.matrix(rocl[[ii]]$specificities)))
      sensitivities = rbind(sensitivities,t(as.matrix(rocl[[ii]]$sensitivities)))
      
      specificitiesl = rbind(specificitiesl,t(as.matrix(roc_pl[[ii]]$specificities)))
      sensitivitiesl = rbind(sensitivitiesl,t(as.matrix(roc_pl[[ii]]$sensitivities)))
    }
    
    # get CI
    specificities  = as.data.frame(specificities)
    sensitivities  = as.data.frame(sensitivities)
    specificitiesl = as.data.frame(specificitiesl)
    sensitivitiesl = as.data.frame(sensitivitiesl)

    # not ci of mean but mean and percentiles over CV rounds
    agk.mean.quantile.c = cmpfun(agk.mean.quantile)
    spec_mci  = lapply(specificities,FUN = agk.mean.quantile.c,lower = 0.025,upper=0.975)
    sens_mci  = lapply(sensitivities,FUN = agk.mean.quantile.c,lower = 0.025,upper=0.975)
    specl_mci = lapply(specificitiesl,FUN = agk.mean.quantile.c,lower = 0.025,upper=0.975)
    sensl_mci = lapply(sensitivitiesl,FUN = agk.mean.quantile.c,lower = 0.025,upper=0.975)
    
    spec_mci     = as.data.frame(matrix(unlist(spec_mci),ncol = 3,byrow = T))
    sens_mci     = as.data.frame(matrix(unlist(sens_mci),ncol = 3,byrow = T))
    specl_mci    = as.data.frame(matrix(unlist(specl_mci),ncol = 3,byrow = T))
    sensl_mci    = as.data.frame(matrix(unlist(sensl_mci),ncol = 3,byrow = T))
    
    # plot
    # real
    plot(spec_mci$V1[order(spec_mci$V1)],sens_mci$V1[order(spec_mci$V1)],xlim = c(1,0),type='l',lty=1,
         xlab = 'specificity', ylab = '\r\nsensitivity',col='blue',lwd=4, font.lab=2,cex.lab = 2)
    
    lines(spec_mci$V2[order(spec_mci$V1)],sens_mci$V2[order(spec_mci$V1)],xlim = c(1,0),col='blue',lty=1)
    lines(spec_mci$V3[order(spec_mci$V1)],sens_mci$V3[order(spec_mci$V1)],xlim = c(1,0),col='blue',lty=1)
    
    # control
    lines(specl_mci$V1[order(specl_mci$V1)],sensl_mci$V1[order(specl_mci$V1)],xlim = c(1,0),type='l',lty=2,lwd=4,col='red')
    lines(specl_mci$V2[order(spec_mci$V1)],sensl_mci$V2[order(spec_mci$V1)],xlim = c(1,0),col='red',lty=2)
    lines(specl_mci$V3[order(spec_mci$V1)],sensl_mci$V3[order(spec_mci$V1)],xlim = c(1,0),col='red',lty=2)
    
    abline(a=1,b=-1,lty=4,lwd=2)
    title('Receiver-Operating-Curve for Behavioral Classifier',cex.main = 1.5)
    legend(1.02, 1.02, c("ROC of classifier with 95% bounds", "ROC of control classifier with 95% bounds", "hypothetical null"),
           col = c('blue', 'red', 'black'),
           text.col = "black", lty = c(1, 2, 2),
           merge = TRUE, bg = "gray90")
    
    par(mar=c(5,5,4.1,2.1))
    
    # density plots with two densities
    # plots also the density of the performance of the classifier
    Ha_auc               = auc
    Ha_auc               = rep_len(Ha_auc,length.out = length(auc_p))
    cur_dat_gl           = data.frame(null_classifier = auc_p,full_classifier = Ha_auc,classifier = 'elastic net')
    cur_dat              = rbind(cur_dat_gl) #rbind(cur_dat_be,cur_dat_gl,cur_dat_sv)
    cur_dat              = melt(cur_dat,id.vars = c('classifier'))
    cur_dat$AUC_ROC      = cur_dat$value
    cur_dat$value        = NULL
    
    # plot
    p = ggplot(cur_dat,aes(x=AUC_ROC, fill=variable)) + geom_density(alpha=0.25)
    p = p + facet_grid(classifier ~ .) + ggtitle('AUC densities for elastic net classifier compared to null-classifier')
    p = p + geom_vline(aes(xintercept = mean(auc)),colour = 'green',size= 1.5)
    p = p + theme_bw()
    p = p + theme(axis.text=element_text(size=14, face = "bold"),
                  axis.title=element_text(size=20,face="bold"))
    print(p)
  }
  
  if (!do_report_feat_only) {
    # winning model frequencies
    disp('These models have been chosen:')
    cur_mod_sel_nooCV = cur_mod_sel_nooCV[1:bootstr_length]
    cur_tab = table(cur_mod_sel_nooCV)
    print(cur_tab)
    
    # barplot the winning models
    # prep data frame
    cur_mod_sel_nooCV_freq  = as.data.frame(table(cur_mod_sel_nooCV))
    
    # get complexity score
    complexity = data.frame(names(fm),unlist(lapply(fm,length)))
    names(complexity) = c('modname','complexity')
    
    all_models_num          = as.character(1:length(fm))
    all_models_str          = names(fm)
    # add models that have 0 freq
    cur_mod_sel_nooCV_freq$cur_mod_sel_nooCV = as.character(cur_mod_sel_nooCV_freq$cur_mod_sel_nooCV)
    for (ii in 1:length(all_models_str)) {
      if (!all_models_str[ii] %in% cur_mod_sel_nooCV_freq$cur_mod_sel_nooCV) {
        cur_mod_sel_nooCV_freq = rbind(cur_mod_sel_nooCV_freq,c(all_models_str[ii],0))
      }
    }
    
    # transform to numeric
    cur_mod_sel_nooCV_freq$Freq = as.numeric(cur_mod_sel_nooCV_freq$Freq)
    
    # get complexity score
    cur_mod_sel_nooCV_freq$complexity = as.numeric(agk.recode.c(cur_mod_sel_nooCV_freq$cur_mod_sel_nooCV,complexity$modname,complexity$complexity))
    
    # sort
    cur_mod_sel_nooCV_freq = cur_mod_sel_nooCV_freq[order(cur_mod_sel_nooCV_freq$complexity),]
    # get the model names
    names_models_orig               = names(fm)
    cur_mod_sel_nooCV_freq$modnames = cur_mod_sel_nooCV_freq$cur_mod_sel_nooCV
    #cur_mod_sel_nooCV_freq$modnames = names_models_orig[as.numeric(cur_mod_sel_nooCV_freq$cur_mod_sel_nooCV)]
    #cur_mod_sel_nooCV_freq$modnames = agk.recode.c(cur_mod_sel_nooCV_freq$modnames,names_models_orig,names_models)
    #cur_mod_sel_nooCV_freq$modnames = factor(cur_mod_sel_nooCV_freq$modnames, levels = cur_mod_sel_nooCV_freq$modnames)
    cur_mod_sel_nooCV_freq$Freq     = as.numeric(cur_mod_sel_nooCV_freq$Freq)
    
    # fixing the order by complexity, which we prepped above
    cur_mod_sel_nooCV_freq$modnames = factor(cur_mod_sel_nooCV_freq$modnames, levels = cur_mod_sel_nooCV_freq$modnames)
    
    p = ggplot(data = cur_mod_sel_nooCV_freq, aes(modnames,Freq))
    p = p+geom_bar(stat="identity") + ylab("Frequency model selected") + xlab('models ordered by model complexity')
    p = p + coord_cartesian(ylim=c(0,700))
    p = p + ggtitle(paste0("Frequency of model selection in ",
                           sum(cur_mod_sel_nooCV_freq$Freq)," cross validation rounds"))
    p = p + theme(text = element_text(size=20),
                  axis.text.x = element_text(angle=90, hjust=1)) 
    p = p + coord_flip()
    print(p)
    
    # get the most often chosen model
    if (!do_report_feat_only) {
      if (length(which(max(cur_tab) == cur_tab)) > 1) {
        stop('At least two models won equally often.')
      }
      winning_mod = names(cur_tab)[which(max(cur_tab) == cur_tab)]
    }
  }
  
  if (do_report_feat_only) {
    # winning model beta values with CIs
    win_mods_distr_c = list_winning_model_c_nooCV_op
    win_mods_distr_l = list_winning_model_l_nooCV_op
  } else {
    # winning model beta values with CIs
    list_winning_model_c_nooCV = list_winning_model_c_nooCV[1:bootstr_length]
    list_winning_model_l_nooCV = list_winning_model_l_nooCV[1:bootstr_length]
    win_mods_distr_c = list_winning_model_c_nooCV[which(winning_mod == cur_mod_sel_nooCV)]
    win_mods_distr_l = list_winning_model_l_nooCV[which(winning_mod == cur_mod_sel_nooCV)]
  }
  
  # make a data frame of it:
  win_mod_coefs        = as.matrix(win_mods_distr_c[[1]])
  win_mod_coefs        = as.data.frame(t(win_mod_coefs))
  names(win_mod_coefs) = win_mods_distr_l[[1]]
  
  for (ii in 2:length(win_mods_distr_c)) {
    cur_win_mod_coefs        = as.matrix(t(win_mods_distr_c[[ii]]))
    cur_win_mod_coefs        = as.data.frame(cur_win_mod_coefs)
    names(cur_win_mod_coefs) = win_mods_distr_l[[ii]]
    win_mod_coefs = rbind.fill(win_mod_coefs,cur_win_mod_coefs)
  }
  imp_0 = function(x) {x[is.na(x)] = 0; return(x)}
  win_mod_coefs = as.data.frame(lapply(win_mod_coefs,FUN=imp_0))
  
  # now we get the mean, upper and lower bootstrapped CI
  # ci_res = lapply(win_mod_coefs,FUN=agk.boot.ci,cur_fun=mean,
  #                 lower=0.025,upper=0.975,R=1000)
  
  ci_res = lapply(win_mod_coefs,FUN=agk.mean.quantile.c,
                  lower=0.025,upper=0.975)
  ci_res = as.data.frame(t(as.data.frame(ci_res)))
  names(ci_res) = c('mean','lower','upper')
  
  # FIRST INTERCEPT IS THE GROUP PRED INTERCEPT, SECOND IS FROM THE BEHAV MODEL
  ci_res$coef                                = row.names(ci_res)
  ci_res$coef[ci_res$coef == 'Intercept']    = 'Int_behav_model'
  ci_res$coef[ci_res$coef == 'X.Intercept.'] = 'Int_classifier'
  
  labels_sources = c("catgambling","catnegative","catpositive","Int_behav_model","Int_classifier","smoking_ftdt") 
  ci_res$coef    = factor(ci_res$coef)
  labels_betas   = agk.recode(levels(ci_res$coef),labels_sources,
                              c("gambling cues","negative cues","positive cues","intercept behavioral model","intercept of classifier","smoking severity"))
  ci_res$coef    = factor(ci_res$coef, levels = levels(ci_res$coef), labels = labels_betas)
  
  if (which_study == 'MRT') {
    # names simpler:
    ci_res$coef = gsub('SS__grp01_noCov_','',ci_res$coef)
    ci_res$coef = gsub('SS__','',ci_res$coef)
    ci_res$coef = gsub('PicGamOnxAcc','PIT',ci_res$coef)
    ci_res$coef = gsub('PicGamOnxacc','PIT',ci_res$coef)
    ci_res$coef = gsub('ROI_LR_','',ci_res$coef)
    ci_res$coef = gsub('_LR','',ci_res$coef)
    ci_res$coef = gsub('noCov_PPI_','',ci_res$coef)
    ci_res$coef = gsub('X','x',ci_res$coef)
  }
  
  p = ggplot(data = ci_res, aes(coef,mean))
  p = p+geom_bar(stat="identity")
  p = p + geom_errorbar(aes(ymin=lower,ymax=upper), size=1.3, color=cbbPalette[4],
                        width = 0) + ylab("mean (95% CI)\n")
  
  p <- p + ggtitle("Estimated regression weights of winning model")
  p = p + theme(text = element_text(size=25),
                axis.text.x = element_text(angle=45, hjust=1)) 
  p = p + xlab("regression weights")
  print(p)
  
  if (which_study == 'MRT') {
    # giving a grouped overview
    grouping                                                  = rep(NA,length(ci_res$coef))
    grouping[grep('^Pic',ci_res$coef)]                        = 'cue_reactivity'
    grouping[grep('^PPI_Amy.*OrG',ci_res$coef)]               = 'Amy_to_OFC'
    grouping[grep('^PPI_Amy.*StrAsCaud',ci_res$coef)]         = 'Amy_to_Striatum'
    grouping[grep('^PPI_Amy.*StrAsPut',ci_res$coef)]          = 'Amy_to_Striatum'
    grouping[grep('^PPI_Amy.*Acc',ci_res$coef)]               = 'Amy_to_Striatum'
    grouping[grep('^PIT',ci_res$coef)]                        = 'PIT'
    grouping[grep('^PPI_Acc_',ci_res$coef)]                   = 'Accumbens_to_Str_Amy'
    ci_res$coef[ci_res$coef == 'x.grp_classifier_intercept.'] = 'Classifier_Icpt'
    grouping[grep('Classifier_Icpt',ci_res$coef)]             = 'cue_reactivity'
    
    ci_res$grouping = factor(grouping,levels = c('cue_reactivity','PIT','Amy_to_OFC','Amy_to_Striatum','Accumbens_to_Str_Amy'),
                             labels = c('cue reactivity','PIT','PPI PIT seed Amygdala to OFC','PPI PIT seed Amygdala to Accumbens','PPI PIT seed Accumbens to Amygdala'))
    
    # shorter names even
    ci_res$coef = gsub('Accumbens','Acc',ci_res$coef)
    ci_res$coef = gsub('Amygdala','Amy',ci_res$coef)
    ci_res$coef = gsub('PPI','c',ci_res$coef)
    ci_res$coef = gsub('PITx','PIT',ci_res$coef)
    ci_res$coef = gsub('c_','c',ci_res$coef)
    ci_res$coef = gsub('StrAso','StrAs',ci_res$coef)
    ci_res$coef = gsub('StrAs','',ci_res$coef)
    ci_res$coef = gsub('ROI_','',ci_res$coef)
    ci_res$coef = gsub('Picgam_','gam_',ci_res$coef)
    ci_res$coef = gsub('Picneg_','neg_',ci_res$coef)
    ci_res$coef = gsub('Picpos_','pos_',ci_res$coef)
    ci_res$coef = gsub('^PITpos_','pos_',ci_res$coef)
    ci_res$coef = gsub('^PITneg_','neg_',ci_res$coef)
    ci_res$coef = gsub('^PITgam_','gam_',ci_res$coef)
    ci_res$coef = gsub('^cAmy','Amy',ci_res$coef)
    ci_res$coef = gsub('^cAcc','Acc',ci_res$coef)
    ci_res$coef = gsub('^AccPITgam','Acc_PITgam',ci_res$coef)
    ci_res$coef = gsub('^AccPITneg','Acc_PITneg',ci_res$coef)
    ci_res$coef = gsub('^AccPITpos','Acc_PITpos',ci_res$coef)
    
    
    # plotting grouped
    message('displaying the set grouped')
    #ci_res_red = ci_res[abs(ci_res$mean) > 1,]
    ci_res_red = ci_res
    
    p = ggplot(data = ci_res_red, aes(coef,mean))
    p = p+geom_bar(stat="identity")
    p = p + geom_errorbar(aes(ymin=lower,ymax=upper), size=1.3, color=cbbPalette[4],
                          width = 0) + ylab("mean (95% CI over CV rounds)\n\n\n")
    
    p <- p + ggtitle("Estimated regression weights with CIs") + theme_bw()
    p = p + theme(text = element_text(size=15),
                  axis.text.x = element_text(angle=45, hjust=1,face='bold')) 
    p = p + facet_wrap(~ grouping, scales = "free_x",nrow=3,ncol=2)
    print(p)
    
    # check how many betas sig
    ci_res_upper_lower = as.matrix(ci_res[c('upper','lower')])
    cur_fun = function(x) {
      return(sign(x[1]) == sign(x[2])) 
    }
    sum(apply(ci_res_upper_lower,MARGIN = 1,FUN = cur_fun))
  }
}


