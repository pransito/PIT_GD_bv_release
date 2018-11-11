# PREAMBLE ====================================================================
# v6.0; stripped down to glmnet
# k-fold CV of the whole procedure
# script which checks for CV score of the whole procedure
# needs before: import_data_pdt.R; select_study.R; group_pred_init.R;

# PREPARATION =================================================================
# get the study selected dat_match which contains demographic info
# group label, smoking, severity information
dat_match = dat_match_bcp_study_selected

if (initial_scale) {
  dat_match[[pred_to_control]] = scale(dat_match[[pred_to_control]])
}

# getting the exp model params
featmod_coefs = featmod_coefs_bcp

# getting the additional features
feature_clusters = feature_clusters_bcp

# getting the folds
# cur_k for k-fold CV (k = n is LOOCV)
cur_lev = unique(dat_match$VPPG)
if (isempty(outer_k)) {
  # LOOCV
  cur_k   = length(cur_lev)
} else {
  # k-fold
  cur_k   = outer_k
}

# when predicting group we have a function getting stratified folds
if (CV == 'wio') {
  flds = f_K_fold_b2g(cur_data = dat_match,k = cur_k,group = "HCPG")
} else if (CV == 'noo') {
  flds = list()
  cur_trte        = list()
  cur_trte$train  = 1:length(cur_lev)
  cur_trte$test   = 1:length(cur_lev)
  flds[[1]]       = cur_trte
} else {
  stop('Variable CV as unexpected value.')
}

# unit test: check if flds is really balanced in test set
if (unit_test_strat_outCV) {
  for (zz in 1:length(flds)) {
    check_half = table(dat_match$HCPG[flds[[zz]]$test])
    if(check_half[1] != check_half[2]) {
      stop('Unbalanced test set for outer CV!')
    }
  }
}

# unit test: check if flds is really balanced in test set for exp. mod
# here one model holds for all
if (unit_test_strat_outCV) {
  for (zz in 1:length(flds)) {
    check_half = table(featmod_coefs[[1]]$HCPG[flds[[zz]]$test])
    if(check_half[1] != check_half[2]) {
      stop('Unbalanced test set for outer CV!')
    }
  }
}

# unit test: check if flds is really balanced in test set for add feat
# here one set holds for all
if (unit_test_strat_outCV) {
  for (zz in 1:length(flds)) {
    check_half = table(featmod_coefs[[1]]$HCPG[flds[[zz]]$test])
    if(check_half[1] != check_half[2]) {
      stop('Unbalanced test set for outer CV!')
    }
  }
}

# MODEL SELECTION =============================================================
total = length(flds)

# get training data for each fold
featmod_coefs_cv     = list()
# if scaling then note down training data's mean and sd
featmod_coefs_cv_msd = list ()
if (scale_for_CV) {
  disp('Scaling training data (predictors).')
} else {
  disp('Not scaling training data (predictors).')
}

for (kk in 1:total) {
  # get the complete feature clusters (all subs)
  tmp        = featmod_coefs
  
  # get the current train indices
  cur_traini     = flds[[kk]]$train
  cur_train_subs = cur_lev[cur_traini]
  tmp_msd        = list()
  # select only trainings subs
  for (ii in 1:length(tmp)) {
    cur_ind   = which(row.names(tmp[[ii]]) %in% cur_train_subs)
    tmp[[ii]] = tmp[[ii]][cur_ind,]
    
    # scale the training data (NEW!)
    if (scale_for_CV) {
      cur_training_dat                                  = scale(tmp[[ii]][grep('HCPG',names(tmp[[ii]]),invert=T)])
      tmp[[ii]][grep('HCPG',names(tmp[[ii]]),invert=T)] = cur_training_dat
      tmp_msd[[ii]]                                     = attributes(cur_training_dat)
    }
  }
  featmod_coefs_cv[[kk]]            = tmp
  featmod_coefs_cv_msd[[kk]]        = tmp_msd
  names(featmod_coefs_cv_msd[[kk]]) = names(featmod_coefs_cv[[kk]])
}

# make the foldid vector for all innerCV folds
cur_fold_id_list = list()
for (kk in 1:total) {
  # make the foldid vector
  if (strat_innerCV) {
    cur_flds    = f_K_fold_b2g(featmod_coefs_cv[[kk]][[1]],'HCPG',cnfolds)
  } else {
    cur_flds    = f_K_fold(length(featmod_coefs_cv[[kk]][[1]]$HCPG),cnfolds)
  }
  cur_fold_id = agk.get.foldid(cur_flds)
  
  # pack
  cur_fold_id_list[[kk]] = cur_fold_id
}

if (use_behav_params & c_mod == F) {
  # perform the feature selection in each fold
  disp('Running model selection in training folds...')
}

# function body
cur.mod.selection.fun = function(kk) {
  # predicting using each cluster of predictors
  # model (params set) selection using ridge (alphas = c(0))
  full_mods_ms = list()
  cvme_vec     = c()
  
  if (!put_all_behav_vars_in) {
    # multiple times doing model comparison
    ms_reps_v = c()
    for (rr in 1:ms_reps) {
      # getting the folds needed
      cur_fold_id = cur_fold_id_list[[kk]]
      
      # get a new folding on every model selection rep round
      # the follwoing two lines were missing!
      cur_flds    = f_K_fold_b2g(featmod_coefs_cv[[kk]][[1]],'HCPG',cnfolds)
      cur_fold_id = agk.get.foldid(cur_flds)
      
      for (ii in 1:length(featmod_coefs_cv[[kk]])) {
        
        # unit test that stratified folds for inner CV
        if (unit_test_strat_innCV) {
          for (zz in unique(cur_fold_id)[order(unique(cur_fold_id))]) {
            check_half = table(featmod_coefs_cv[[kk]][[ii]][['HCPG']][cur_fold_id == zz])
            if(check_half[1] != check_half[2]) {
              stop('Unbalanced test set for inner CV!')
            }
          }
        }
        
        cur_data = featmod_coefs_cv[[kk]][[ii]]
        
        # unit test check that data is in correct order
        cur_sub_grp_matching = subset(sub_grp_matching, row.names(sub_grp_matching) %in% row.names(cur_data))
        cur_data_pm          = cur_data[c('HCPG')]
        cur_data_pm          = cur_data_pm[order(row.names(cur_data_pm)),]
        cur_sub_grp_matching = cur_sub_grp_matching[order(row.names(cur_sub_grp_matching)),]
        stopifnot(all(row.names(cur_data_pm) == row.names(cur_sub_grp_matching)))
        stopifnot(all(cur_data_pm == cur_sub_grp_matching))
        
        # check if a single-predictor case (a model)
        if (length(cur_data) == 2) {
          # adding some random noise variable
          cur_data$err = randn(length(cur_data[,1]),1)
        } else {
          # cur data can stay as it is
        }
        
        # make formula
        cur_formula = as.formula('HCPG ~ .')
        
        out = tryCatch(
          # work around for bug
          # https://github.com/lmweber/glmnet-error-example/blob/master/glmnet_error_example.R 
          {
            cur_lam            = NULL
            cur_cvmoda         = glmnetUtils::cva.glmnet(cur_formula,family=cur_family,
                                                         data=cur_data,type.measure=type_measure,
                                                         alpha = 0,foldid = cur_fold_id,
                                                         use.model.frame=useModelFrame,grouped=doGrouped,lambda = cur_lam,
                                                         standardize = des_stand_glmnet)
          },
          error=function(cond) {
            message("Had to resort to predefined lambdas")
            cur_lam = des_lambdas
            cur_cvmoda         = glmnetUtils::cva.glmnet(cur_formula,family=cur_family,
                                                         data=cur_data,type.measure=type_measure,
                                                         alpha = 0,foldid = cur_fold_id,lambda=cur_lam,
                                                         use.model.frame=useModelFrame,grouped=doGrouped,
                                                         standardize = des_stand_glmnet)
          }
        )
        cur_cvmoda = out
        
        full_mods_ms[[ii]] = agk.glmnetUtils.cvalpha.getwinner(cur_cvmoda)
        
        # unit test: predictors before and after need to align
        cur_coef = coef(full_mods_ms[[ii]]$winning_model,s='lambda.min')
        cur_coef = row.names(cur_coef)
        cur_coef = cur_coef[grep('(Intercept)',cur_coef,fixed=T,invert = T)]
        cur_fprd = names(cur_data)[grep('HCPG',names(cur_data),invert = T)]
        if (!all(cur_coef == cur_fprd)) {
          stop('Behav model selection: formula predictors and predictors in resulting model do not align!')
        }
        cvme_vec[ii]       = full_mods_ms[[ii]]$cvme
      }
      
      # picking the best model and extracting
      # the survived features (in case lasso is used)
      cur_mod_sel   = which(min(cvme_vec) == cvme_vec)
      
      if (length(cur_mod_sel) > 1) {
        warning('We have a tie in performance for exp params extraction. Using simpler model.')
        cur_mod_sel = cur_mod_sel[1]
      }
      ms_reps_v[rr] = cur_mod_sel
    }
    
    # take the modus cur_mod_sel, if tie, it takes lower complexity model
    cur_mod_sel   = as.numeric(names(sort(table(ms_reps_v),decreasing=TRUE)[1]))
    cur_mod_sel_n = names(featmod_coefs_cv[[kk]])[cur_mod_sel]
  } else {
    cur_mod_sel   = 'kitchen_sink'
    cur_mod_sel_n = 9999 
  }
  
  # packing the winning model
  featmod_coefs_sel = list()
  
  if (put_all_behav_vars_in) {
    # this is experimental to check what happens if we do not do model selection
    tmp = agk.merge.dfs(featmod_coefs_cv[[kk]],'HCPG', dat_match)
    # survived coefs (taking all even if some are 0; it is about the whole model)
    cur_surv = names(tmp)[grep('HCPG',names(tmp),invert = T)]
  } else {
    # extracting the survived features
    # cur featmod coefs
    tmp = featmod_coefs_cv[[kk]][[cur_mod_sel_n]]
    # survived coefs
    # survived coefs (taking all even if some are 0; it is about the whole model)
    cur_surv = names(tmp)[grep('HCPG',names(tmp),invert = T)]
  }
  
  # pack the winning model
  featmod_coefs_sel[[1]] = tmp

  # from selected features
  for (ii in 1:length(featmod_coefs_sel)) {
    if (ii == 1) {
      dat_sel_feat = featmod_coefs_sel[[1]]
    } else {
      dat_sel_feat = agk.merge.df.by.row.names(dat_sel_feat,featmod_coefs_sel[[ii]])
    }
  }
  
  # packing the result of this cv fold
  cur_featmod_coefs_sel_cv = list()
  cur_featmod_coefs_sel_cv$cur_featmod_coefs_sel_cv = dat_sel_feat
  #cur_featmod_coefs_sel_cv$cur_mod_sel              = cur_mod_sel # should only work with cur_mod_sel_n
  cur_featmod_coefs_sel_cv$cur_mod_sel_n            = cur_mod_sel_n
  cur_featmod_coefs_sel_cv$fit_mod_sel              = full_mods_ms[[cur_mod_sel]]
  # return
  return(cur_featmod_coefs_sel_cv)
}

# compile the function
cur.mod.selection.fun.c = cmpfun(cur.mod.selection.fun)

if (use_behav_params & c_mod == F) {
  if (CV == 'wio') {
    # run it using multiple processes
    total = length(flds)
    cl    = parallel::makeCluster((detectCores()-1))
    registerDoSNOW(cl)
    featmod_coefs_sel_cv = foreach(kk=1:total, .packages=c('glmnetUtils','cvTools','boot','pROC','pracma'),.verbose=T,.export = c()) %dopar% {
      cur.mod.selection.fun.c(kk)
    } 
    stopCluster(cl)
  } else {
    featmod_coefs_sel_cv      = list()
    featmod_coefs_sel_cv[[1]] = cur.mod.selection.fun.c(1)
  }
  
  # unpack the featmod_results and make a vector which tells me which model was selected
  featmod_coefs_sel_cv_list = list()
  #cur_mod_sel_vec           = c()
  cur_mod_sel_n_vec         = c()
  scaled_train_convar       = list()
  cur_mod_sel_fit           = list() # good to make directly predictions from it
  for (kk in 1:length(featmod_coefs_sel_cv)) {
    #cur_mod_sel_vec[kk]             = featmod_coefs_sel_cv[[kk]]$cur_mod_sel
    cur_mod_sel_n_vec[kk]           = featmod_coefs_sel_cv[[kk]]$cur_mod_sel_n
    cur_mod_sel_fit[[kk]]           = featmod_coefs_sel_cv[[kk]]$fit_mod_sel
    featmod_coefs_sel_cv_list[[kk]] = featmod_coefs_sel_cv[[kk]]$cur_featmod_coefs_sel_cv
  }
  featmod_coefs_sel_cv = featmod_coefs_sel_cv_list
}

# FEATURE COLLECTION AND SCALING ==============================================
# getting all the to-be innerCV's feature clusters ready
total                   = length(flds)
# get training data for each fold
feature_clusters_cv     = list()
# saving scaling info as well
feature_clusters_cv_msd = list()
if (scale_for_CV) {
  disp('Scaling training data of additional features (predictors)')
} else {
  disp('Not scaling training data of additional features (predictors)')
}
for (kk in 1:total) {
  # get the complete feature clusters (all subs)
  tmp        = feature_clusters
  # prepare saving scaling info
  tmp_msd    = list()
  
  # get the current train indices
  cur_traini     = flds[[kk]]$train
  cur_train_subs = cur_lev[cur_traini]
  if (add_cr_pp | add_cr_ra) {
    # select only trainings subs in additional features
    # parameters of training have already been selected for training sample
    for (ii in 2:length(tmp)) {
      cur_ind   = which(row.names(tmp[[ii]]) %in% cur_train_subs)
      tmp[[ii]] = tmp[[ii]][cur_ind,]
      
      # scale the training data (and noting down scaling info for applying later to test data)
      if (scale_for_CV) {
        cur_df                                            = tmp[[ii]][grep('HCPG',names(tmp[[ii]]),invert=T)]
        cur_df                                            = scale(cur_df)
        cur_training_dat                                  = cur_df
        tmp[[ii]][grep('HCPG',names(tmp[[ii]]),invert=T)] = cur_training_dat
        tmp_msd[[ii]]                                     = attributes(cur_training_dat)
      }
    }
  }
  if (use_behav_params & c_mod == F) {
    # REPLACE the first feature cluster by innerCV-selected model coefs
    # which already are only the respective subs of training fold
    # only if behav mods been used and model selection has been done
    tmp[[1]] = featmod_coefs_sel_cv[[kk]]
  }
  feature_clusters_cv[[kk]]     = tmp
  feature_clusters_cv_msd[[kk]] = tmp_msd
}

feature_clusters_sel_cv = feature_clusters_cv

# PREP FOR BUILDING THE COMPLETE MODEL ========================================
if (c_mod == F) {
  feature_clusters_sel_cv_noNULL = list()
  for (kk in 1:length(feature_clusters_sel_cv)) {
    # merging selected features
    # first get rid of NULL
    feature_clusters_sel_cv_noNULL[[kk]] = feature_clusters_sel_cv[[kk]][unlist(lapply(feature_clusters_sel_cv[[kk]],is.null)) == FALSE]
    for (ii in 1:length(feature_clusters_sel_cv_noNULL[[kk]])) {
      if (ii == 1) {
        dat_sel_feat = feature_clusters_sel_cv_noNULL[[kk]][[1]]
      } else {
        dat_sel_feat = agk.merge.df.by.row.names(dat_sel_feat,feature_clusters_sel_cv_noNULL[[kk]][[ii]],x1 = 'HCPG',x2 = 'HCPG')
      }
    }
    feature_clusters_sel_cv[[kk]] = dat_sel_feat
  }
}

# BUILDING THE COMPLETE MODEL =================================================
# (using selected model features and additional features features to predict)
if (c_mod == F) {
  disp('Using selected model (and additional) features to predict in every training fold')
} else if (c_mod == T) {
  disp('Running the control model.')
}

# function body
cur.compl.model.fun = function(kk) {
  # function to do the complete model CV
  # is a function of the environment
  # list2env(el, .GlobalEnv)
  
  # the current selected features (so only training data)
  if(c_mod == F) {
    dat_sel_feat = feature_clusters_sel_cv[[kk]]
    f = function(x) {return(any(is.na(x)))}
    if(any(unlist(lapply(dat_sel_feat,FUN = f)))) {
      stop('NAs in dat_sel_feat!!')
    }
  } else if (c_mod == T) {
    dat_sel_feat = featmod_coefs_cv[[kk]][[1]]
    dat_sel_feat = dat_sel_feat[c('HCPG')]
  } else {
    stop('dat_sel_feat is empty!')
  }
  
  # get the complete data frame with all selected features
  tmp = dat_sel_feat
  
  if (do_con) {
    # add control variables
    cdm       = dat_match[c('VPPG',pred_to_control)]

    # new_names = paste0('pred_',names(cdm)[grep(x=names(cdm),pattern = 'VPPG',invert = T)])
    # # change to pred_ variable
    # names(cdm)[grep(x=names(cdm),pattern = 'VPPG',invert = T)] = new_names
    # ad the con variables
    tmp$subject             = row.names(tmp)
    tmp                     = merge(tmp,cdm,by.x = 'subject',by.y = 'VPPG')
    dat_sel_feat            = tmp
    row.names(dat_sel_feat) = dat_sel_feat$subject
    tmp$subject             = NULL
    dat_sel_feat$subject    = NULL
    
    # scaling
    if (scale_for_CV) {
      cur_scaled_train_convar                                 = scale(dat_sel_feat[grep('HCPG',names(dat_sel_feat),invert=T)])
      dat_sel_feat[grep('HCPG',names(dat_sel_feat),invert=T)] = cur_scaled_train_convar
    } else {
      cur_scaled_train_convar = list()
    }
    
  }
  
  # make the complete formula
  form_cmpl = as.formula('HCPG ~ .')
  
  # getting the folds needed
  cur_fold_id = cur_fold_id_list[[kk]]
  
  # unit test that stratified folds for inner CV
  if (unit_test_strat_innCV) {
    for (zz in unique(cur_fold_id)[order(unique(cur_fold_id))]) {
      check_half = table(dat_sel_feat$HCPG[cur_fold_id == zz])
      if(check_half[1] != check_half[2]) {
        stop('Unbalanced test set for inner CV!')
      }
    }
  }
  
  # unit test check that data is in correct order
  cur_data             = dat_sel_feat
  cur_data$subject     = gsub(pattern = '\\.[1-9]','',row.names(cur_data))
  cur_data             = cur_data[!duplicated(cur_data$subject),]
  cur_sub_grp_matching = subset(sub_grp_matching, row.names(sub_grp_matching) %in% cur_data$subject)
  cur_sub_grp_matching$subject = row.names(cur_sub_grp_matching)
  cur_data             = cur_data[c('HCPG','subject')]
  cur_data             = cur_data[order(cur_data$subject),]
  cur_sub_grp_matching = cur_sub_grp_matching[order(row.names(cur_sub_grp_matching)),]
  stopifnot(all(cur_data$subject == cur_sub_grp_matching$subject))
  stopifnot(all(cur_data == cur_sub_grp_matching))
  
  if (which_ML == 'ML') {
    # adjust prediction model using an ML (glmnet) technique
    # alpha is free here and will be cv'd
    if (length(dat_sel_feat) > 2) {
      
      all_winning_models = list()
      for (ll in 1:fullm_reps) {
        print(ll)
        cur_flds      = f_K_fold_b2g(dat_sel_feat,'HCPG',cnfolds)
        cur_fold_id   = agk.get.foldid(cur_flds)
        
        # learning
        cur_cvmoda         = glmnetUtils::cva.glmnet(form_cmpl,family=cur_family,
                                                     data=dat_sel_feat,type.measure=type_measure,
                                                     alpha = alphas,foldid = cur_fold_id,use.model.frame=useModelFrame,
                                                     grouped=doGrouped,standardize = des_stand_glmnet)
        
        winning_model = agk.glmnetUtils.cvalpha.getwinner(cur_cvmoda)
        all_winning_models[[ll]] = winning_model
      }
      
      cur_full_mod_cv   = agk.get.mean.model(all_winning_models)
      
      # unit test: predictors before and after need to align
      #cur_coef = coef(cur_full_mod_cv$winning_model,s='lambda.min')
      cur_coef = names(cur_full_mod_cv)
      cur_coef = cur_coef[grep('(Intercept)',cur_coef,invert = T,fixed = T)]
      cur_fprd = as.character(names(dat_sel_feat))
      if (!all(cur_coef %in% cur_fprd)) {
        stop('Complete model: formula predictors and predictors in resulting model do not align!')
      }
    } else if (length(dat_sel_feat) == 1) {
      # no model if there are no predictors
      cur_full_mod_cv = 'no_model'
    } else if (length(dat_sel_feat) == 2) {
      # default to non-ML method
      cur_full_mod_cv = robust::glmRob(form_cmpl,family=cur_family,data=dat_sel_feat)
    } else {
      stop('dat_sel_feat seems to be empty')
    }
  } else if (which_ML == 'SVM') {
    cur_cost        = c(0.001,0.01,2)
    cur_gamm        = 2^(-1:1)
    cur_tune        = tune.svm(form_cmpl,data=dat_sel_feat,
                               cost = cur_cost,gamma= cur_gamm,
                               tunecontrol = tune.control(nrepeat = fullm_reps))
    best_svm        = cur_tune$best.model
    cur_full_mod_cv = best_svm
  }
  
  # return
  if (do_con) {
    return(list(cur_full_mod_cv,cur_scaled_train_convar))
  } else {
    return(list(cur_full_mod_cv))
  }
}

# compile the fun
cur.compl.model.fun.c = cmpfun(cur.compl.model.fun)

# run the complete model CV using parallel processing
if (all_alphas == F | (all_alphas == T & c_mod == T) | (use_behav_params == F & all_alphas == F)) {
  
  if (CV == 'wio') {
    total = length(flds)
    cl    = parallel::makeCluster(detectCores()-1)
    registerDoSNOW(cl)
    full_mod_cv = foreach(kk=1:total, .packages=c('glmnetUtils','pracma','robustbase','e1071'),
                          .verbose=T,.export = c()) %dopar% {
                            cur.compl.model.fun.c(kk)
                          } 
    stopCluster(cl)
  } else {
    full_mod_cv      = list()
    full_mod_cv[[1]] = cur.compl.model.fun.c(1)
  }
  
  # unpack
  actual_full_mod_cv         = list()
  scaled_train_convar        = list()
  for (ll in 1:length(full_mod_cv)) {
    actual_full_mod_cv[[ll]]  = full_mod_cv[[ll]][[1]]
    if (do_con) {
      scaled_train_convar[[ll]] = full_mod_cv[[ll]][[2]]
    }
  }
  full_mod_cv = actual_full_mod_cv
} else {
  # in case we have already done model selection and fitting of complete model in one step
  full_mod_cv = cur_mod_sel_fit
}

# GET THE TEST DATA ===========================================================
# prep: get a data frame which holds all additional features from all subs
if (length(feature_clusters) > 1) {
  for (ii in 2:length(feature_clusters)) {
    if (ii == 2) {
      dat_all_add_feat = feature_clusters[[ii]]
    } else {
      stop('length(feature_clusters) > 3 not implemented; need to do merg.by.row.names')
      dat_all_add_feat = merge(dat_all_add_feat,feature_clusters[[ii]],
                               by = c('subject','HCPG'),all.x = T)
    }
  }
}

test_data_coefs     = list()
for (kk in 1:length(flds)) {
  # get the test data parameters
  # get the test subjects(s)
  cur_testi             = flds[[kk]]$test
  cur_test_sub          = cur_lev[cur_testi]
  
  if (use_behav_params & c_mod == F) {
    # get the current chosen params extraction model
    # and the pertinent parameters of the test sub
    cur_mod_sel_n            = cur_mod_sel_n_vec[kk]
    if (!put_all_behav_vars_in) {
      cur_df                   = fm[[cur_mod_sel_n]]
      # select the correct test subjects
      cur_modpe_coef           = subset(cur_df,row.names(cur_df) %in% cur_test_sub)
      cur_modpe_coef           = agk.clean.intercept.name(cur_modpe_coef)
    } else {
      cur_df                   = agk.merge.dfs(fm,'HCPG',dat_match)
      cur_modpe_coef           = subset(cur_df,row.names(cur_df) %in% cur_test_sub)
    }

  }
  
  if (add_cr_pp || add_cr_ra) {
    # getting the additional features coefs
    cur_test_data_coefs   = dat_all_add_feat[row.names(dat_all_add_feat) %in% cur_test_sub,]
  } else if (c_mod == F) {
    cur_test_data_coefs   = cur_modpe_coef
  } else {
    cur_test_data_coefs            = data.frame(subject = cur_test_sub)
    cur_test_data_coefs$HCPG       = agk.recode.c(cur_test_data_coefs$subject,dat_match$VPPG,dat_match$HCPG)
    row.names(cur_test_data_coefs) = cur_test_data_coefs$subject 
  }
  
  if (add_cr_pp || add_cr_ra) {
    cur_chosen_feats      = names(feature_clusters_sel_cv[[kk]]) 
    cur_chosen_feats      = cur_chosen_feats[grep('HCPG',cur_chosen_feats,invert=T)]
    cur_test_data_coefs   = cur_test_data_coefs[colnames(cur_test_data_coefs) %in% cur_chosen_feats]
  }
  
  if ((add_cr_pp || add_cr_ra) & use_behav_params) {
    # merging extracted model params and additional features params
    #cur_test_data_coefs   = agk.merge.df.by.row.names(cur_modpe_coef,cur_test_data_coefs,dat_match)
    cur_test_data_coefs   = agk.merge.df.by.row.names(cur_modpe_coef,cur_test_data_coefs)
  }
  
  if (do_con) {
    # adding the control variable(s)
    cur_con_vars            = dat_match[pred_to_control]
    cur_rn                  = dat_match$VPPG
    cur_con_vars            = subset(cur_con_vars,dat_match$VPPG %in% cur_test_sub)
    cur_rn                  = cur_rn[which(dat_match$VPPG %in% cur_test_sub)]
    cur_con_vars            = data.frame(cur_con_vars)
    row.names(cur_con_vars) = cur_rn
    names(cur_con_vars)     = pred_to_control
    if (c_mod == F) {
      cur_test_data_coefs     = agk.merge.df.by.row.names(cur_test_data_coefs,cur_con_vars)
    } else {
      cur_test_data_coefs     = cur_con_vars
    }
  }
  
  # name correction (changing : to _of_)
  names(cur_test_data_coefs) = gsub(':','_of_',names(cur_test_data_coefs))
  
  if(scale_for_CV) {
    # correcting test data for training data's scaling
    # behav data
    if (c_mod == F & use_behav_params) {
      if (!put_all_behav_vars_in) {
        cur_scaling_info                   = featmod_coefs_cv_msd[[kk]][[cur_mod_sel_n]]
        cur_sd                             = cur_scaling_info$`scaled:scale`
        cur_m                              = cur_scaling_info$`scaled:center`
        cur_sd                             = agk.clean.intercept.name(cur_sd)
        
      } else {
        cur_scaling_info                   = featmod_coefs_cv_msd[[kk]][[1]]
        cur_sd                             = cur_scaling_info$`scaled:scale`
        cur_m                              = cur_scaling_info$`scaled:center`
        names(cur_sd)                      = paste0(names(featmod_coefs_cv_msd[[kk]])[1],'_',names(cur_sd))
        names(cur_m)                      = paste0(names(featmod_coefs_cv_msd[[kk]])[1],'_',names(cur_m))
        for (tt in 2:length(featmod_coefs_cv_msd[[kk]])) {
          cur_scaling_info                   = featmod_coefs_cv_msd[[kk]][[tt]]
          sd_tmp                             = cur_scaling_info$`scaled:scale`
          mn_tmp                             = cur_scaling_info$`scaled:center`
          names(sd_tmp)    = paste0(names(featmod_coefs_cv_msd[[kk]])[tt],'_',names(sd_tmp))
          names(mn_tmp)    = paste0(names(featmod_coefs_cv_msd[[kk]])[tt],'_',names(mn_tmp))
          cur_sd                             = c(cur_sd,sd_tmp)
          cur_m                              = c(cur_m,mn_tmp)
        }
        names(cur_sd) = gsub('.','_of_',names(cur_sd),fixed=T)
        names(cur_m)  = gsub('.','_of_',names(cur_m),fixed=T)
        names(cur_test_data_coefs) = gsub('.','_of_',names(cur_test_data_coefs),fixed=T)
      }
      cur_test_data_coefs[names(cur_sd)] = scale(cur_test_data_coefs[names(cur_sd)],center = cur_m,scale=cur_sd)
    }
    
    # add feat data
    if (c_mod == F & add_cr_pp == T) {
      cur_scaling_info = feature_clusters_cv_msd[[kk]][[2]]
      cur_sd           = cur_scaling_info$`scaled:scale`
      cur_m            = cur_scaling_info$`scaled:center`
      cur_test_data_coefs[names(cur_sd)] = scale(cur_test_data_coefs[names(cur_sd)],center = cur_m,scale=cur_sd)
    }
    
    # convar data
    if (do_con) {
      cur_scaling_info = scaled_train_convar[[kk]]
      cur_scaling_info = attributes(cur_scaling_info)
      cur_sd           = cur_scaling_info$`scaled:scale`
      cur_m            = cur_scaling_info$`scaled:center`
      cur_sd           = cur_sd[names(cur_sd) %in% pred_to_control]
      cur_m            = cur_m[names(cur_m) %in% pred_to_control]
      cur_test_data_coefs[names(cur_sd)] = scale(cur_test_data_coefs[names(cur_sd)],center = cur_m,scale=cur_sd)
    }
  }
  
  stopifnot(is.null(row.names(cur_test_data_coefs)) == F)
  test_data_coefs[[kk]] = cur_test_data_coefs
}

# MAKING PREDICTIONS ==========================================================
predictions        = c()
responses          = c()
predictions_list   = list() # to tally after each CV-fold (no pred pooling)
responses_list     = list()
predictions_orig   = c()
predictions_orig_l = list()
truth              = c()
truth_list         = list() # to tally after each CV-fold (no pred pooling)
for (kk in 1:length(flds)) {
  cur_test_params      = test_data_coefs[[kk]]
  cur_test_params$HCPG = NULL
  disp('Trying prediction now!')
  if (which_ML == 'ML') {
    # linear glmnet case
    if (length(full_mod_cv[[kk]]) == 1 & is.character(full_mod_cv[[kk]][1]) & full_mod_cv[[kk]][1] == 'no_model') {
      cur_sub      = row.names(cur_test_params)
      cur_response = randn(1,length(cur_sub))
      cur_predic   = cur_response
    } else if (any(names(full_mod_cv[[kk]]) %in% 'winning_model')) {
      # glmnet used to predict from exp params
      stopifnot('Not implemented correctly: glmnet used to predict from exp params')
      cur_full_mod    = full_mod_cv[[kk]]$winning_model
      tmp_x           = as.matrix(cur_test_params[grep(names(cur_test_params),pattern="pred_")])
      cur_predic      = predict.cv.glmnet(cur_full_mod, s="lambda.min",newx=tmp_x)
      cur_response    = as.numeric(cur_predic)
    } else {
      cur_full_mod    = full_mod_cv[[kk]]
      if (any(class(cur_full_mod) == 'glmRob')) {
        # default glm (no ML)
        # glmnet used to predict from exp params
        cur_predic      = predict.lmRob(cur_full_mod,newdata=cur_test_params) 
      } else if (is.numeric(cur_full_mod)) {
        test_data     = as.matrix(cur_test_params)
        if (length(cur_full_mod) == 0) {
          cur_full_mod = t(zeros(1,length(test_data[1,])+1))
        }
        cur_predic    = t(as.matrix(cur_full_mod)) %*% t(cbind(ones(length(test_data[,1]),1),test_data))
      } else {
        # default glm (no ML)
        stop('Predict from glm does not work right now.')
        cur_predic      = predict.glm(cur_full_mod,newdata=cur_test_params)
      }
      cur_response    = as.numeric(cur_predic)
    }
    
    # predicting the test data's group
    cur_predic_orig = cur_predic
    cur_predic      = ifelse(cur_predic > 0,"PG","HC")
    predictions = c(predictions,cur_predic)
  } else if (which_ML == 'SVM') {
    cur_predic   = predict(full_mod_cv[[kk]],newdata = cur_test_params)
    predictions  = c(predictions,as.character(cur_predic))
    cur_response = predict(full_mod_cv[[kk]],newdata = cur_test_params, decision.values = T)
    cur_response = attributes(cur_response)
    cur_response = as.numeric(cur_response$decision.values)
  }
  
  # collect the responses (for ROC)
  responses = c(responses,cur_response)
  
  # the truth
  cur_sub   = row.names(cur_test_params)
  criterion = "HCPG"
  cur_truth = as.character(dat_match[dat_match$VPPG %in% cur_sub,criterion])
  truth     = c(truth,cur_truth)
  
  # collect for by-fold calc of confusion matrix
  responses_list[[kk]]   = cur_response
  predictions_list[[kk]] = cur_predic
  truth_list[[kk]]       = cur_truth
}

# EVALUATING: GETTING CV SCORES OUTER  CV =====================================
disp('Using pooled predictions and responses. Better for AUC.')
res_list   = list()
res_list_m = list()

# prep the pred, truth list
pd_list = list()
for (ff in 1:length(truth_list)) {
  cur_entry = list()
  cur_entry$predictions = as.character(predictions_list[[ff]])
  cur_entry$truth       = as.character(truth_list[[ff]])
  pd_list[[ff]]         = cur_entry
}

# unit test pd_list's truth should always be 50/50
for (ff in 1:length(pd_list)) {
  check_half = table(pd_list[[ff]]$truth)
  if (check_half[1] != check_half[2]) {
    stop('Truth in current fold is not 50/50!!')
  }
}

# function to get confusion matrix
f = function(x) {
  # factor
  x[[1]] = factor(x[[1]],levels=c('HC','PG'),labels=c('HC','PG'))
  x[[2]] = factor(x[[2]],levels=c('HC','PG'),labels=c('HC','PG'))
  
  cm = confusionMatrix(x[[1]],x[[2]])
  return(c(cm$overall[1],cm$byClass[1],cm$byClass[2]))
}
cm_collect    = lapply(pd_list,FUN = f)
cm_collect    = unlist(cm_collect)
cm_collect    = matrix(cm_collect,ncol = 3,byrow = T)
accura        = round(mean(cm_collect[,1]), digits = 4)
sensitiv      = round(mean(cm_collect[,2]), digits = 4)
specif        = round(mean(cm_collect[,3]), digits = 4)

criterion     = "HCPG"
accuracy      = round(mean(truth==predictions), digits = 4)
# unit test that pooled and not pooled accura is the same
if (accura != accuracy) {
  warning('Pooled accuracy is not the same as aggregated accuracy!')
  pool_ok = 0
}
cur_xtab      = data.frame(truth,predictions)
cur_xtab      = xtabs(~predictions+truth,data=cur_xtab)
cur_xtab_m    = as.matrix(cur_xtab)
if (mean(dim(cur_xtab_m) == c(2,2)) != 1) {
  acc = data.frame(accuracy)
} else {
  cur_sens    = round(sensitivity(cur_xtab), digits = 4)
  if (cur_sens != sensitiv) {
    warning('Pooled sensitivity is not the same as aggregated sensitivity!')
  }
  cur_spec    = round(specificity(cur_xtab), digits = 4)
  if (specif != cur_spec) {
    warning('Pooled specificity is not the same as aggregated specificity!')
  }
  disp('returning unpooled acc, spec, sens, auc')
  acc         = data.frame(accuracy,cur_sens,cur_spec)
}

# AUC
tp_list  = list()
cur_roc  = roc(truth,responses)

tp              = data.frame(truth, predictions)
res_list$xtab   = cur_xtab
res_list$acc    = acc
res_list$tp     = tp
res_list$xtab_m = cur_xtab_m
res_list$roc    = cur_roc
#res_list$roc    = rocs
res_list$auc    = cur_roc$auc

CV_res          = res_list

# REPORTING ===================================================================
# reporting
disp("#############")
disp("# CV REPORT #")
disp("#############")
disp("")
disp("CV of sev_compl and pred_grp yielded...")
disp("For full model...")
print(CV_res$acc)
disp("AUC:")
print(round(CV_res$auc,digits = 4))
