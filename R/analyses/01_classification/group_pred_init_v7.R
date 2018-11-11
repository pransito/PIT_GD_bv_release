agk.group.pred.init = function() {
  ## PREAMBLE ===================================================================
  # Version 6.0
  # predict group (HC vs. PG)
  # using all LA models with or without cat modulation
  # selection of which model to extract exp params is
  # part of training (inner CV will be used)
  # this script initializes the data; gets the behav model parameters
  
  # initialization script to set parameters and get data
  # here using on lmlist to get the single-subject exp paramters using different
  # single-subject modeling strategies
  
  # author: Alexander Genauck
  # email:  alexander.genauck@charite.de
  # date:   09.11.2018
  
  ## PROCESS PREPS ==============================================================
  # estimate models backup
  est_models_bcp = est_models
  
  # getting data and subsetting if desired
  data_pdt = data_pdt_bcp_study_selected
  dat_match = dat_match_bcp_study_selected
  
  # plot ratings
  setwd(root_wd)
  source("01_classification/plot_ratings_betting_behav.R")
  
  # prep if peripheral physiology or ratings should be added
  if ((add_cr_pp  == 1 || add_cr_ra  == 1) & which_study != 'MRI') {
    # excluding subjects because of missing in pp or ra
    # prepping data frames for ra and pp
    cur_path = getwd()
    setwd(root_wd)
    source("01_classification/get_phys_and_rating_params_and_plots_v2.R")
  } else if (add_cr_pp  == 1 & which_study == 'MRI') {
    # data is already gathered and in workspace
  } else {
    # do nothing
  }
  
  if (which_study == 'MRI' & add_cr_pp_ma) {
    #cr_agg_pp = cr_agg_pp_readin
    cr_agg_pp = cr_agg_pp_r_MRI
    # clean out subject variable
    row.names(cr_agg_pp) = cr_agg_pp$subject
    cr_agg_pp$subject    = NULL
    
    # standardize within subject
    cr_agg_pp = as.data.frame(t(scale(t(cr_agg_pp))))
    
    # reduce set further by computing means
    cr_agg_pp_m = cr_agg_pp
    # SS cue reactivity variables; mean of L and R (so mean of two variables)
    all_cr_names   = names(cr_agg_pp_m)[grep('SS__grp01_',names(cr_agg_pp_m))]
    all_cr_names_L = all_cr_names[grep('_Left_',all_cr_names)]
    for (nn in 1:length(all_cr_names_L)) {
      cur_name                     = all_cr_names_L[nn]
      cur_R                        = gsub('_Left_','_Right_',cur_name)
      cr_agg_pp_m[[cur_name]]      = (cr_agg_pp_m[[cur_name]] + cr_agg_pp_m[[cur_R]]) /2
      cr_agg_pp_m[cur_R]           = NULL
      # get the name we need to change
      cur_ind                     = which(names(cr_agg_pp_m) == cur_name)
      names(cr_agg_pp_m)[cur_ind] = gsub('_Left_','_LR_',names(cr_agg_pp_m[cur_name]))
    }
    # SS gPPI PIT variables; mean of L and R (so mean of four variables)
    all_gppi_names         = names(cr_agg_pp_m)[grep('SS__PPI_',names(cr_agg_pp_m))]
    all_gppi_names_L       = all_gppi_names[grep('_PPI_L',all_gppi_names)]
    all_gppi_names_L_L_tar = all_gppi_names[grep('_ROI_L',all_gppi_names_L)]
    for (nn in 1:length(all_gppi_names_L_L_tar)) {
      cur_name                     = all_gppi_names_L_L_tar[nn]             # left source left target
      cur_R_source                 = gsub('_PPI_L_','_PPI_R_',cur_name)     # right source left target
      cur_R_source_L_tar           = gsub('_ROI_L_','_ROI_R_',cur_R_source) # right source right target
      cur_L_source_L_tar           = gsub('_ROI_L_','_ROI_R_',cur_name)     # left source right target
      
      # new name and calculation of mean
      new_name                     = gsub('_PPI_L_','_PPI_LR_',cur_name)
      new_name                     = gsub('_ROI_L_','_ROI_LR_',new_name)
      cr_agg_pp_m[[new_name]]      = (cr_agg_pp_m[[cur_name]] + cr_agg_pp_m[[cur_R_source]] + cr_agg_pp_m[[cur_R_source_L_tar]] + cr_agg_pp_m[[cur_L_source_L_tar]])/4
      
      # delete unneeded variables
      cr_agg_pp_m[cur_name]            = NULL
      cr_agg_pp_m[cur_R_source]        = NULL
      cr_agg_pp_m[cur_R_source_L_tar]  = NULL
      cr_agg_pp_m[cur_L_source_L_tar]  = NULL
    }
    # some by hand
    if (fmri_extr == 'ngm' | fmri_extr == 'glc') {
      cr_agg_pp_m$SS__PPI_LR_Amy_noCov_PPI_PicGamOnxaccXgam_ROI_R_MOrG = (cr_agg_pp_m$SS__PPI_L_Amy_noCov_PPI_PicGamOnxaccXgam_ROI_R_MOrG + cr_agg_pp_m$SS__PPI_R_Amy_noCov_PPI_PicGamOnxaccXgam_ROI_R_MOrG)/2
      cr_agg_pp_m$SS__PPI_LR_Amy_noCov_PPI_PicGamOnxaccXneg_ROI_R_MOrG = (cr_agg_pp_m$SS__PPI_L_Amy_noCov_PPI_PicGamOnxaccXneg_ROI_R_MOrG + cr_agg_pp_m$SS__PPI_R_Amy_noCov_PPI_PicGamOnxaccXneg_ROI_R_MOrG)/2
      cr_agg_pp_m$SS__PPI_LR_Amy_noCov_PPI_PicGamOnxaccXpos_ROI_R_MOrG  = (cr_agg_pp_m$SS__PPI_L_Amy_noCov_PPI_PicGamOnxaccXpos_ROI_R_MOrG + cr_agg_pp_m$SS__PPI_R_Amy_noCov_PPI_PicGamOnxaccXpos_ROI_R_MOrG)/2
    } else if (fmri_extr == 'val') {
      cr_agg_pp_m$SS__PPI_LR_Amy_noCov_PPI_PicGamOnxvalXgam_ROI_R_MOrG = (cr_agg_pp_m$SS__PPI_L_Amy_noCov_PPI_PicGamOnxvalXgam_ROI_R_MOrG + cr_agg_pp_m$SS__PPI_R_Amy_noCov_PPI_PicGamOnxvalXgam_ROI_R_MOrG)/2
      cr_agg_pp_m$SS__PPI_LR_Amy_noCov_PPI_PicGamOnxvalXneg_ROI_R_MOrG = (cr_agg_pp_m$SS__PPI_L_Amy_noCov_PPI_PicGamOnxvalXneg_ROI_R_MOrG + cr_agg_pp_m$SS__PPI_R_Amy_noCov_PPI_PicGamOnxvalXneg_ROI_R_MOrG)/2
      cr_agg_pp_m$SS__PPI_LR_Amy_noCov_PPI_PicGamOnxvalXpos_ROI_R_MOrG  = (cr_agg_pp_m$SS__PPI_L_Amy_noCov_PPI_PicGamOnxvalXpos_ROI_R_MOrG + cr_agg_pp_m$SS__PPI_R_Amy_noCov_PPI_PicGamOnxvalXpos_ROI_R_MOrG)/2
    }
    
    # delete unneeded variables
    cr_agg_pp_m[names(cr_agg_pp_m)[grep('PPI_._',names(cr_agg_pp_m))]] = NULL
    cr_agg_pp = cr_agg_pp_m
  }
  
  # cleaning cr_agg_pp
  if (add_cr_pp_ma & which_study == 'MRI') {
    if (regress_out_covs) {
      if (!exists('cr_agg_pp_cleaned')) {
        vars_to_cov       = pred_to_clean
        res               = agk.clean.vars(cr_agg_pp,dat_match,vars_to_cov,clean_inf_crit)
        cr_agg_pp_cleaned = res$cr_agg_pp_cleaned
        cr_agg_pp         = res$cr_agg_pp
      } else {
        cr_agg_pp = cr_agg_pp_cleaned
      }
    }
  }
  
  # set deviance measure
  cur_family = "binomial"
  type_measure = type_measure_binomial
  
  # cross-valid results
  CV_res = list()
  
  # cons
  if (do_data_inv == 0) {
    contrasts(data_pdt$cat) = contrasts(data_pdt$cat)
  }
  
  # control variables
  if (length(pred_to_control) == 1) {
    # exactly one thing to control for
    # then in ridging need to add a random regressor
    do_con = 1
  } else if (length(pred_to_control) > 1) {
    # more than one to control for
    do_con = 2
  } else {
    # will no con for anything
    do_con = 0
  }
  
  ## CHECK AVAILABLE DATA =======================================================
  # checking data which is there
  # data pdt
  disp("In current data_pdt I have")
  all_subs = unique(data_pdt$subject)
  grp_var  = agk.recode.c(all_subs,data_pdt$subject,data_pdt$HCPG)
  subgrp   = data.frame(subject = all_subs,group = grp_var)
  print(table(subgrp$group))
  
  # dat_match
  disp(paste("Of the", length(all_subs), "subjects, I have dat_match info on",
             sum(dat_match$VPPG %in% all_subs), "subjects"))
  
  # dat_match KFG, SOGS, BIG check
  cur_dat_match = dat_match[dat_match$VPPG %in% all_subs,]
  cur_dat_match = cur_dat_match[c("KFG","SOGS",'BIG')]
  disp("Are there any NA's in KFG or SOGS or BIG?")
  print(table(is.na(cur_dat_match)))
  
  ## MODEL PARAMS PER SUB =======================================================
  # report to user
  disp('')
  disp('Getting the behavioral model parameters per subject.')
  
  # function for charpentier model
  # maybe not use anymore cause can be done in the general get lambda function
  transform_coefs = function(lml) {
    if (!is.data.frame(lml)) {
      lml      = coef(lml)
    }
    lml$mu   = lml$gain
    lml$LA   = -(lml$loss/lml$gain)
    lml$gain = NULL
    lml$loss = NULL
    return(lml)
  }
  
  
  # for charpentier per category
  data_pdt_neu       = data_pdt[data_pdt$cat == "neutral",]
  data_pdt_gam       = data_pdt[data_pdt$cat == "gambling",]
  data_pdt_neg       = data_pdt[data_pdt$cat == "negative",]
  data_pdt_pos       = data_pdt[data_pdt$cat == "positive",]
  
  # estimate models for params extraction
  if (acc_num == 1) {
    needed_family = 'gaussian'
  } else {
    needed_family = 'binomial'
  }
  
  if (est_models == 1) {
    
    if (!mod_physio_val) {
      # legacy; but needs to stay for now
      clmList = lmList
      
      # lmlist
      a         = clmList(accept_reject ~ 1 | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      ac        = clmList(accept_reject ~ cat | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      laec      = clmList(accept_reject ~ gain+loss+ed_abs+cat | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      lac       = clmList(accept_reject ~ gain+loss+cat | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      laeci     = clmList(accept_reject ~ (gain+loss+ed_abs)*cat | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      laci      = clmList(accept_reject ~ (gain+loss)*cat | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      lae       = clmList(accept_reject ~ gain+loss+ed_abs | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      la        = clmList(accept_reject ~ gain+loss | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      # ratio lmlist (no mu, please! drop it! it's enough!)
      lar       = clmList(accept_reject ~ ratio | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      larc      = clmList(accept_reject ~ ratio + cat | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      larci     = clmList(accept_reject ~ ratio*cat | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      
      # charpentier lmlist (same as la but without intercept, then mu and lambda can be analytically computed)
      laCh      = clmList(accept_reject ~ 0 + gain + loss | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      # Charpentier per category
      laChci    = clmList(accept_reject ~ 0 + gain + loss + gain.catgambling + gain.catnegative + gain.catpositive +       
                            loss.catgambling + loss.catnegative + loss.catpositive | subject,
                          data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      
      # packing (in order of complexity)
      fm = list()
      fm$a       = a
      fm$ac      = ac
      fm$lar     = lar
      fm$laCh    = laCh
      fm$la      = la
      fm$lae     = lae
      fm$larc    = larc
      fm$lac     = lac
      fm$laec    = laec
      fm$larci   = larci
      fm$laci    = laci
      fm$laeci   = laeci
      fm$laChci  = laChci
      
    } else {
      # with mod_physio_val
      # so valence and arousal can modulate choice behavior
      
      # legacy; but needs to stay for now
      clmList = lmList
      
      # lmlist
      a         = clmList(accept_reject ~ 1 | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      ap        = clmList(accept_reject ~ scr_arousal + cozy_valence | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      laep      = clmList(accept_reject ~ gain+loss+ed_abs+scr_arousal + cozy_valence | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      lap       = clmList(accept_reject ~ gain+loss+scr_arousal + cozy_valence | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      laepi     = clmList(accept_reject ~ (gain+loss+ed_abs)*(scr_arousal+cozy_valence) | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      lapi      = clmList(accept_reject ~ (gain+loss)*(scr_arousal+cozy_valence) | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      lae       = clmList(accept_reject ~ gain+loss+ed_abs | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      la        = clmList(accept_reject ~ gain+loss | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      # ratio lmlist (no mu, please! drop it! it's enough!)
      lar       = clmList(accept_reject ~ ratio | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      larp      = clmList(accept_reject ~ ratio + (scr_arousal+cozy_valence) | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      larpi     = clmList(accept_reject ~ ratio*(scr_arousal+cozy_valence) | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      
      # charpentier lmlist (same as la but without intercept, then mu and lambda can be analytically computed)
      laCh      = clmList(accept_reject ~ 0 + gain + loss | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      # Charpentier per category [TODO!]
      #laChci    = clmList(accept_reject ~ 0 + gain + loss + gain.catgambling + gain.catnegative + gain.catpositive +       
      #                      loss.catgambling + loss.catnegative + loss.catpositive | subject,
      #                    data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
      
      # packing (in order of complexity)
      fm = list()
      fm$a       = a
      fm$ap      = ap
      fm$lar     = lar
      fm$laCh    = laCh
      fm$la      = la
      fm$lae     = lae
      fm$larp    = larp
      fm$lap     = lap
      fm$laep    = laep
      fm$larpi   = larpi
      fm$lapi    = lapi
      fm$laepi   = laepi
      #fm$laChci  = laChpi
      
    }
  }
  
  # backup the fm
  fm_bcp = fm
  
  ## RIDGE THE EXP BEHAV MODELS =================================================
  if (ridge_behav_models == T) {
    if (ridge_behav_models_anew == T) {
      # warning
      message('update glmnet to perhaps avoid having to provide fixed lambdas directly')
      # what alphas to use?
      cur_alpha = ridge_bm_alphas
      
      disp('Ridging the extracted exp. parameters for all models within every subject')
      # body of function to parallelize
      cur.mod.ridging = function(kk) {
        # decide whether this model is ridgeable
        cur_coefs = coef(fm[[kk]])
        if (!isempty(grep('Intercept',names(cur_coefs)))) {min_length = 3} else {min_length = 2} 
        if (length(names(cur_coefs)) < min_length) {
          fm[[kk]] = coef(fm[[kk]])
          next
        }
        
        disp(paste('Ridging model...',names(fm)[kk]))
        cur_mod    = lapply(fm[[kk]],FUN=agk.glmnet.lm.cvalpha,fam_arg='binomial',
                            type_measure='auc',lambda=des_lambdas,nfolds=5,
                            alphas = cur_alpha,verbosity = 1)
        coef_listr = cur_mod[[1]]$coef
        
        for (ii in 2:length(cur_mod)) {
          coef_listr = rbind(coef_listr,cur_mod[[ii]]$coef)
        }
        row.names(coef_listr) = NULL
        coef_listr = data.frame(coef_listr)
        colnames(coef_listr) = names(coef(fm[[kk]]))
        
        # packing back
        return(coef_listr)
      }
      
      total = length(fm)
      cl    = parallel::makeCluster((detectCores()-1))
      registerDoSNOW(cl)
      fm_ridge = foreach(kk=1:total, .packages=c('glmnetUtils','cvTools'),.verbose=T,.export = c()) %dopar% {
        cur.mod.ridging(kk)
      } 
      stopCluster(cl)
      
      # packing back, saving
      fm = fm_ridge
      save(file = 'ridged_behav_models.RData', list = c('fm'))
    } else {
      disp('loading the ridged behav params')
      load('ridged_behav_models.RData')
    }
  }
  
  # getting all the subs
  all_subs = names(fm[[1]])
  
  # Compute and ammend lambdas (as extra models)
  # state which models should get LA appended
  LA_append = c("la","lae","lac","laec","laci","laeci","laCh","laChci") 
  fm_LA = list()
  for (ii in which(names(fm) %in% LA_append)) {
    fm_LA[[ii]]      = agk.get.lambda(fm[[ii]],'gain','loss','(:|\\.)')
    names(fm_LA)[ii] = names(fm)[ii]
  }
  
  # drop models which haven't got LA appended
  fm_LA = fm_LA[!unlist(lapply(fm_LA,is.null))]
  
  # bind to none_LA models
  names(fm_LA) = paste0(names(fm_LA),'_LA')
  fm = c(fm,fm_LA)
  cur_is_null = as.logical(unlist(lapply(fm,FUN=is.null))) == FALSE
  fm          = fm[cur_is_null]
  
  # coerce all fm to data frame (extract the coefs)
  for (ii in 1:length(fm)) {
    df = fm[[ii]]
    if (is.data.frame(df)) {
      df = df
    } else {
      df = coef(df)
    }
    if (!is.data.frame(df)) {
      stop('Cur df is not a data frame nor can I extract coefficients data frame from df object.')
    }
    fm[[ii]] = df
  }
  
  # put back the row.names
  for (ff in 1:length(fm)) {
    row.names(fm[[ff]]) = row.names(fm[[1]])
  }
  
  # Count and order by complexity
  disp('Ordering models by complexity.')
  complexity = unlist(lapply(fm,length))
  fm         = fm[order(complexity)]
  
  # names of models
  names_models = names(fm)
  
  # Scaling
  if (initial_scale) {
    disp('Scaling the model parameters in fm.')
    fm = lapply(fm,FUN=agk.scale.ifpossible)
  }
  
  # Within subject z-scoring of feature vectors
  if (behav_within_sub_z) {
    for (ff in 1:length(fm)) {
      if (length(fm[[ff]]) >= 2) {
        fm[[ff]] = as.data.frame(t(scale(t(fm[[ff]]))))
      }
    }
  }
  
  # Regress out covariates
  if (regress_out_covs) {
    for (ff in 1:length(fm)) {
      vars_to_cov       = pred_to_clean
      res               = agk.clean.vars(fm[[ff]],dat_match,vars_to_cov,clean_inf_crit)
      fm[[ff]]          = res$cr_agg_pp
    }
  }
  
  # getting the coefficient data frames per model: packing into featmod_coefs
  disp('Packing fm into featmod_coefs')
  featmod_coefs = list()
  
  for (ii in 1:length(fm)) {
    # lmlist coef extraction of complete model
    if (is.data.frame(fm[[1]])) {
      tmp = fm[[ii]]
    } else {
      tmp = coef(fm[[ii]])
    }
    
    # feature: exp params
    if (names(tmp)[1] == '(Intercept)') {
      names(tmp)[1] = 'Intercept'
    }
    tmp_exp       = tmp
    
    # packing
    featmod_coefs[[ii]] = tmp_exp
  }
  
  # give a name
  names(featmod_coefs) = names(fm)
  
  # prep the coefs tables
  for (ii in 1:length(featmod_coefs)) {
    tmp                 = featmod_coefs[[ii]]
    names(tmp)          = gsub(pattern = ":",replacement = "_of_",names(tmp))
    #names(tmp)          = paste0("pred_",names(tmp))
    #tmp$subject         = row.names(tmp)
    featmod_coefs[[ii]] = tmp
  }
  
  # getting the dependent variables ((PCA on) KFG and SOGS, BIG, HCPG)
  for (ii in 1:length(featmod_coefs)) {
    tmp      = featmod_coefs[[ii]]
    #tmp$KFG  = as.numeric(agk.recode.c(tmp$subject,dat_match$VPPG,dat_match$KFG))
    tmp$HCPG = as.factor(agk.recode.c(row.names(tmp),dat_match$VPPG,dat_match$HCPG))
    featmod_coefs[[ii]] = tmp
  }
  
  # make a back up of the featmod_coefs
  featmod_coefs_bcp = featmod_coefs
  
  # # make the formulas for model selection
  # featmod_sel_forms = list()
  # for (ii in 1:length(featmod_coefs)) {
  #   tmp        = featmod_coefs[[ii]]
  #   pred_vars  = names(tmp)[grep("pred_",names(tmp))]
  #   form_cmpl  = as.formula(paste('HCPG','~', paste(pred_vars, collapse=" + ")))
  #   featmod_sel_forms[[ii]] = form_cmpl
  # }
  
  # getting other features ready for analysis
  # scaling (not done in 'with outer CV')
  if (initial_scale) {
    disp('Scaling additional features overall (NOT RECOMMENDED!).')
  } else {
    disp('Not scaling additional features overall.')
  }
  
  if (add_cr_pp  == 1 || add_cr_ra  == 1) {
    # ADDING OTHER FEATURES (E.G. PHYSIO)
    # getting the rat and phys feature cluster
    tmp         = tmp_exp
    exp_params  = names(tmp)
    tmp_rat     = agk.merge.df.by.row.names(tmp,cr_agg_ra)
    tmp_rat     = tmp_rat[!names(tmp_rat) %in% exp_params]
    rat_params  = names(tmp_rat)
    rat_params  = rat_params[!rat_params %in% c('subject')]
    tmp_rat$subject = NULL
    tmp_phys    = agk.merge.df.by.row.names(tmp_rat,cr_agg_pp)
    tmp_phys    = tmp_phys[!names(tmp_phys) %in% rat_params]
    tmp_phys$subject    = NULL
    
    # scaling (not done in with outer CV)
    if (initial_scale) {
      tmp_rat             = as.data.frame(scale(tmp_rat))
      tmp_phys            = as.data.frame(scale(tmp_phys))
    }
    
    # make sure that only used subjects are in additional pp
    tmp_phys = subset(tmp_phys, row.names(tmp_phys) %in% unique(data_pdt$subject))
    tmp_rat  = subset(tmp_rat, row.names(tmp_rat) %in% unique(data_pdt$subject))
  }
  
  # packing selected model coefs plus other feature clusters
  feature_clusters      = list()
  feature_clusters[[1]] = NULL
  ct                    = 2
  if (add_cr_ra) {
    feature_clusters[[ct]] = tmp_rat
    ct                     = ct + 1
  }
  if (add_cr_pp) {
    feature_clusters[[ct]] = tmp_phys
    ct                     = ct + 1
  }
  
  # prep the coef table
  for (ii in 1:length(feature_clusters)) {
    if (ii > 1) {
      # ii == 1 is already done
      tmp = feature_clusters[[ii]]
      names(tmp)    = gsub(pattern = ":",replacement = "_of_",names(tmp))
      #names(tmp)    = paste0("pred_",names(tmp))
      #tmp$subject   = row.names(tmp)
      feature_clusters[[ii]] = tmp
    }
  }
  
  # getting the dependent variable (PCA on KFG and SOGS and BIG)
  for (ii in 1:length(feature_clusters)) {
    if (ii > 1) {
      tmp = feature_clusters[[ii]]
      #tmp$KFG  = as.numeric(agk.recode.c(tmp$subject,dat_match$VPPG,dat_match$KFG))
      tmp$HCPG = as.factor(agk.recode.c(row.names(tmp),dat_match$VPPG,as.character(dat_match$HCPG)))
      feature_clusters[[ii]] = tmp
    }
  }
  
  # scale the control variables
  if (initial_scale & !isempty(pred_to_control)) {
    dat_match[[pred_to_control]] = scale(dat_match[[pred_to_control]])
  }
  
  # back up the feature_clusters list
  feature_clusters_bcp = feature_clusters
  
  # return everything
  return(environment())
}