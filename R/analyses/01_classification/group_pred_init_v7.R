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
  
  # add novelty factor
  add_novelty_factor = function(data_pdt) {
    
    data_pdt$novel = 1
    
    all_subs = unique(data_pdt$subject)
    
    for (cur_sub in all_subs) {
      novelity_info = as.numeric(base::duplicated(data_pdt$stim[data_pdt$subject == cur_sub]) == TRUE)
      data_pdt$novel[data_pdt$subject == cur_sub] = novelity_info
    }
    
    data_pdt$novel = as.factor(data_pdt$novel)
    return(data_pdt)
  }
  
  data_pdt = add_novelty_factor(data_pdt)
  
  # cut gambling cues
  cut_gam_cues = function(data_pdt) {
    
    all_subs = unique(data_pdt$subject)
    res = list()
    for (ii in 1:length(all_subs)) {
      cur_sub = all_subs[ii]
      cur_dat = subset(data_pdt,subject == cur_sub)
      cur_dat_gam = subset(cur_dat, cat == 'gambling' )
      cur_dat_oth = subset(cur_dat, cat != 'gambling' )
      cur_dat_gam = cur_dat_gam[sample(seq(1,length(cur_dat_gam[,1])),45),]
      cur_dat = rbind(cur_dat_gam,cur_dat_oth)
      res[[ii]] = cur_dat
    }
    
    res_all = res[[1]]
    for (ii in 2:length(res)) {
      res_all = rbind(res_all,res[[ii]])
    }
    
    return(res_all)
  }
  
  
  if (cut_gambling_cues) {
    message('Cutting gambling cues to 45.')
    data_pdt = cut_gam_cues(data_pdt)
  }
  
  # plot ratings
  setwd(root_wd)
  source("01_classification/plot_ratings_betting_behav.R")
  
  # set deviance measure
  cur_family = "binomial"
  type_measure = type_measure_binomial
  
  # cross-valid results
  CV_res = list()
  
  # cons
  contrasts(data_pdt$cat) = contrasts(data_pdt$cat)
  
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
      
      if (!add_novelty) {
      
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
      
      } else {
        message('All single-subject models controlled for novelty')
        # adding the novelty factor
        # lmlist
        a         = clmList(accept_reject ~ novel | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
        ac        = clmList(accept_reject ~ cat + novel | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
        laec      = clmList(accept_reject ~ gain+loss+ed_abs+cat+novel | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
        lac       = clmList(accept_reject ~ gain+loss+cat+novel | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
        laeci     = clmList(accept_reject ~ (gain+loss+ed_abs)*cat+novel | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
        laci      = clmList(accept_reject ~ (gain+loss)*cat+novel | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
        lae       = clmList(accept_reject ~ gain+loss+ed_abs+novel | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
        la        = clmList(accept_reject ~ gain+loss+novel | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
        # ratio lmlist (no mu, please! drop it! it's enough!)
        lar       = clmList(accept_reject ~ ratio+novel | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
        larc      = clmList(accept_reject ~ ratio + cat +novel | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
        larci     = clmList(accept_reject ~ ratio*cat+novel | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
        
        # charpentier lmlist (same as la but without intercept, then mu and lambda can be analytically computed)
        laCh      = clmList(accept_reject ~ 0 + gain + loss+novel | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
        # Charpentier per category
        laChci    = clmList(accept_reject ~ 0 + gain + loss + gain.catgambling + gain.catnegative + gain.catpositive +       
                              loss.catgambling + loss.catnegative + loss.catpositive+novel | subject,
                            data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
        
        
      }
      
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
      warning('NO ADDING OF NOVELTY OF STIMULUS IMPLEMENTED.')
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
  
  # take out the cov of no-interest: novelty of stimulus
  if (add_novelty) {
    cur_fun = function(x) {
      x$novel1 = NULL
      x$novel0 = NULL
      return(x)
    }
    
    fm = lapply(fm,FUN =cur_fun)
  }
  
  # Within subject z-scoring of feature vectors
  if (behav_within_sub_z) {
    for (ff in 1:length(fm)) {
      if (length(fm[[ff]]) >= 2) {
        fm[[ff]] = as.data.frame(t(scale(t(fm[[ff]]))))
      }
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
  
  # packing selected model coefs plus other feature clusters
  feature_clusters      = list()
  feature_clusters[[1]] = NULL
  ct                    = 2

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