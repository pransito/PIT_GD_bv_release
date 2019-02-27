## PREAMBLE ==================================================================
# script to do the reporting for estimated classification
# run before: (import_data.R); select_study.R
# in classifier_setting.R: set the basic parameters for reporting

## GET REPORTING SETTINGS =====================================================
setwd(root_wd)
setwd('01_classification/')
source('classifier_settings.R')
stopifnot(do_report_no_added_feat | do_report_with_added_feat | do_report_feat_only)

## REPORTING: PREPARATION =====================================================
# loads the results files
cur_home = getwd()
setwd(path_res_classif)
setwd('results')
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

# REPORTING: p-values =========================================================
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
  Ha_auc                 = auc
  Ha_auc                 = rep_len(Ha_auc,length.out = length(auc_p))
  cur_dat_gl             = data.frame(null_classifier = auc_p,full_classifier = Ha_auc,classifier = 'elastic net')
  cur_dat                = rbind(cur_dat_gl) #rbind(cur_dat_be,cur_dat_gl,cur_dat_sv)
  cur_dat                = melt(cur_dat,id.vars = c('classifier'))
  cur_dat$AUC_ROC        = cur_dat$value
  cur_dat$value          = NULL
  cur_dat$algorithm      = cur_dat$classifier
  cur_dat$classifier     = cur_dat$variable  
  cur_dat$classifier     = agk.recode.c(cur_dat$classifier,c('full_classifier','null_classifier'),c('full','null'))
  
  # plot
  p = ggplot(cur_dat,aes(x=AUC_ROC, fill=classifier)) + geom_density(alpha=0.25)
  p = p + ggtitle('AUC densities for elastic net classifier compared to null-classifier')
  #p = p + facet_grid(algorithm ~ .)
  p = p + geom_vline(aes(xintercept = mean(auc)),colour = 'green',size= 1.5)
  p = p + theme_bw()
  p = p + theme(axis.text=element_text(size=14, face = "bold"),
                axis.title=element_text(size=20,face="bold"))
  p = p + theme(plot.title = element_text(size=22))
  p = p + theme(legend.text = element_text(size=18))
  p = p + theme(legend.title= element_text(size=18))
  print(p)
}

## REPORTING: INSPECTION OF CLASSIFIER ========================================
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

labels_sources = c("catgambling","catnegative","catpositive","Int_behav_model","Int_classifier","edu_years") 
ci_res$coef    = factor(ci_res$coef)
labels_betas   = agk.recode(levels(ci_res$coef),labels_sources,
                            c("gambling cues","negative cues","positive cues","intercept behavioral model","intercept of classifier","years of education"))
ci_res$coef    = factor(ci_res$coef, levels = levels(ci_res$coef), labels = labels_betas)

p = ggplot(data = ci_res, aes(coef,mean))
p = p+geom_bar(stat="identity")
p = p + geom_errorbar(aes(ymin=lower,ymax=upper), size=1.3, color=cbbPalette[4],
                      width = 0) + ylab("mean (95% CI)\n")

p <- p + ggtitle("Estimated regression weights of winning model")
p = p + theme(text = element_text(size=25),
              axis.text.x = element_text(angle=45, hjust=1)) 
p = p + xlab("regression weights")
print(p)

## HAUFE TRANFORMATION OF CLASSIFIER ==========================================
# this is an experimental part

agk.haufe.transformation = function(X,win_mod_coefs) {
  # function to perform the Haus tranformation from linear classifier to intepretable univariate group differences
  cov_X = cov(X)
  
  # using all weight vectors
  all_ws              = win_mod_coefs
  all_ws$X.Intercept. = NULL
  all_ws              = as.data.frame(lapply(all_ws,as.numeric))
  stopifnot(colnames(cov_X) == names(all_ws))
  
  # get the importance at every channel (predictor)
  all_As = cov_X %*% as.matrix(as.numeric(all_ws[1,]))
  for (aa in 2:length(all_ws[,1])) {
    all_As = cbind(all_As,cov_X %*% as.matrix(as.numeric(all_ws[aa,])))
  }
  
  # get the importance mean, and ci bootstrapped (quantile!)
  all_As      = data.frame(t(all_As))
  all_As      = lapply(all_As,agk.mean.quantile.c,lower = 0.025,upper=0.975)
  all_As      = data.frame(all_As)
  all_As      = t(all_As)
  all_As      = as.data.frame(all_As)
  all_As$coef = row.names(all_As)
  
  return(all_As)
}

if (which_study == 'MRI') {
  
  # short names and grouping
  ci_res = agk.mri.shorter.and.grouped.names(ci_res)
  
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
  
  # the top betas
  ci_res_ordered = ci_res[order(abs(ci_res$mean),decreasing = T),]
  message('The strongest betas are:')
  print(ci_res_ordered[1:4,])
  
  # get the mixing matrix A to interpret not W (the unmixing vector), but the actual channel behavior
  # Haufe et al. 2014, Eq. 7
  X           = feature_clusters[[2]]
  X$edu_years = as.numeric(agk.recode.c(row.names(X),dat_match$VPPG,dat_match$edu_years))
  X$HCPG      = NULL
  X           = agk.scale.ifpossible(X)
  
  # core Haufe transformation
  all_As = agk.haufe.transformation(X,win_mod_coefs)
  
  # shorter names and grouping
  all_As = agk.mri.shorter.and.grouped.names(all_As)
  
  # plot
  p = ggplot(data = all_As, aes(coef,mean))
  p = p+geom_bar(stat="identity")
  p = p + geom_errorbar(aes(ymin=lower,ymax=upper), size=1.3, color=cbbPalette[4],
                        width = 0) + ylab("mean (95% CI over CV rounds)\n\n\n")
  
  p <- p + ggtitle("Estimated predictor importance with 95% quantiles") + theme_bw()
  p = p + theme(text = element_text(size=15),
                axis.text.x = element_text(angle=45, hjust=1,face='bold')) 
  p = p + facet_wrap(~ grouping, scales = "free_x",nrow=3,ncol=2)
  print(p)
  
  # the top betas
  all_As_ordered = all_As[order(abs(all_As$mean),decreasing = T),]
  message('The strongest betas are:')
  print(all_As_ordered[1:4,])
} else if (which_study == 'POSTPILOT_HCPG') {
  # display the Haufe transformed A
  X              = featmod_coefs[[winning_mod]]
  X$smoking_ftdt = as.numeric(agk.recode.c(row.names(X),dat_match$VPPG,dat_match$smoking_ftdt))
  X$HCPG         = NULL
  X              = agk.scale.ifpossible(X)
  
  # core Haufe transformation
  all_As = agk.haufe.transformation(X,win_mod_coefs)
  
  # plot
  p = ggplot(data = all_As, aes(coef,mean))
  p = p+geom_bar(stat="identity")
  p = p + geom_errorbar(aes(ymin=lower,ymax=upper), size=1.3, color=cbbPalette[4],
                        width = 0) + ylab("mean (95% CI over CV rounds)\n\n\n")
  
  p <- p + ggtitle("Estimated predictor importance with 95% quantiles") + theme_bw()
  p = p + theme(text = element_text(size=15),
                axis.text.x = element_text(angle=45, hjust=1,face='bold')) 
  print(p)
  
  # the top betas
  all_As_ordered = all_As[order(abs(all_As$mean),decreasing = T),]
  message('The strongest betas are:')
  print(all_As_ordered[1:4,])
}
