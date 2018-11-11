## PREAMBLE ===================================================================
# script to apply PIT_GD_behav classifier to MRI data
# which_study must be set to MRI in select_study.R first and run
if (init_done == F) {
  stop('You have to first run the select_study.R script with which_study set to "MRI" ')
}

stopifnot(which_study == 'MRI') 

# functions
get.truth = function() {
  sample(c(rep('HC',3),rep('PG',3)))
}

get.truth.4 = function() {
  sample(c(rep('HC',4),rep('PG',4)))
}

get.truth.2 = function() {
  sample(c(rep('HC',2),rep('PG',2)))
}

# set runs of random classification
runs0 = 10000

# under 0
# pooled
all_aucs  = c()
all_aucsl = list()
all_accs  = c()
all_accsl = list()
all_sens  = c()
all_spec  = c()
for (ii in 1:runs0) {
  print(ii)
  inner_truths = c()
  inner_resps  = c()
  # 3
  for (jj in 1:10) {
    # get truth
    inner_truths = c(inner_truths,as.character(get.truth()))
    # get response
    inner_resps  = c(inner_resps,as.numeric(randn(1,6)*10))
  }
  
  # cur_auc
  cur_roc         = roc(inner_truths,inner_resps)
  all_aucs[ii]    = cur_roc$auc
  all_aucsl[[ii]] = cur_roc
  
  # accuracy
  inner_preds     = ifelse(inner_resps<0,'HC','PG')
  all_accsl[[ii]] = inner_truths == inner_preds
  all_accs[ii]    = mean(all_accsl[[ii]])
  
  # sens and spec
  cur_cm = caret::confusionMatrix(table(inner_truths,inner_preds))
  all_sens[ii]    = cur_cm$byClass[1]
  all_spec[ii]    = cur_cm$byClass[2]
}

## use a consensus of ALL models from PDT behav ===============================
setwd(root_wd)
setwd('01_classification/results/1010/')
load('POSTPILOT_HCPG_predGrp1_rounds_noo_noaddfeat.RData')

# get the standardization
# THIS CODE JUST FOR DOCUMENTATION; HAS BEEN DONE BEFORE AND RESULT SAVED
# # first load postpilot data [prep for publication a new workspace]
# setwd('C:/Users/genaucka/Google Drive/Library/01_Projects/PIT_GD/R/analyses/01_classification/results/1010')
# win_mods = agk.recode(win_mods,c('acc'),c('ac'))
# for (mm in 1:length(win_mods)) {
#   pp_b_dat = featmod_coefs[[win_mods[mm]]]
#   pp_b_dat = pp_b_dat[,grep('HCPG',names(pp_b_dat),invert=TRUE)]
#   pp_b_dat = data.frame(pp_b_dat,pred_smoking_ftdt=dat_match$smoking_ftdt)
#   pp_b_dat = scale(pp_b_dat)
#   save(pp_b_dat, file=paste0('POSTPILOT_',win_mods[mm],'_stand.RData'))
# }

# predict with each model
responses = list()
cur_mod_sel_nooCV = agk.recode(cur_mod_sel_nooCV,c('acc'),c('ac'))
for (mm in 1:length(list_winning_model_c_nooCV)) {
  cur_c = list_winning_model_c_nooCV[[mm]]
  cur_l = list_winning_model_l_nooCV[[mm]]
  cur_m = cur_mod_sel_nooCV[mm]
  
  mr_b_dat = featmod_coefs[[cur_m]]
  mr_b_dat = mr_b_dat[,grep('HCPG',names(mr_b_dat),invert=T)]

  # apply the standardization and get decision value
  load(paste0('POSTPILOT_', cur_m,'_stand.RData'))
  pp_scale = attributes(pp_b_dat)

  mr_b_dat = data.frame(mr_b_dat,pred_smoking_ftdt = dat_match$smoking_ftdt)
  mr_b_dat = scale(mr_b_dat,center = pp_scale$`scaled:center`, scale = pp_scale$`scaled:scale`)
  mr_b_dat = data.frame(ones(length(mr_b_dat[,1]),1),mr_b_dat)
  
  # prediction
  responses[[mm]] = t(as.matrix(cur_c)) %*% t(as.matrix(mr_b_dat))
}

# consensus (weighted sum of decision values)
weighted_responses = responses[[1]]
for (mm in 2:length(responses)) {
  weighted_responses = weighted_responses + responses[[mm]]
}
preds     = ifelse(weighted_responses <= 0, 'HC','PG')

acc = mean(preds == dat_match$HCPG)
roc = pROC::roc(as.numeric(dat_match$HCPG),predictor=as.numeric(weighted_responses))
auc = roc$auc
cm  = confusionMatrix(as.factor(preds),dat_match$HCPG)
sen = cm$byClass[1]
spe = cm$byClass[2]

# test
message('p-values for accuracy, AUC, sensitivity, specificity:')
message(' ')
message(1-agk.density_p.c(all_accs,acc))
message(' ')
message(1-agk.density_p.c(all_aucs,auc))
message(' ')
message(1-agk.density_p.c(all_sens,sen))
message(' ')
message(1-agk.density_p.c(all_spec,spe))
message(' ')

# weighted models auc
message('Applying the PIT GD classifiers to the PIT  GD MRI behav data the AUC is:')
message(' ')
all_mod_mean_auc = auc
message(all_mod_mean_auc)

# get auc per classifier
real_aucs = c()
for (aa in 1:length(responses)) {
  real_aucs[aa] = as.numeric(auc(roc(as.numeric(dat_match$HCPG),as.numeric(responses[[aa]]))))
}

## density plots ==============================================================
cur_dat_be = data.frame(random_classifier = all_aucs,mean_auc = rep(all_mod_mean_auc,length(all_aucs)),classifier = 'prev_behav_glmnet')


cur_dat              = rbind(cur_dat_be) #,cur_dat_gl,cur_dat_sv)
cur_dat              = melt(cur_dat,id.vars = c('classifier'))
cur_dat_H_0          = subset(cur_dat,variable == 'random_classifier')
cur_dat_H_0$mean_auc = cur_dat$value[cur_dat$variable == 'mean_auc']
cur_dat              = cur_dat_H_0
cur_dat$AUC_ROC      = cur_dat$value
cur_dat$value        = NULL
p = ggplot(cur_dat,aes(x=AUC_ROC, fill=variable)) + geom_density(alpha=0.25)
#p = p + facet_grid(classifier ~ .) + ggtitle('AUC densities for estimated classifier compared to random classifier')
p = p + ggtitle('AUC densities for estimated classifier compared to random classifier')
p = p + geom_vline(aes(xintercept = mean_auc),colour = 'green',size= 1.5)
p = p + coord_cartesian(xlim = c(0.42, 0.8)) 
p = p + theme_bw()
p = p + theme(axis.text=element_text(size=14, face = "bold"),
              axis.title=element_text(size=20,face="bold"))
print(p)

## two density plots ==============================================================
agk.make.as.long = function(x1,x2) {
  if (length(x1) < length(x2)) {
    x1 = rep_len(x1,length.out = length(x2))
  } else {
    x2 = rep_len(x2,length.out = length(x1))
  }
  return(list(x1 = x1, x2 = x2))
}

Ha_auc               = real_aucs
cur_res              = agk.make.as.long(all_aucs,Ha_auc)
all_aucs             = cur_res$x1
Ha_auc               = cur_res$x2
cur_dat_gl           = data.frame(H0 = all_aucs,Ha_auc = Ha_auc,classifier = 'full_classifier')
cur_dat              = rbind(cur_dat_gl) #rbind(cur_dat_be,cur_dat_gl,cur_dat_sv)
cur_dat              = melt(cur_dat,id.vars = c('classifier'))

# plot
p = ggplot(cur_dat,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
p = p + facet_grid(classifier ~ .) + ggtitle('AUC densities for glmnet classifier on MRI behav data compared to random classifier')
p = p + geom_vline(aes(xintercept = as.numeric(all_mod_mean_auc)),colour = 'green',size= 1.5)
p = p + theme_bw()
p = p + theme(axis.text=element_text(size=14, face = "bold"),
              axis.title=element_text(size=20,face="bold"))
p = p + coord_cartesian(xlim = c(0.4, 0.8)) 
print(p)