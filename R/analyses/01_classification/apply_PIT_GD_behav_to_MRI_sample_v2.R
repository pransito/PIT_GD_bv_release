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

get.truth.1 = function() {
  sample(c('HC','PG'),size = 1)
}

# set runs of random classification
runs0 = 4000

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
  # get 30 in each group
  for (jj in 1:10) {
    # get truth
    inner_truths = c(inner_truths,as.character(get.truth()))
    # get response
    inner_resps  = c(inner_resps,as.numeric(randn(1,6)*10))
  }
  
  # # add 4
  # inner_truths = c(inner_truths,as.character(get.truth.2()))
  # # get response
  # inner_resps  = c(inner_resps,as.numeric(randn(1,4)*10))
  # 
  # # add 1
  # inner_truths = c(inner_truths,as.character(get.truth.1()))
  # # get response
  # inner_resps  = c(inner_resps,as.numeric(randn(1,1)*10))
  
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
setwd(path_res_classif)
setwd('results/1008/')
load('POSTPILOT_HCPG_predGrp1_rounds_noo_noaddfeat.RData')

# #get the standardization
# #THIS CODE JUST FOR DOCUMENTATION; HAS BEEN DONE BEFORE AND RESULT SAVED
# # first load postpilot data [prep for publication a new workspace]
# setwd('C:/Users/genaucka/GitHub/PIT_GD_bv_release/R/analyses/01_classification/results/1003')
# win_mods = cur_mod_sel_nooCV
# win_mods = agk.recode(win_mods,c('acc'),c('ac'))
# for (mm in 1:length(win_mods)) {
#  pp_b_dat = featmod_coefs[[win_mods[mm]]]
#  pp_b_dat = pp_b_dat[,grep('HCPG',names(pp_b_dat),invert=TRUE)]
#  pp_b_dat = data.frame(pp_b_dat,pred_smoking_ftdt=dat_match$smoking_ftdt)
#  pp_b_dat = scale(pp_b_dat)
#  save(pp_b_dat, file=paste0('POSTPILOT_',win_mods[mm],'_stand.RData'))
# }

# predict with each model
responses = list()
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


cur_dat                = rbind(cur_dat_be) #,cur_dat_gl,cur_dat_sv)
cur_dat                = melt(cur_dat,id.vars = c('classifier'))
cur_dat_H_0            = subset(cur_dat,variable == 'random_classifier')
cur_dat_H_0$mean_auc   = cur_dat$value[cur_dat$variable == 'mean_auc']
cur_dat                = cur_dat_H_0
cur_dat$AUC_ROC        = cur_dat$value
cur_dat$value          = NULL
cur_dat$algorithm      = cur_dat$classifier
cur_dat$classifier     = cur_dat$variable
cur_dat$classifier = agk.recode.c(cur_dat$classifier,'random_classifier','random')

p = ggplot(cur_dat,aes(x=AUC_ROC, fill=classifier)) + geom_density(alpha=0.25)
#p = p + facet_grid(classifier ~ .) + ggtitle('AUC densities for estimated classifier compared to random classifier')
p = p + ggtitle('AUC of trained classifier compared to random classifier')
p = p + geom_vline(aes(xintercept = mean_auc),colour = 'green',size= 1.5)
p = p + coord_cartesian(xlim = c(0.42, 0.8)) 
p = p + theme_bw()
p = p + theme(axis.text=element_text(size=30, face = "bold"),
              axis.title=element_text(size=30,face="bold"))
p = p + theme(plot.title = element_text(size=25,face = 'bold'))
p = p + theme(legend.text = element_text(size=25))
p = p + theme(legend.title= element_text(size=25))
p = p + xlab('AUC ROC')
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
p = p + facet_grid(classifier ~ .) + ggtitle('AUC densities for elastic net classifier and random classifier')
p = p + geom_vline(aes(xintercept = as.numeric(all_mod_mean_auc)),colour = 'green',size= 1.5)
p = p + theme_bw()
p = p + theme(axis.text=element_text(size=14, face = "bold"),
              axis.title=element_text(size=20,face="bold"))
p = p + coord_cartesian(xlim = c(0.4, 0.8)) 
p = p + theme(plot.title = element_text(size=22))
p = p + xlab('AUC ROC')
print(p)



# plots also the margin value
mean_response = weighted_responses/1001
mean_response_p = rand(length(mean_response),1)-0.5
real_group = dat_match$HCPG

HCmarg = boot::inv.logit(simpleboot::one.boot(mean_response[real_group == 'HC'], median,50000)$t)
GDmarg = boot::inv.logit(simpleboot::one.boot(mean_response[real_group == 'PG'], median,50000)$t)
HCmarg_p = boot::inv.logit(simpleboot::one.boot(mean_response_p[real_group == 'HC'], median,50000)$t)
GDmarg_p = boot::inv.logit(simpleboot::one.boot(mean_response_p[real_group == 'PG'], median,50000)$t)


cur_dat_marg_cl        = data.frame(HCmarg = HCmarg,GDmarg = GDmarg,classifier = 'elastic net')
cur_dat_marg_bl        = data.frame(HCmarg = HCmarg_p,GDmarg = GDmarg_p,classifier = 'random')
cur_dat                = rbind(cur_dat_marg_cl,cur_dat_marg_bl) #rbind(cur_dat_be,cur_dat_gl,cur_dat_sv)
cur_dat                = melt(cur_dat,id.vars = c('classifier'))
cur_dat$group = cur_dat$variable 
cur_dat$variable = NULL

# density plot
p = ggplot(cur_dat,aes(x=value, fill=group)) + geom_density(aes(x=value, y=..scaled.., fill=group),kernel = 'gaussian', n = 512,alpha=0.25)
p = p + ggtitle('Probability distributions of median GD and HC margins on valid. sample test data\nfor elastic net classifier and baseline classifier over 1000 rounds of classifier estimation')
p = p + scale_x_continuous(lim = c(0, 1))
p = p + geom_vline(aes(xintercept = 0.5),colour = 'blue',size= 1.0)
p = p + theme_bw()
p = p + facet_wrap(~ classifier, scales = "free_x",nrow=3,ncol=2)
p = p + theme(strip.text.x = element_text(size = 15, face='bold'))
p = p + theme(axis.text=element_text(size=20, face = "bold"),
              axis.title=element_text(size=20,face="bold"))
p = p + theme(plot.title = element_text(size=25, face='bold'))
p = p + theme(legend.text = element_text(size=25))
p = p + theme(legend.title= element_text(size=25))
p = p + xlab('median margin of classifier [probability(Group = GD)]')
p = p + ylab('probability of median margin value')
print(p)

# print the median margin values
message('Median marginal values for GD, HC and then GD, HC under H0:')
print(median(GDmarg))
print(median(HCmarg))

print(median(GDmarg_p))
print(median(HCmarg_p))

## EXTRA NEW STUFF
## do a paired test ===============================================================
# shorten to see if we can get it with less samples
Ha_auc_short = Ha_auc[1:1000] 
improvement = c()
ct          = 0
for (cur_auc in all_aucs){
  ct = ct + 1
  improvement[ct] = sample(Ha_auc_short,size = 1) - cur_auc
}
agk.density_p(improvement,0)

## take mean of single predictions ================================================
message('p-value for mean(AUC):')
message(' ')
message(1-agk.density_p.c(all_aucs,mean(Ha_auc)))
