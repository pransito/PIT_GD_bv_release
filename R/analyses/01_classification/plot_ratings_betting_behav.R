## plot ratings, collect info on betting behavior =============================
## SETTINGS
cur_control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE,
                          check.conv.grad="ignore",
                          check.conv.singular="ignore",
                          check.conv.hess="ignore",
                          optCtrl=list(optimizer = "nloptwrap",maxfun=500))

## Plots means and CIs of ratings variables (stim)=============================
# means and CIs
des_vars    = c('valence','arousal','dominance',"imageRating1s","imageRating2s","imageRating3s","imageRating4s")
cur_map     = c('Valence', 'Arousal','Dominance', "Elicits_Craving",
                "Representative_for_Gambles","Representative_for_Negative","Representative_for_Positive")
des_vars_ra = des_vars

cfint_list = list()
meansra    = list()
for (ii in 1:length(des_vars)) {
  
  cur_form_str = paste(des_vars[ii],'~ subject + cat')
  cur_form     = as.formula(cur_form_str)
  cfint        = aggregate(cur_form,data = data_pdt,FUN=mean)
  # unit test no missings
  if (any(is.na(cfint[[3]]))) {
    stop('Missings!')
  }
  # pack for predictions
  cur_form      = paste(gsub('subject + ','',cur_form_str,fixed=T),'|subject')
  cur_lml       = lmList(as.formula(cur_form),data=data_pdt,pool = F)
  cur_lml       = coef(cur_lml)
  meansra[[ii]] = cur_lml
  
  if (plot_ratings) {
    # more for plotting
    cfint$HCPG = agk.recode.c(cfint$subject,dat_match$VPPG,dat_match$HCPG)
    cur_var    = names(cfint)[3]
    cur_form   = as.formula(paste(cur_var,'~ HCPG + cat'))
    cfint$HCPG = agk.recode(cfint$HCPG,c('PG'),c('GD'))
    cfint$HCPG = factor(cfint$HCPG, levels = c('HC','GD'))
    cfint      = aggregate(cur_form,data = cfint,FUN=agk.boot.ci,
                           cur_fun = mean,lower=0.025, upper=0.975,R=5000)
    
    cfint = as.data.frame(cbind(cfint[c('HCPG','cat')],cfint[[3]]))
    names(cfint)[c(3:5)] = c('mean','lower','upper')
    cfint$variable = agk.recode(des_vars[ii],des_vars,cur_map)
    cfint_list[[ii]] = cfint
  }
}

if (plot_ratings) {
  cfint = cfint_list[[1]]
  for(ii in 2:length(cfint_list)) {
    cfint = rbind(cfint,cfint_list[[ii]])
  }
  
  # change the order
  cfint$variable = factor(cfint$variable,
                          levels =  c('Valence', 'Arousal', 'Dominance',
                                      'Elicits_Craving','Representative_for_Gambles','Representative_for_Negative','Representative_for_Positive'),
                          labels = c('Valence', 'Arousal', 'Dominance',
                                     'Elicits_Craving','Representative_for_Gambles','Representative_for_Negative','Representative_for_Positive'))
  
  # plotting (faceting by category)
  mRat  = ggplot(cfint, aes(HCPG, mean,fill=variable))
  mRat  = mRat + labs(x='Group', y=paste('Mean (',0.95*100,'% CI, bootstrapped)'))
  mRat  = mRat + ggtitle("Ratings of different categories of cues")
  mRat  = mRat + geom_bar(position="dodge", stat="identity")
  dodge = position_dodge(width=0.9)
  mRat  = mRat + geom_bar(position=dodge, stat="identity")
  mRat  = mRat + geom_errorbar(aes(ymin = lower, ymax = upper), position=dodge, width=0.25) 
  mRat  = mRat + facet_grid(cat ~ .)
  mRatg = mRat + theme_bw()
  print(mRatg)
}

## get info which gambles they do by SOGS
des_vars = c('SOGS_1[a-z]')
des_vars = names(dat_match)[grep(des_vars,names(dat_match))]
cur_map = c('cards', 'bet_animals','bet_sports', "dice","casinos",
            "lotteries","bingo",'stock_market','slot_machines','played_sports')

cfint_list = list()
for (ii in 1:length(des_vars)) {
  x1      = dat_match[dat_match$HCPG == 'HC',des_vars[ii]]
  x2      = dat_match[dat_match$HCPG == 'PG',des_vars[ii]]
  
  # unit test no missings
  if (any(is.na(x1)) | any(is.na(x2))) {
    stop('Missings!')
  }
  
  cfint   = t.test(x1,x2,'less')
  cur_res = list()
  cur_res$meanHC = mean(x1)
  cur_res$meanPG = mean(x2)
  cur_res$p      = round(cfint$p.value,digits = 3)
  
  cfint_list[[cur_map[ii]]] = cur_res
}

cfint = as.data.frame(cfint_list[[1]])
cfint$gamble = names(cfint_list)[1]
for(ii in 2:length(cfint_list)) {
  cur_res        = as.data.frame(cfint_list[[ii]])
  cur_res$gamble = names(cfint_list)[ii]
  cfint = rbind(cfint,cur_res)
}
cfint = cfint[c(4,1,2,3)]

message('Kind of gambles HC and GD play, with p-values')
print(cfint)

# pack ratings for prediction
all_des_vars = c(des_vars_ra)
crdf        = agk.clean.intercept.name(meansra[[1]])
names(crdf) = paste0(names(crdf),'_',all_des_vars[1])
for (ll in 2:length(meansra)) {
  cur_df        = agk.clean.intercept.name(meansra[[ll]])
  names(cur_df) = paste0(names(cur_df),'_',all_des_vars[ll])
  crdf          = data.frame(crdf,cur_df)
}
cr_agg_ra         = crdf
cr_agg_ra$subject = row.names(cr_agg_ra)