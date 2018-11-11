# explorative correlations
# LA, PIT variables and BDI, BIS, AUDIT, KFG, GBQ, edu, ftnd

# run import_data
# run select_study (which_study = "MRT_and_POSTPILOT")
# run group_pred_init_v6.R with everything under "what to run" set to 0

# prepping
row.names(dat_match)  = dat_match$VPPG
pit_vars              = fm$ac
pit_vars              = agk.clean.intercept.name(pit_vars)
names(pit_vars)       = paste0('ac_',names(pit_vars))
pit_vars$ac_Intercept = NULL
la_vars               = fm$laCh_LA
la_vars               = la_vars[c('gain','loss','LA')]
names(la_vars)        = paste0('laCh_',names(la_vars))


# merging
cur_merge = agk.merge.df.by.row.names(dat_match,pit_vars)
cur_merge = agk.merge.df.by.row.names(cur_merge,la_vars)
cur_merge = subset(cur_merge, HCPG == 'PG')

# selecting vars
behav_params  = cur_merge[c(grep('^ac_',names(cur_merge)),grep('^laCh_',names(cur_merge)))]
des_params    = c('BDI', 'BIS','BIS_AmBasImp','BIS_NPlanImp','BIS_Mot', 'audit', 'KFG', 'GBQ','GBQ_illus','GBQ_persi' ,'edu_years', 'smoking_ftdt')
constr_params = cur_merge[des_params]

# mean imputation for BDI
constr_params$BDI = agk.impute.mean(constr_params$BDI)

# correlation
r_sp   = corr.test(behav_params,constr_params,adjust = 'none',method = 'spearman')
r_ps   = corr.test(behav_params,constr_params,adjust = 'none',method = 'pearson')
r_ps_p = r_ps$p
r_sp_p = r_sp$p

# plot function if sig
agk.plot.cor.if.sig = function(r_ps_p,r_sp_p,behav_params, constr_params) {
  stopifnot(all(dim(r_ps_p) == dim(r_sp_p)))
  for (ii in 1:length(r_ps_p[,1])) {
    for (jj in 1:length(r_ps_p[1,])) {
      cur_p_ps = r_ps_p[ii,jj]
      cur_p_sp = r_sp_p[ii,jj]
      if (cur_p_sp < 0.05 & cur_p_ps < 0.05) {
        plot(behav_params[[ii]],constr_params[[jj]],
             ylab=names(constr_params)[jj],
             xlab=names(behav_params)[ii],
             xlim=c(-5,5))
      }
    }
  }
}

agk.plot.cor.if.sig(r_sp$p,r_ps$p,behav_params,constr_params)
