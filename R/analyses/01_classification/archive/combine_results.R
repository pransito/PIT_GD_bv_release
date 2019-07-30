# put together results

# POSTPILOT_HCPG_predGrp1_rounds_noo_noaddfeat.RData
e_300 = new.env()
load('POSTPILOT_HCPG_predGrp1_rounds_noo_noaddfeat_300.RData',envir = e_300)

e_703 = new.env()
load('POSTPILOT_HCPG_predGrp1_rounds_noo_noaddfeat_703.RData',envir = e_703)

cur_mod_sel_nooCV = c(e_300$cur_mod_sel_nooCV,e_703$cur_mod_sel_nooCV)
CVnoo_res_list = c(e_300$CVnoo_res_list,e_703$CVnoo_res_list)
des_seed = e_703$des_seed
fm = e_703$fm
list_winning_model_c_nooCV = c(e_300$list_winning_model_c_nooCV,e_703$list_winning_model_c_nooCV)
list_winning_model_l_nooCV = c(e_300$list_winning_model_l_nooCV,e_703$list_winning_model_l_nooCV)

save(file='POSTPILOT_HCPG_predGrp1_rounds_noo_noaddfeat.RData',list = c('cur_mod_sel_nooCV','CVnoo_res_list',
                                                                        'des_seed','fm','list_winning_model_c_nooCV','list_winning_model_l_nooCV'))


# POSTPILOT_HCPG_predGrp1_rounds_wio_conmod.RData
e_300 = new.env()
load('POSTPILOT_HCPG_predGrp1_rounds_wio_conmod_300.RData',envir = e_300)

e_703 = new.env()
load('POSTPILOT_HCPG_predGrp1_rounds_wio_conmod_703.RData',envir = e_703)

CVcm_res_list = c(e_300$CVcm_res_list,e_703$CVcm_res_list)
fm = e_703$fm
des_seed = e_703$des_seed

save(file='POSTPILOT_HCPG_predGrp1_rounds_wio_conmod.RData',list = c('CVcm_res_list','fm',
                                                                        'des_seed'))

# POSTPILOT_HCPG_predGrp1_rounds_wio_noaddfeat.RData
e_300 = new.env()
load('POSTPILOT_HCPG_predGrp1_rounds_wio_noaddfeat_300.RData',envir = e_300)

e_703 = new.env()
load('POSTPILOT_HCPG_predGrp1_rounds_wio_noaddfeat_703.RData',envir = e_703)

CV_res_list = c(e_300$CV_res_list,e_703$CV_res_list)
fm = e_703$fm
des_seed = e_703$des_seed

save(file='POSTPILOT_HCPG_predGrp1_rounds_wio_noaddfeat.RData',list = c('CV_res_list','fm',
                                                                     'des_seed'))



