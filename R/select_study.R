## PREAMBLE ===================================================================
# script that selects from data_pdt and dat_match the data of desired study
# the study (cohort) you need
# run data_import.R before or load the .RData file (foreign run)
# will output data_pdt, dat_match to be used for further analysis

# set to T, if you have this script from GitHub
FOREIGN_RUN = T

# PREPARATION FOR FOREIGN RUN =================================================
if (FOREIGN_RUN) {
  # find our the roo directory
  rm(list=ls())
  root_wd  = paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/analyses/')
  setwd(root_wd)
  load('.RData')
}

# get the paths
res      = agk.get.working.location()
base_gd  = res$base_gd
path_ghb = res$path_ghb

# foreign run
path_res_classif = file.path(root_wd,'01_classification')

# get the original data_pdt from data_import.R
data_pdt     = data_pdt_bcp
data_pdt_inv = data_pdt

## PARAMETER SETTINGS =========================================================
# which study to look at (Cohorts)? ===========================================
which_study = "POSTPILOT_HCPG"
#which_study = "MRI"


# add novelity of cue as factor
# will control all single subject models for novelty of cue (duplicate or not) 
add_novelty = F

# cut gambling cues from 67 to 45?
cut_gambling_cues = F

## PREPARATIONS ===============================================================
if (physio_sum_fun == 'mean') {
  data_pdt$corr = data_pdt$corr_auc
  data_pdt$eda  = data_pdt$eda_auc
  data_pdt$zygo = data_pdt$zygo_auc
} else if (physio_sum_fun == 'max') {
  data_pdt$corr = data_pdt$corr_max
  data_pdt$eda  = data_pdt$eda_max
  data_pdt$zygo = data_pdt$zygo_max
} else if (physio_sum_fun == 'median') {
  data_pdt$corr = data_pdt$corr_median
  data_pdt$eda  = data_pdt$eda_median
  data_pdt$zygo = data_pdt$zygo_median
} else if (physio_sum_fun == 'all') {
  # do nothing
  # here all of the stuff will be used
} else {
  stop('No proper physio_sum_fun provided.')
}

## FUNCTIONS ==================================================================
# summary function
cur_summary_function = function(x) median(x, na.rm=TRUE)
# needs the f.difftest function from the import_data file
f = function(x) {
  tmp <- t.test(x)
  return(as.matrix(c(tmp$p.value,mean(x,na.rm = T))))
}

agk.load.ifnot.install <- function(package_name){
  if(require(package_name,character.only = T,quietly = T)){
    print(paste (package_name,"is loaded correctly"))
  } else {
    print(paste("trying to install", package_name))
    install.packages(pkgs = c(package_name))
    if(require(package_name,character.only = T)){
      print(paste(package_name,"installed and loaded"))
    } else {
      stop(paste("could not install",package_name))
    }
  }
}

## PACKAGES ===================================================================
agk.load.ifnot.install('reshape2')

## SUBSETTING DATA ============================================================
data_pdt$Cohort[data_pdt$Cohort ==  "PhysioIAPS"] = "sanity"
data_pdt$Cohort[is.na(data_pdt$Cohort)]           = "Pretest"

if (which_study == "Prestudy") {
  data_pdt = subset(data_pdt, Cohort == "Pretest" | Cohort == "PhysioPilot")
} else if (which_study == "POSTPILOT_HCPG" | which_study == "POSTPILOT_HC" | 
           which_study == "POSTPILOT_PG" | which_study == "POSTPILOT_PGxGENDER") {
  data_pdt = subset(data_pdt, Cohort == "POSTPILOT" | Cohort == "PGPilot")
} else if (which_study == "MRI" | which_study == "MRI_HC" | which_study == "MRI_PG"| which_study == "MRI_LB") {
  data_pdt = subset(data_pdt,Cohort == "MRI")
} else if (which_study == "PhysioPilot") {
  data_pdt = subset(data_pdt,Cohort == "PhysioPilot")
} else if (which_study == "MRI_LB") {
  includelist=read.table("E:/MATLAB/info_mri_selection.csv")
  dat_match_MRI_only=subset(dat_match, dat_match$Cohort=="MRI")
  dat_match_MRI_only=subset(dat_match_MRI_only, dat_match_MRI_only$Einschluss==1)
  data_pdt= subset(data_pdt,data_pdt$subject %in% dat_match_MRI_only$VPPG)
  data_pdt= subset(data_pdt,data_MRI_only$subject %in% includelist$V1)
} else if (which_study == "sanity") {
  data_pdt = subset(data_pdt,Cohort == "sanity")
} else if (which_study == "MRI_and_POSTPILOT") {
  data_pdt = subset(data_pdt,Cohort == "POSTPILOT" | Cohort == "PGPilot" | Cohort == "MRI")
} else if (which_study == "TEST") {
  data_pdt = subset(data_pdt, Cohort == "TEST")
} else {
  stop('No valid cohort selected with var which_study!')
}

# get a HCPG variable
data_pdt$HCPG   = as.factor(as.character(data_pdt$HCPG))

# only one group, i.e. HC or PG?
if ((length(grep(which_study, pattern = "HC")) != 0) & (length(grep(which_study, pattern = "PG")) == 0)) {
  data_pdt = subset(data_pdt,HCPG == "HC")
} else if ((length(grep(which_study, pattern = "PG")) != 0) & (length(grep(which_study, pattern = "HC")) == 0)) {
  data_pdt = subset(data_pdt,HCPG == "PG")
}

if ((length(grep(which_study, pattern = "HC")) != 0) & (length(grep(which_study, pattern = "PG")) == 0)) {
  data_pdt = subset(data_pdt,HCPG == "HC")
} else if ((length(grep(which_study, pattern = "PG")) != 0) & (length(grep(which_study, pattern = "HC")) == 0)) {
  data_pdt = subset(data_pdt,HCPG == "PG")
}

# CATEGORY LABELS =============================================================
# Main effect of the final experimental categories: gam, pos, neg, neu_aw
data_pdt_finCat = data_pdt
if (which_study == "Prestudy") {
  data_pdt_finCat$cat = agk.recode.c(as.character(data_pdt_finCat$cat),c("1","2","3"),c("1","2","3"))
  data_pdt_finCat$cat = factor(as.numeric(as.character(data_pdt_finCat$cat)),levels = c(1,2,3),
                               labels = c('gambling','negative', 'positive')) ## CAREFUL: I TOOK OUT ALL NEUTRAL PICTURES HERE!
} else if (which_study == "sanity") {
  data_pdt_finCat$cat = factor(as.numeric(as.character(data_pdt_finCat$cat)),levels = c(6,2,3,7,8),
                               labels = c('neutral_IAPS','negative_VPPG', 'positive_VPPG','negative_IAPS','positive_IAPS'))
} else if (which_study == "POSTPILOT_PG" | which_study == "POSTPILOT_PGxGENDER" | 
           which_study == "POSTPILOT_HC" | which_study == "MRI_HC" | 
           which_study == "MRI_PG" | which_study == "POSTPILOT_HCPG" | 
           which_study == "MRI" | which_study == "MRI_and_POSTPILOT" |
           which_study == "TEST") {
  data_pdt_finCat$cat = factor(as.numeric(as.character(data_pdt_finCat$cat)),levels = c(6,1,2,3),
                               labels = c('neutral','gambling','negative', 'positive'))
}

if(sum(is.na(data_pdt$cat))) {
  stop('There are NAs in the data_pdt$cat variable!')
}

## functions ==================================================================
getmode <- function(v) {
  # mode function
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## VARIABLE TRANSFORMATIONS ===================================================
data_pdt_finCat$valence_log       = get.log(data_pdt_finCat$valence)
data_pdt_finCat$imageRating2s_log = get.log.base(data_pdt_finCat$imageRating2s,10)
data_pdt_finCat$imageRating1s_log = get.log(data_pdt_finCat$imageRating1s)
data_pdt_finCat$imageRating4s_log = get.log(data_pdt_finCat$imageRating4s)
data_pdt_finCat$imageRating3s_log = get.log(data_pdt_finCat$imageRating3s)


## GET CAT LABELS AND TRANSFORMATIONS INTO DATA_PDT ===========================
# DO THIS BFORE STARTING ANY ANALYSIS OF BEHAVIORAL PDT TASK DATA
data_pdt         = data_pdt_finCat
data_pdt$subject = droplevels(data_pdt$subject)

## UNIT TEST GROUP VAR ========================================================
# test of group variable
if (which_study != 'sanity' & which_study != 'TEST') {
  data_pdt = data_pdt[!(data_pdt$HCPG != "PG" & data_pdt$HCPG != "HC"),]
  if (length(levels(data_pdt$HCPG))>2) {
    stop("WOULD NEED TO DROP SUBS DUE TO NO GROUP INFO!")
  }
}

# SUBSET DAT_MATCH ============================================================
# also select dat_match
dat_match = dat_match_bcp
if (which_study == "POSTPILOT_HCPG" | which_study == "POSTPILOT_HC" | 
    which_study == "POSTPILOT_PG" | which_study == "POSTPILOT_PGxGENDER") {
  dat_match = subset(dat_match, Cohort == "POSTPILOT" | Cohort == "PGPilot")
}
# or better yet, align
dat_match = dat_match[dat_match$VPPG %in% data_pdt$subject,]

## REPORTING ON MISSINGS AND THEN DROPPING ALL MISSINGS =======================
missing_trials = xtabs(is.na(data_pdt$accept_reject) ~ data_pdt$subject + data_pdt$HCPG)
missing_trials = melt(missing_trials)
names(missing_trials) = c('subject','HCPG','num_missing')
print(summary(lm(num_missing ~ HCPG,missing_trials)))

data_pdt = data_pdt[!is.na(data_pdt$accept_reject),]

## ADD CATEGORY VARIABLES FOR laCh ============================================
all_subs = unique(data_pdt$subject)
enh_dpdt = list()

for (ss in 1:length(all_subs)) {
  # get data and make model matrix
  cur_dat = data_pdt[data_pdt$subject == all_subs[ss],]
  cur_mm  = model.matrix.lm(accept_reject ~ (gain + loss)*cat - cat,data = cur_dat)
  cur_mm  = data.frame(cur_mm)
  
  # dropping some unneeded columns
  undes = c('X.Intercept.','gain','loss')
  cur_mm  = cur_mm[-which(names(cur_mm) %in% undes)]
  
  # attaching
  cur_dat        = data.frame(cur_dat,cur_mm)
  enh_dpdt[[ss]] = cur_dat 
}

# making new data_pdt
data_pdt = enh_dpdt[[1]]
for (ss in 2:length(enh_dpdt)) {
  data_pdt = rbind(data_pdt,enh_dpdt[[ss]])
}

## MRI DATA LOADING ===========================================================
# load the MRI extracts
# excluding subjects because of missing in pp
# prepping data frames for pp
if (which_study == 'MRI') {
  cr_agg_pp        = cr_agg_pp_r_MRI
}

## SAVE TO WORKSPACE ==========================================================
# saving this result
data_pdt_bcp_study_selected  = data_pdt
dat_match_bcp_study_selected = dat_match

## initialize for all analyses ================================================
## initialization settings [DEFAULT, DO NOT CHANGE] ===========================
# just the behavioral parameter sets
outer_cv_noaddfeat      = 0 # with outer CV, getting generalization error, Ha
noout_cv_noaddfeat      = 0 # no outer CV, get complete model on whole sample

# behavior plus peripheral-physiological stuff
outer_cv_wiaddfeat      = 0 # adding physio, Ha
noout_cv_wiaddfeat      = 0 # adding physio, get complete model

# only peripheral-physiological / MRI / rating (all saved under "phys")
outer_cv_addfeaton      = 0 # Ha only, i.e. physio/MRI  
noout_cv_addfeaton      = 0 # to get the complete model 

# control model
outer_cv_c_model        = 0 # control model/null-model for classification; predict with covariate
# not needed for MRI case (p-value comp in dfferent script, using random classification)

# what to report
do_report                 = 0
do_report_no_added_feat   = 0
do_report_with_added_feat = 0
do_report_feat_only       = 0

# Any reporting of p-values against null? Set to F if you do that in a separate script.
report_CV_p = T

# we never deal with physio or fmri data
add_cr_pp_ma = F

# master add cue reactivity: ratings
# should never be done, cause ratings are post-experiment
add_cr_ra_ma         = F

# run the initializations
setwd('..')
setwd('analyses/01_classification/')
init_run = T
source('group_pred_loop_v7.R')
init_run     = F
init_done    = T
