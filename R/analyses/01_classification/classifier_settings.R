## PREAMBLE ===================================================================
# script what to run, already set like in paper;
# tip: set runs to a smaller number than 1010 if you do not want to wait hours
# number of runs to get the CV results distribution, >=1000 recommended
# runs also will be the name of the results folder
# available results in results folder and their contents described further down
# for a result (try 10 for starters)

## WHAT TO RUN ================================================================
# just the behavioral parameter sets
outer_cv_noaddfeat      = T # with outer CV, getting generalization error, Ha
noout_cv_noaddfeat      = T # no outer CV, get complete model on whole sample

# control model
outer_cv_c_model        = T # control model/null-model for classification; predict with covariate

# what to report
do_report               = F
do_report_no_added_feat = F

# how many runs?
# set to 1010: behav predictions (PIT GD behav paper, for reporting) 
runs                    = 10

# advanced settings for other studies =========================================
# [cannot be used in PIT GD behav release] 
if (which_study == 'MRI') {
  # Any reporting of p-values against null? Set to F if you do that in a separate script.
  report_CV_p = F
} else {
  # Any reporting of p-values against null? Set to F if you do that in a separate script.
  report_CV_p = T
}

# no other features, only behavior
# master add cue reactivity: peripheral physiology or MRI
add_cr_pp_ma         = F
# master add cue reactivity: ratings
# should never be done, cause ratings are post-experiment
add_cr_ra_ma         = F

