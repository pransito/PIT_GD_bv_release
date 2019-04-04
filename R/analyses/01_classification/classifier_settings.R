## PREAMBLE ===================================================================
# script what to run, already set like in paper;
# tip: set runs to a smaller number than 1008 if you do not want to wait hours
# number of runs to get the CV results distribution, >=1000 recommended
# runs also will be the name of the results folder
# available results in results folder and their contents described further down
# for a result (try 10 for starters)

## WHAT TO RUN ================================================================
# just the behavioral parameter sets
outer_cv_noaddfeat      = F # with outer CV, getting generalization error, Ha
noout_cv_noaddfeat      = F # no outer CV, get complete model on whole sample

# control model
outer_cv_c_model        = F # baseline model for classif.: predict with covariate

# what to report
do_report_no_added_feat = T # set to 1 to see the report of the behavioral classifier

# number of runs to get the CV results distribution, >=1000 recommended
# 1008: for reporting of results as in paper
runs = 1008

# advanced settings (do not change) ===========================================
report_CV_p = T

# no other features, only behavior
# master add cue reactivity: peripheral physiology or MRI
if (outer_cv_noaddfeat == T | noout_cv_noaddfeat == T | do_report_no_added_feat == T) {
  add_cr_pp_ma = F
} else {
  add_cr_pp_ma = T
}

# master add cue reactivity: ratings
# should never be done, cause ratings are post-experiment
add_cr_ra_ma         = F


