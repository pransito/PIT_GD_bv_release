R_PDT
-----

Code for publication of "Cue-induced effects on decision-making distinguish gambling disorder subjects from healthy controls"
by: Alexander Genauck (2018)


How to use
----------

Fork the whole repository or download the zip and put it into some working directory.
The .RData file has all the data and functions you will need. Especially data_pdt (choice data) and dat_match (questionnaire, demographic data). You will need R (https://cran.r-project.org/bin/windows/base/) and Rstudio (https://www.rstudio.com/products/rstudio/#Desktop) for this. Both are freely available software.

1)
Run the "R/select_study.R" script. It selects the Cohort (POSTPILOT_HCPG) and initializes all the analyses. Do this before anything else. fm (a list holding all the estimated parameters from all models for each subject) will be created. 

Careful: The script "R/analyses/01_classification/severity_pred_loop_v7.R" (called by "R/select_study.R") installs and loads many R packages (see in the beginning of the script). They are all useful and should not hurt your installation of R or bother your other code. However, revise the list before you run the code and decide
if you would like to continue.

2)
The machine learning part is started with the script "R/analyses/01_classification/severity_pred_loop_v7.R" Before running it, you may adjust the settings in: "R/analyses/01_classification/classifier_settings.R" [set "runs" to a low number like 10 for it is an intense script which takes a long time to run]  If you want to see results you set everything to FALSE, except the reporting section to TRUE (all of them) Set runs to "1010" to see the results reported in the paper. 

3)
Producing p-values for applying the classifier on the validation data you use "R/analyses/01_classification/apply_PIT_GD_behav_to_MRI_sample_v2.R"; Before running the script, run "R/select_study.R" with 'which_study = "MRI"' 

4)
The classical group analysis part using hierarchical regression (lme4) is  with the
"/02_classical_group_analyses/glmer_accRate_la_cat_v3.R script"; Before running it, make sure to run "R/select_study.R" first with which_study set to "POSTPILOT_HCPG". It is set such that the glmer models are run (10 to 20 minutes) but the confidence intervals for the fixed effects of beta_gain, beta_loss and loss aversion per group are just loaded; running it takes about an hour or longer; 

5)
Making the ratings graph. It is automatically produced when running "severity_pred_loop_v7.R". 

6)
Ratings: statistical tests. Run the script "R/analyses/03_image_adequacy/ratings_analysis_for_paper.R". Check the instructions at the top of the script. 

7)
Explorative correlations between loss aversion/PIT parameters and questionnaire data, is  with the script: "/02_classical_group_analyses/correlations_PIT_LA_questionnaires.R" Check the instructions at the top of the script.


Enjoy! Any questions, or bugs please send by raising an issue in this GitHub repository.