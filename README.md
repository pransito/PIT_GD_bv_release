R Code for the PIT Gambling Disorder paper (behavioral data)
-----------------------------------------------------------

Code for publication of "Cue-induced effects on decision-making distinguish subjects with gambling disorder from healthy controls"

see here: https://doi.org/10.1111/adb.12841; https://doi.org/10.1101/564781

by: Alexander Genauck et al. (2019)


How to use
----------

Fork the whole repository or download the zip and put it into some working directory.
The .RData file has all the data and functions you will need. Especially data_pdt (choice data) and dat_match (questionnaire, demographic data). You will need R (https://cran.r-project.org/bin/windows/base/) and Rstudio (https://www.rstudio.com/products/rstudio/#Desktop) for this. Both are freely available software.

Version of R this project has been developed with: 3.5.0 (2018-04-23)
Version of RStudio this project has been developed with: 1.1.453

1. Set-up
Run the "R/select_study.R" script (press "source" button in upper right corner). It selects the Cohort (POSTPILOT_HCPG) and initializes all the analyses. Do this before anything else. fm (a list holding all the estimated parameters from all models for each subject) will be created. 

Careful: The script "R/analyses/01_classification/group_pred_loop_v7.R" (called by "R/select_study.R") installs and loads many R packages (see in the beginning of the script). They are all useful and should not hurt your installation of R or bother your other code. However, revise the list before you run the code and decide
if you would like to continue.

2. Get the results
If you want to see results (p-values, classifier performance, graphs) you set everything to FALSE in "WHAT TO RUN" section in "R/analyses/01_classification/classifier_settings.R", except "do_report_no_added_feat" to TRUE. Set runs to "1008" to see the results reported in the paper, then run the script: "R/analyses/01_classification/reporting_v7.R"

3. Run the algorithm
The machine learning part is started with the script "R/analyses/01_classification/group_pred_loop_v7.R" Before running it, you may adjust the settings in: "R/analyses/01_classification/classifier_settings.R" [set "runs" to a low number like 10 for it is an intense script which takes a long time to run]. Set everything under "WHAT TO RUN" to TRUE except "do_report_no_added_feat" to FALSE, unless you want to get a report on your runs (then runs should be at least 50-100 to get sensible results and you would have to run "R/analyses/01_classification/reporting_v7.R" again). 

4. Get the performance on external data set
Producing p-values, AUC-values, graphs for applying the classifier on the validation data, first run "R/select_study.R" with 'which_study = "MRI"' (line 34) to switch to the other study. Then run: "R/analyses/01_classification/apply_PIT_GD_behav_to_MRI_sample_v2.R"

5. Classical group analyses
The classical group analysis part using hierarchical regression (lme4) is run with the
"/02_classical_group_analyses/glmer_accRate_la_cat_v3.R script"; Before running it, make sure to run "R/select_study.R" first with which_study set to "POSTPILOT_HCPG". It is set such that the glmer models are run (10 to 20 minutes) but the confidence intervals for the fixed effects of beta_gain, beta_loss and loss aversion per group are just loaded; running it takes about an hour or longer; 

6. Ratings graph
Making the ratings graph. It is automatically produced when running "select_study.R" with which_study = "POSTPILOT_HCPG". 

7. Analysis of ratings data
Ratings: statistical tests. Run the script "R/analyses/03_image_adequacy/ratings_analysis_for_paper.R". Check the instructions at the top of the script. 


Enjoy! Any questions, or bugs please send by raising an issue in this GitHub repository.
