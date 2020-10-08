# BrazilDCI_R
Source code for Safeguarding migratory fish via strategic planning of future small hydropower in Brazil. 
Thiago B. A. Couto, Mathis L. Messager & Julian D. Olden

##Set-up
All analysis was perform in R 4.0.  
We used a project structure with dependency (package) management using renv (See [renv reference] {https://rstudio.github.io/renv/articles/renv.html}).  
When one launches this project, renv should automatically bootstrap itself, 
thereby downloading and installing the appropriate version of renv into the 
project library. After this has completed, one can then use renv::restore() 
to restore the project library locally on their machine.

The only files needed for the core of the R analysis are:   
- results/dci.gdb/networkattributes  
- results/dci.gdb/damattributes  
Both are available in the study's figshare permanent repository (add link).

For ease of use, codes were numbered in the order that they need to be run for 
reproducing the analysis.  
Despite this order, note that all codes can be run in separate sessions 
i.e. if one has run all previous code beforehand, 
they do not need to re-use the environment from previous scripts
(e.g. one can run 1.DCIAnalysis_New, clear their environment, 
and then run 2.DCI_fragmentTime_New)	

##Workflow
00.DCI_functions.R and 00.DCI_packages.R: do not need to be run manually,
they are sourced from subsequent codes. 

To proceed with the analysis, run the scripts in the following order
(run for both DCIp for potadromous and DCIi for diadromous species):  
- 0.DCI_Format: read in network and dam data (results/dci.gdb/networkattributes and results/dci.gdb/damattributes), and format them  
- 1.DCIAnalysis_New: compute DCI for current and future scenarios; for SHP-only, LHP-only, and all dams  
- 2.DCI_fragmentTime_New.R: compute DCI over time, generate Figure 1, get statistics for manuscript  
- 3.DCI_SamplingIndividualDams_NewFast.R: compute DCI for every possible dam-building scenario by basin, compute DCI loss with and without individual dams  
- 4.DCI_SamplingPareto_New.R: Run sampling based pareto-front analysis. Sample five million future dam portfolios  .  
- 5.DCI_Analysis_plot.R: plot figure 2. get statistics for manuscript.  
- 6.MigratorySpp_New.R: analyze data for migratory species. produce figure 4 and statistics.  
- 7.DCI_map.R: plot figure 3 (map)  
- 8.DCI_SensitivityAnalysis.R: run sensitivity analysis regarding dam passability of SHP vs LHP  
- ANEEL_HydroSHEDS_QAQC_summary.R: check accuracy of dam locations  
- DCI_IndividualsDams_plot.R: produce figure 5  
- DCI_Pareto_plot.R: produce figure 6  
- DCI_SensitivityAnalysis_plot: produce sensitivity analysis figure  