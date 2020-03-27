library(rprojroot)
library(data.table)
library(fst)
library(plyr)
library(ggplot2)

# Import network and dams dataset (Mathis folder structure)
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
outdir_sensitivity = file.path(resdir, 'outsensitivity_passscenarios')

#Import sensitivity analysis files
DCIsensitivity <-  Map(read.fst,
                     file.path(outdir_sensitivity,
                               list.files(path=outdir_sensitivity, pattern="DCIsensitivity_.*[.]fst"))) %>%
  do.call(rbind, .) %>%
  setDT %>%
  .[, DCI := as.integer(round(DCI))]

#Plot distribution of variables to check for odd values
summary(DCIsensitivity)

#Identify basins that are currently free of hydropower, currently regulated, and those with both SHPs and LHPs in the future
DCIsensitivity[, `:=`(prevfree = +(DAMBAS_ID08ext %in% 
                                     DCIsensitivity[scenario=='Allcurrent' & DCI==100, 
                                                    DAMBAS_ID08ext]),
                      bothtypes = +(DAMBAS_ID08ext %in% 
                                      DCIsensitivity[(scenario=='SHPfuture' & DCI!=100) | (scenario=='LHPfuture' & DCI!=100), 
                                                     DAMBAS_ID08ext]))]
  
#Filter results 
               
      
#Compute average brazil-wide DCI for 
DCIsensitivity_mean <- DCIsensitivty[mean(DCI), by=.(passpar, scenario)]



################################################################################################################################
################################################################################################################################

## Calculate differences
DCILoss_All <- resuDCI$All_fut - resuDCI$All_curr
DCILoss_SHP <- resuDCI$SHP_fut - resuDCI$SHP_curr
DCILoss_LHP <- resuDCI$LHP_fut - resuDCI$LHP_curr

DCILoss_Both_All <- resuDCI_Both$All_fut - resuDCI_Both$All_curr
DCILoss_Both_SHP <- resuDCI_Both$SHP_fut - resuDCI_Both$SHP_curr
DCILoss_Both_LHP <- resuDCI_Both$LHP_fut - resuDCI_Both$LHP_curr

DCILoss_Free_All <- resuDCI_Free$All_fut - resuDCI_Free$All_curr
DCILoss_Free_SHP <- resuDCI_Free$SHP_fut - resuDCI_Free$SHP_curr
DCILoss_Free_LHP <- resuDCI_Free$LHP_fut - resuDCI_Free$LHP_curr

DCILoss_Regul_All <- resuDCI_Regul$All_fut - resuDCI_Regul$All_curr
DCILoss_Regul_SHP <- resuDCI_Regul$SHP_fut - resuDCI_Regul$SHP_curr
DCILoss_Regul_LHP <- resuDCI_Regul$LHP_fut - resuDCI_Regul$LHP_curr


## Calculate % change (new - old / old * 100)
PercLoss_All <- DCILoss_All/resuDCI$All_curr * 100
PercLoss_SHP <- DCILoss_SHP/resuDCI$SHP_curr * 100
PercLoss_LHP <- DCILoss_LHP/resuDCI$LHP_curr * 100

PercLoss_Both_All <- DCILoss_Both_All/resuDCI_Both$All_curr * 100
PercLoss_Both_SHP <- DCILoss_Both_SHP/resuDCI_Both$SHP_curr * 100
PercLoss_Both_LHP <- DCILoss_Both_LHP/resuDCI_Both$LHP_curr * 100

PercLoss_Free_All <- DCILoss_Free_All/resuDCI_Free$All_curr * 100
PercLoss_Free_SHP <- DCILoss_Free_SHP/resuDCI_Free$SHP_curr * 100
PercLoss_Free_LHP <- DCILoss_Free_LHP/resuDCI_Free$LHP_curr * 100

PercLoss_Regul_All <- DCILoss_Regul_All/resuDCI_Regul$All_curr * 100
PercLoss_Regul_SHP <- DCILoss_Regul_SHP/resuDCI_Regul$SHP_curr * 100
PercLoss_Regul_LHP <- DCILoss_Regul_LHP/resuDCI_Regul$LHP_curr * 100




points(x = treatment[1] * 1, y = mean(PercLoss_Free_SHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 2, y = mean(PercLoss_Free_LHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 3, y = mean(PercLoss_Free_All), pch = "-", col = "black", cex = 13.0)