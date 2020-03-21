#############         Script 3 to run DCI analysis - September 2019                ###############################
#############        Based on distances calculated in Arcgis - PythonCode         ##############################
#############              Manuscript in preparation - PNAS                       ##############################

################################################################################################################
###########       FRAGMENTATION BY PLANNED PROJECTS - BRAZIL      #############################################
###########                         (INDIVIDUALLY)                ############################################
###########       RUN DCI for all possible future portifolios     ###########################################
###########            (Assign it back to each planned dam)       ##########################################
###########################################################################################################
#If run into memoery errors: see https://stackoverflow.com/questions/51248293/error-vector-memory-exhausted-limit-reached-r-3-5-0-macos


## Choose what type of DCI function to run (DCIp or DCIi)
DCIfunc <- DCIp_opti3
# DCIfunc <- DCId_opti

## Packages
require(fst)
require(tictoc)
require(plyr)
require(bigstatsr)
require(parallel)
require(doParallel)
require(Rcpp)
require(rccpcomb)
require(devtools)
devtools::source_gist("https://gist.github.com/r2evans/e5531cbab8cf421d14ed", filename = "lazyExpandGrid.R") #Get expand.grid version that won't run out of memory when > 25 dams
require(ggplot2)
library(tcltk)

# # Import network and dams dataset (Mathis folder structure)
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
dcigdb <- file.path(resdir, 'dci.gdb')

## Import network and dams dataset (alternative)
# rootdir <- find_root(has_dir("PythonOutputs"))
# datadir <- file.path(rootdir, "PythonOutputs")
# dcigdb <- file.path(datadir, 'dci.gdb')

#Import formatted data
DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% setDT
NetworkBRAZIL <- read.fst(file.path(resdir, 'NetworkBRAZIL.fst')) %>% setDT


## Remove the dams that have NAs from the dataset
DamAttributesCrude <- DamAttributesCrude[!is.na(DamAttributesCrude$HYBAS_ID08),]


#------ FORMAT DATA ------
#Estimate number of iterations required to process all dams
permuttable <- setorder(DamAttributes[ESTAGIO_1 == "Planned", 
                                      list(ndams = .N, 
                                           npermut =  do.call(lazyExpandGrid, rep(list(c(0.1, 1)), .N))$n), 
                                      by=DAMBAS_ID08ext], 'ndams') %>%
  .[,list(.N, sum=sum(npermut)), by=ndams] %>%
  .[, list(ndams = ndams, nbasins=N, permut=sum, cumpermut = cumsum(sum))]


ggplot(permuttable, aes(x=ndams, y=cumpermut)) + 
  geom_bar(stat='identity') +
  scale_y_log10()

## Create a vector with unique basin IDs
basinList <- DamAttributes[, .N, by=DAMBAS_ID08ext] %>% 
  setorder(N) %>% 
  .[,DAMBAS_ID08ext]

#Remove the basins that will have more than X future dams if needed
# X <- 10
# basinList <- DamAttributes[ESTAGIO_1 == "Planned", .N, by=DAMBAS_ID08ext][N<X, DAMBAS_ID08ext]

#### Parameters to run permutations-based DCI ####
#minimum number of dams in basin at which the function starts sampling permutation scenarios (for speed/memory sake)
#22 could be handled with 16GB of RAM
minlazy <- 20 
#number of seed scenarios to sample from lazy permutation. The actual number of scenario pairs that will actually 
#be evluated will often be at least >5-10 times more than that number
nsamples <- 10000

############## LAUNCH ANALYSIS ####################
tic("total")
#Initialize parallel analysis
cl <- parallel::makeCluster(bigstatsr::nb_cores()) #make cluster based on recommended number of cores
on.exit(stopCluster(cl))
doParallel::registerDoParallel(cl)

#60808065901 has 17 dams
DCIscens <- foreach(j=basinList, ## Loop over basins
                    .packages = c("data.table", "Rcpp", "RcppAlgos", "magrittr", 'plyr','tcltk')) %dopar% {
                      
                      ## filter attributes of the basin j         
                      NetX <- NetworkBRAZIL[HYBAS_ID08ext == j, ]
                      DamX <- DamAttributes[DAMBAS_ID08ext == j, ]       
                      DAMIDvec <- as.character(DamX[ESTAGIO_1 == "Planned", DAMID])
                      
                      #Create progress bar for long processes
                      if (length(DAMIDvec) > minlazy) {
                        mypb <- tkProgressBar(title = "R progress bar", label = "",
                                              min = 0, max = 1, initial = 0, width = 300)
                        setTkProgressBar(mypb, 1, title = j, label = NULL)
                      }
                      
                      
                      #Get all permutation scenarios
                      if(length(DAMIDvec) <= minlazy) {
                        #Create all combinations of future dams permeabilities
                        permut <- DamX[ESTAGIO_1 == "Planned", do.call(CJ, rep(list(c(0.1, 1)), .N))] #The equivalent of expand.grid but faster
                        colnames(permut) <- DAMIDvec
                        
                        #Add all existing dams
                        permut[, as.character(DamX[ESTAGIO_1 == "Operation", as.character(DAMID)]) := 0.1]
                        
                      } else { #If > minlazy dams, the number of possible permutations leads to memory errors
                        permutlazy <- DamX[ESTAGIO_1 == "Planned", do.call(lazyExpandGrid, rep(list(c(0.1, 1)), .N))] #Get lazy permutation
                        permutsample <- as.data.table(permutlazy$nextItems(sample(permutlazy$n, size=nsamples))) #Sample 10000 permutations
                        colnames(permutsample) <- DAMIDvec
                        #For each dam, select all scenarios where dam is absent and create paired scenario where all other dams are unchanged, but dam is present
                        permut <- ldply(DAMIDvec, function(dam) {
                          rbind(permutsample[get(dam) == 1,],permutsample[get(dam) == 1,]) %>%
                            .[seq((.N/2+1),(.N)), (dam) := 0.1]
                        }) %>% 
                          setDT %>%
                          unique #Remove duplicate scenarios
                        
                        #Add all existing dams
                        permut[, as.character(DamX[ESTAGIO_1 == "Operation", as.character(DAMID)]) := 0.1] #will return a warning if currently no dams 
                      }
                      
                      #Melt and assign all combinations to a given basin-scenario ID to run with data.table
                      scenarios_reg <- melt(cbind(permut, 
                                                  DAMBAS_ID08ext_scenario = paste(j, seq_len(nrow(permut)), sep='_')), 
                                            id.var="DAMBAS_ID08ext_scenario") %>%
                        setnames('variable', 'DAMID') %>%
                        merge(DamX[, .(DAMID, DownSeg, UpSeg)], by='DAMID', allow.cartesian=F)
                      
                      #Run DCI for each scenario
                      DCIX<- scenarios_reg[, 
                                           list(DCI = DCIfunc(
                                             d3=NetX[, list(id=as.character(SEGID),
                                                            l=Shape_Length)],
                                             d2=.SD[, list(id1 = DownSeg,
                                                           id2 = UpSeg,
                                                           pass = value)],
                                             print = F)),
                                           by=.(DAMBAS_ID08ext_scenario)]
                      #For each dam, get DCI diff for all scenarios with and without it (could replace ldply with a data.table solution -TBD)
                      #This just orders rows for each dam's permeability value in the scenario grid and then divides the output array in two (by rows)
                      #The first half of the rows have passability = 0.1 for the ordered dam, the second half has pasability =1
                      #By computing the row-wise difference between the two halves, it gives us the difference for all scenarios at once
                      DCIXpermut <- cbind(DCIX, permut)
                      
                      DCIdiffstats <- ldply(1:length(DAMIDvec), function(i) {
                        colorder <- c(i, (1:length(DAMIDvec))[-i])
                        setorderv(DCIXpermut, cols = DAMIDvec[colorder], na.last=FALSE)
                        DCIdiff <- DCIXpermut[(1+nrow(DCIXpermut)/2):(nrow(DCIXpermut)), DCI] - DCIXpermut[1:(nrow(DCIXpermut)/2), DCI] 
                        return(data.table(
                          DAMIDvec[i],
                          mean(DCIdiff), 
                          max(DCIdiff),
                          min(DCIdiff), 
                          quantile(DCIdiff, 0.975), 
                          quantile(DCIdiff, 0.025),
                          length(DCIdiff)*2))
                      }) %>% #Add dam attributes
                        cbind(DamX[ESTAGIO_1 == "Planned", .(Tipo_1, 
                                                             ESTAGIO_1, 
                                                             POT_KW/1000, 
                                                             NOME, 
                                                             DAMBAS_ID08ext)])
                      
                      colnames(DCIdiffstats) <- c("DAMID", "DCIMeanDiff", "DCIUppLim", "DCIDownLim", "DCIUpCI", "DCIDownCI", "Nscenarios",
                                                  "Type", "Situation", "Capacity","Name", "Basin")

                      #Close progress bar
                      if (length(DAMIDvec) >= minlazy) {close(mypb)}
                      
                      return(DCIdiffstats)
                    }
toc()
big_data <- as.data.table(do.call("rbind", DCIscens))

#Write big_data to fst
write.fst(big_data, file.path(resdir, 'big_data.fst'))