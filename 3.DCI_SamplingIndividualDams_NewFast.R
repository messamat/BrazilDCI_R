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

#Import packages
source('00.DCI_packages.R')
#Import directory structure and functions
source('00.DCI_functions.R')

#Import formatted data
DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% setDT
NetworkBRAZIL <- read.fst(file.path(resdir, 'NetworkBRAZIL.fst')) %>% setDT

sample_indivdams <- function(DamAttributes, NetworkBRAZIL, DCIfunc, DCIname, 
                             minlazy, nsamples) {
  
  #Create output directory
  outdir_permut = file.path(resdir, paste0('outpermut_basins', DCIname))
  if (!dir.exists(outdir_permut)) {
    print(paste0('Create ', outdir_permut, '...'))
    dir.create(outdir_permut)
  }
  
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
  basinList <- DamAttributes[ESTAGIO_1 == "Planned", .N, by=DAMBAS_ID08ext] %>% 
    setorder(N) %>% 
    .[,DAMBAS_ID08ext]
  
  #Remove the basins that will have more than X future dams if needed
  # X <- 10
  # basinList <- DamAttributes[ESTAGIO_1 == "Planned", .N, by=DAMBAS_ID08ext][N<X, DAMBAS_ID08ext]
  
  #### Parameters to run permutations-based DCI ####
  #minimum number of dams in basin at which the function starts sampling permutation scenarios (for speed/memory sake)
  #minlazy
  #number of seed scenarios to sample from lazy permutation. The actual number of scenario pairs that will actually 
  #be evluated will often be at least >5-20 times more than that number
  #nsamples <- 10000
  
  ############## LAUNCH ANALYSIS ####################
  tic("total")
  #Initialize parallel analysis
  cl <- parallel::makeCluster(bigstatsr::nb_cores()) #make cluster based on recommended number of cores
  on.exit(stopCluster(cl))
  doParallel::registerDoParallel(cl)
  
  #60808111301 has 22 planned dams
  DCIscens <-
    foreach(j=basinList, #basinList, ## Loop over basins
            .packages = c("data.table", "Rcpp", "RcppAlgos",
                          "magrittr", 'plyr','tcltk', 'fst'),
            .export = 'lazyExpandGrid') %dopar% {
              
              ## filter attributes of the basin j         
              NetX <- NetworkBRAZIL[HYBAS_ID08ext == j, ]
              DamX <- DamAttributes[DAMBAS_ID08ext == j, ]
              DAMIDvec <- as.character(DamX[ESTAGIO_1 == "Planned", DAMID])
              
              #Was this basin free flowing before?
              prevfreebas <- DamX[ESTAGIO_1=='Operation', .N]==0
              
              #-------------------------------------------------------------------
              #Get all permutation scenarios
              if(length(DAMIDvec) <= minlazy) {
                #Create all combinations of future dams permeabilities
                permut <- DamX[ESTAGIO_1 == "Planned",
                               do.call(CJ, rep(list(c(0.1, 1)), .N))] #The equivalent of expand.grid but faster
                colnames(permut) <- DAMIDvec
                
                #Add all existing dams
                permut[, as.character(DamX[ESTAGIO_1 == "Operation", 
                                           as.character(DAMID)]) := 0.1]
                
              } else { #If > minlazy dams, the number of possible permutations leads to memory errors
                #Get lazy permutation — like expand.grid but for cases where  too many possibilities to possibly thold in memory
                permutlazy <- DamX[ESTAGIO_1 == "Planned",
                                   do.call(lazyExpandGrid, rep(list(c(0.1, 1)), .N))] 
                permutsample <- as.data.table(permutlazy$nextItems(
                  sample(permutlazy$n, size=nsamples))) #Sample 10000 permutations
                colnames(permutsample) <- DAMIDvec
                #For each dam, select scenarios where dam is absent & create paired scen. where other dams are unchanged, but dam is present
                permut <- ldply(DAMIDvec, function(dam) {
                  rbind(permutsample[get(dam) == 1,],
                        permutsample[get(dam) == 1,]) %>%
                    .[seq((.N/2+1),(.N)), (dam) := 0.1]
                }) %>% 
                  setDT %>%
                  unique #Remove duplicate scenarios
                
                #Add all existing dams
                permut[, as.character(DamX[ESTAGIO_1 == "Operation", 
                                           as.character(DAMID)]) := 0.1] #will return warning if currently no dams 
              }
              
              #-------------------------------------------------------------------
              #Melt and assign all combinations to a given basin-scenario ID to run with data.table
              scenarios_reg <- melt(cbind(permut, 
                                          scenbasin = seq_len(nrow(permut))), 
                                    id.var="scenbasin") %>%
                setnames('variable', 'DAMID') %>%
                merge(DamX[, .(DAMID, DownSeg, UpSeg, Tipo_1, POT_KW, ESTAGIO_1)], 
                      by='DAMID', allow.cartesian=F)
              
              #Run DCI for each scenario
              DCIX<- scenarios_reg[
                , list(DCI = DCIfunc(d3=NetX[, list(id=as.character(SEGID), l=Shape_Length)],
                                     d2=.SD[, list(id1 = DownSeg, id2 = UpSeg, pass = value)],
                                     print = F)),
                by=.(scenbasin)]
              
              #-------------------------------------------------------------------
              #Get statistics on scenario for portfolio analysis
              scenarios_reg[ESTAGIO_1=='Planned' & value==0.1,
                            list(DAMBAS_ID08ext=j,
                                 ndams=.N,
                                 POT_KWbasin = sum(POT_KW/1000),
                                 Guarantee = sum(fifelse(Tipo_1 == 'LHP', 0.55, 0.6)*POT_KW/1000),
                                 SHPnum = sum(Tipo_1=='SHP'),
                                 LHPnum = sum(Tipo_1=='LHP'),
                                 DamIDs = toString(DAMID),
                                 prevfree = prevfreebas),
                            by=.(scenbasin)] %>%
                merge(DCIX[, DAMBAS_ID08ext := j], on='scenbasin', all.y=T) %>%
                .[is.na(ndams), `:=`(ndams=0, POT_KWbasin=0, SHPnum=0, LHPnum=0)] %>%
                write.fst(file.path(outdir_permut, paste0('permut_', j, '.fst')))
              
              #-------------------------------------------------------------------
              #For each dam, get DCI diff for all scenarios with and without it (could replace ldply with a data.table solution -TBD)
              #This just orders rows for each dam's permeability value in the scenario grid and then divides the output array in two (by rows)
              #The first half of the rows have passability = 0.1 for the ordered dam, the second half has pasability =1
              #By computing the row-wise difference between the two halves, it gives us the difference for all scenarios at once
              DCIXpermut <- cbind(DCIX, permut)
              
              DCIdiffstats <- ldply(1:length(DAMIDvec), function(i) {
                colorder <- c(i, (1:length(DAMIDvec))[-i])
                setorderv(DCIXpermut, cols = DAMIDvec[colorder], na.last=FALSE)
                DCIdiff <- (DCIXpermut[(1+nrow(DCIXpermut)/2):(nrow(DCIXpermut)), DCI] - 
                              DCIXpermut[1:(nrow(DCIXpermut)/2), DCI])
                
                return(data.table(
                  DAMIDvec[i],
                  mean(DCIdiff), 
                  max(DCIdiff),
                  min(DCIdiff), 
                  quantile(DCIdiff, 0.975), 
                  quantile(DCIdiff, 0.025),
                  length(DCIdiff)*2))
              }) %>% #Add dam attributes
                cbind(DamX[ESTAGIO_1 == "Planned", .(NOME, 
                                                     DAMBAS_ID08ext,
                                                     Tipo_1, 
                                                     ESTAGIO_1, 
                                                     POT_KW/1000)]) %>%
                setDT
              
              colnames(DCIdiffstats) <- c("DAMID", "DCIMeanDiff", "DCIUppLim", 
                                          "DCIDownLim", "DCIUpCI", "DCIDownCI", 
                                          "Nscenarios", "Name", "Basin",
                                          "Type", "Situation", "Capacity")
              
              DCIdiffstats[, Guarantee :=
                             fifelse(Type == 'LHP', 0.55, 0.6)*Capacity]
              
              return(DCIdiffstats)
            }
  stopCluster(cl)
  toc()
  big_data <- as.data.table(do.call("rbind", DCIscens))
  
  #Write big_data to fst
  out_fst = file.path(resdir, 
                      paste0('SamplingIndividualDams_results_', DCIname, 
                             Sys.Date(), '.fst'))
  write.fst(big_data, out_fst)
}

sample_indivdams(DamAttributes = DamAttributes,
                 NetworkBRAZIL = NetworkBRAZIL,
                 DCIfunc = DCIp_opti5,
                 DCIname = 'DCIp',
                 minlazy=21,
                 nsamples=10000)

sample_indivdams(DamAttributes = DamAttributes,
                 NetworkBRAZIL = NetworkBRAZIL,
                 DCIfunc = DCIi_opti,
                 DCIname = 'DCIi',
                 minlazy=21,
                 nsamples=10000)

