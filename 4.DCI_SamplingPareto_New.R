#############         Script 4 to run DCI analysis - September 2019                ###############################
#############        Based on distances calculated in Arcgis - PythonCode         ###############################
#############              Manuscript in preparation - PNAS                       ##############################

####################################################################################################################
###########        Sampling-based PRIORITIZATION ANALYSIS - BRAZIL         ########################################
###########                         (PARETO-FRONT Analysis)                #######################################
###########                        Sample future portifolios               ######################################
###########            (National-level multi-objective optimization)       #####################################
###############################################################################################################
## Choose the total number of scenarios to sample (includes already sampled scenarios)
numbScen <- 1000000

#Import packages
source('00.DCI_packages.R')
#Import directory structure and functions
source('00.DCI_functions.R')

#Import formatted data
DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% setDT
NetworkBRAZIL <- read.fst(file.path(resdir, 'NetworkBRAZIL.fst')) %>% setDT


sample_pareto <- function(DamAttributes, NetworkBRAZIL, DCIfunc, DCIname,
                          checkprevious=F) {
  #Run DCI for current scenario
  #Compute DCI for Allcurrent scenario
  DCI_L8_current <- NetworkBRAZIL[,
                                  list(DCI = DCIfunc(
                                    d2= DamAttributes[DAMBAS_ID08ext == HYBAS_ID08ext,
                                                      list(
                                                        id1 = DownSeg,
                                                        id2 = UpSeg,
                                                        pass = Allcurrent)],
                                    d3 = .SD[, list(id=as.character(SEGID),
                                                    l=Shape_Length)],
                                    print = F)),
                                  by=.(HYBAS_ID08ext)]
  
  ## Get the total number of future dams
  FutDams <- DamAttributes[ESTAGIO_1 == "Planned", .(DAMBAS_ID08ext, DAMID)] %>% setorder(DAMID)
  MaxFutDams <- FutDams[,.N]
  
  ## Get the total number of possible portfolios for each basin
  basinpermutlen <- DamAttributes[ESTAGIO_1 == "Planned", 2^(.N), by=DAMBAS_ID08ext]
  print('Total number of possible portfolios in Brazil:')
  sum(log10(basinpermutlen$V1)) #Total number of possible portfolios in brazil 10^682
  
  #Import all scenarios and sort dam names
  outdir_permut = file.path(resdir, paste0('outpermut_basins', DCIname))
  allscens <- lapply(file.path(outdir_permut, list.files(outdir_permut)), read.fst) %>% #Read all scenarios for each basin
    do.call(rbind, .) %>% #Join them together
    setDT %>%
    .[!is.na(DamIDs), DamIDs := toString(sort(str_split(string=DamIDs, pattern='[,]\\s*',simplify=T))), #Order dam names
      by=.(scenbasin, DAMBAS_ID08ext)] %>%
    setkey(DamIDs) #Set key
  
  #Prepare current DCI DT for merging with future sampled scenarios
  DCI_L8_current[, `:=`(scenbasin = 0,
                        ndams=0,
                        POT_KWbasin = 0,
                        Guarantee = 0,
                        SHPnum = 0,
                        LHPnum = 0,
                        DamIDs = NA,
                        prevfree = NA)] %>%
    setnames('HYBAS_ID08ext', 'DAMBAS_ID08ext') 
  
  
  #---------------- Check whether scenarios already exist and adjust number to sample ------
  if (checkprevious) {
    nationalfiles <- list.files(path=resdir,
                                pattern= paste0("NationalScen_",
                                                DCIname, "_.*[.]fst"))
    if (length(nationalfiles) > 0) {
      for (nfile in nationalfiles) {
        NationalScenarios <- read.fst(file.path(resdir, nfile)) %>%
          setDT %>%
          unique #Only keep unique scenarios (mostly removes the multiple scenarios with all dams or no dams)
        numbScen = numbScen - NationalScenarios[substr(scenIDs, 1,2) != 'NA', .N]
        remove(NationalScenarios)
        gc()
      }
    }
  }

  
  #----------------- Launch sampling -------------------------------------------------------
  #Sample number of scenarios
  scen_numdams<- sample(x = 1:MaxFutDams, size = numbScen, replace=T)
  scens_sampledbas <- allscens[DAMBAS_ID08ext=='60808111301',] #All scenarios that have DCI for basin ID '60808111301'
  
  #First sample number of dams, then sample dams
  tic()
  cl <- parallel::makeCluster(1) #make cluster based on recommended number of cores
  on.exit(stopCluster(cl))
  doParallel::registerDoParallel(cl)
  
  DCIscens <- foreach(i=scen_numdams, 
                      .packages = c("data.table", "magrittr", "stringr")) %dopar% {
                        
                        scendams <- FutDams[, .SD[sample(.N, i, FALSE, NULL)]] %>% #Sample nsamp dams throughout brazil
                          setorder(DAMID) %>%
                          .[, .(DamIDs = toString(DAMID)) , by=DAMBAS_ID08ext] %>%
                          setkey(DamIDs) 
                        
                        ndams_sampledbas <- scendams[DAMBAS_ID08ext == '60808111301', 
                                                     length(str_split(string=DamIDs, pattern='[,]\\s*', simplify=T))]
                        
                        selectscen_sampledbas <- scens_sampledbas[ndams==ndams_sampledbas,][sample(.N, 1), DamIDs]
                        
                        if (length(selectscen_sampledbas) != 0 | ndams_sampledbas == 0) {
                          scendams[DAMBAS_ID08ext == '60808111301',
                                   DamIDs := selectscen_sampledbas] %>%
                            .[, DAMBAS_ID08ext := NULL] %>%
                            setkey(DamIDs) 
                          
                          #head(allscens)
                          return(
                            allscens[scendams] %>% #Match with specific scenario DCI based on list of dams 
                              rbind(DCI_L8_current[!(DAMBAS_ID08ext %in% .$DAMBAS_ID08ext),]) %>% #Add other basins without new dams
                              setorder(DAMBAS_ID08ext) %>%
                              .[, list(NatAverageDCI = mean(DCI, na.rm=T), #Get statistics
                                       AddCapacity = sum(POT_KWbasin, na.rm=T),
                                       AddGuarantee = sum(Guarantee, na.rm=T),
                                       NFutDams = sum(ndams, na.rm=T), 
                                       NFutSHP = sum(SHPnum, na.rm=T),
                                       NFutLHP = sum(LHPnum, na.rm=T),
                                       NFreeDammed = sum(prevfree, na.rm=T),
                                       scenIDs = paste(scenbasin, collapse=','))]
                          )
                        } 
                        else {
                          return(NULL)
                        }
                      } %>%
    do.call(rbind, .) %>% setDT
  stopCluster(cl)
  toc()
  
  write.fst(DCIscens, file.path(resdir, paste0("NationalScen_", 
                                               DCIname, "_", 
                                               format(Sys.time(), "%Y%m%d_%Hh"),
                                               ".fst")))
}

for (i in 1:5) {
  sample_pareto(DamAttributes = DamAttributes,
                NetworkBRAZIL = NetworkBRAZIL,
                DCIfunc = DCIp_opti5,
                DCIname = 'DCIp',
                checkprevious=F)
  gc()
}


sample_pareto(DamAttributes = DamAttributes,
              NetworkBRAZIL = NetworkBRAZIL,
              DCIfunc = DCIi_opti,
              DCIname = 'DCIi')