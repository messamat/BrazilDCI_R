#############         Script 4 to run DCI analysis - September 2019                ###############################
#############        Based on distances calculated in Arcgis - PythonCode         ###############################
#############              Manuscript in preparation - PNAS                       ##############################

####################################################################################################################
###########        Sampling-based PRIORITIZATION ANALYSIS - BRAZIL         ########################################
###########                         (PARETO-FRONT Analysis)                #######################################
###########                        Sample future portifolios               ######################################
###########            (National-level multi-objective optimization)       #####################################
###############################################################################################################
## Choose the number of scenarios to sample
numbScen <- 10000

## Packages
require(tictoc)
require(plyr)
require(bigstatsr)
require(parallel)
require(doParallel)
library(fst)

# # Import network and dams dataset (Mathis folder structure)
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
dcigdb <- file.path(resdir, 'dci.gdb')
outdir_permut = file.path(resdir, 'outpermut_basins')

#Source functions
source(file.path(rootdir, 'src', 'BrazilDCI_R', 'DCIAnalysis.R'))
DCIfunc <- DCIp_opti5 ## Choose what type of DCI function to run (DCIp or DCIi)

#Import formatted data
DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% setDT
NetworkBRAZIL <- read.fst(file.path(resdir, 'NetworkBRAZIL.fst')) %>% setDT

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
sum(log10(basinpermutlen$V1)) #Total number of possible portfolios in brazil 10^682

#Import all scenarios and sort dam names
outdir_permut = file.path(resdir, 'outpermut_basins')
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
                      SHPnum = 0,
                      LHPnum = 0,
                      DamIDs = NA,
                      prevfree = NA)] %>%
  setnames('HYBAS_ID08ext', 'DAMBAS_ID08ext') 

#Sample number of scenarios
scen_numdams<- sample(x = 1:MaxFutDams, size = numbScen, replace=T)

#First sample number of dams, then sample dams
tic()
cl <- parallel::makeCluster(bigstatsr::nb_cores()) #make cluster based on recommended number of cores
on.exit(stopCluster(cl))
doParallel::registerDoParallel(cl)

DCIscens <- foreach(i=scen_numdams, 
                    .packages = c("data.table", "magrittr")) %dopar% {

                      scenbasin <- FutDams[, .SD[sample(.N, i, FALSE, NULL)]] %>% #Sample nsamp dams throughout brazil
                        setorder(DAMID) %>%
                        .[, .(DamIDs = toString(DAMID)) , by=DAMBAS_ID08ext] %>%
                        .[, DAMBAS_ID08ext := NULL] %>%
                        setkey(DamIDs) 
                      
                      return(
                        allscens[scenbasin] %>% #Match with specific scenario DCI based on list of dams 
                          rbind(DCI_L8_current[!(DAMBAS_ID08ext %in% .$DAMBAS_ID08ext),]) %>% #Add other basins without new dams
                          setorder(DAMBAS_ID08ext) %>%
                          .[, list(NatAverageDCI = mean(DCI, na.rm=T), #Get statistics
                                   AddCapacity = sum(POT_KWbasin, na.rm=T),
                                   NFutDams = sum(ndams, na.rm=T), 
                                   NFutSHP = sum(SHPnum, na.rm=T),
                                   NFutLHP = sum(LHPnum, na.rm=T),
                                   NFreeDammed = sum(prevfree, na.rm=T),
                                   scenIDs = toString(scenbasin))]
                      )} %>%
  do.call(rbind, .) %>% setDT
stopCluster(cl)
toc()

write.fst(DCIscens, file.path(resdir, paste0("NationalScen_DCIp_", Sys.Date(), ".fst")))
