#Import packages
source('00.DCI_packages.R')
#Import directory structure and functions
source('00.DCI_functions.R')
DCIfunc <- DCIp_opti5 ## Choose what type of DCI function to run (DCIp or DCIi)

#Import formatted data
DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% setDT
NetworkBRAZIL <- read.fst(file.path(resdir, 'NetworkBRAZIL.fst')) %>% setDT

#Create output directory
outdir_sensitivity = file.path(resdir, 'outsensitivity_passscenarios')
if (!dir.exists(outdir_sensitivity)) {
  print(paste0('Create ', outdir_sensitivity, '...'))
  dir.create(outdir_sensitivity)
}

#Functions
dampassformat <- function(damdt, pasvec = c(0.1, 0.1)) {
  SHPpas <- pasvec[1]
  LHPpas <- pasvec[2]
  damdt[, `:=`(Allcurrent = ifelse(ESTAGIO_1 == 'Operation', 
                                   ifelse(Tipo_1 == 'SHP', SHPpas, LHPpas),
                                   1),
               Allfuture = ifelse(Tipo_1 == 'SHP', SHPpas, LHPpas),
               SHPcurrent = ifelse(ESTAGIO_1 == 'Operation' & Tipo_1 == 'SHP', 
                                   SHPpas, 1),
               LHPcurrent = ifelse(ESTAGIO_1 == 'Operation' & Tipo_1 == 'LHP', 
                                   LHPpas, 1),
               SHPfuture = ifelse(Tipo_1 == 'SHP', SHPpas, 1),
               LHPfuture = ifelse(Tipo_1 == 'LHP', LHPpas, 1),
               passpar = toString(pasvec)
  )] %>%
    melt(id.vars=c('DAMID', 'UpSeg', 'DownSeg', 'SEGID', 'DAMBAS_ID08ext', 'passpar'),
         measure.vars=c('Allcurrent', 'Allfuture', 'SHPcurrent', 'LHPcurrent', 'SHPfuture', 'LHPfuture')) %>%
    setnames(c('variable', 'value'), c('scenario', 'passability'))
}

DCI_bybasscen <- function(netdt, damformatdt, scencol = 'scenario') {
  damformatdt[,
              list(DCI = DCIfunc(d3=netdt[HYBAS_ID08ext==DAMBAS_ID08ext,
                                          list(id=as.character(SEGID), 
                                               l=Shape_Length)],
                                 d2=.SD[, list(id1 = DownSeg, 
                                               id2 = UpSeg, 
                                               pass = passability)],
                                 print = F)),
              by=c('DAMBAS_ID08ext', 'passpar', scencol)]
}


#Generates all combinations of passability for LHP and SHP for SHP<=LHP passability
allpass <- do.call(CJ, rep(list(seq(0, 0.5, 0.02)), 2)) %>%
  setnames(c('V1', 'V2'), c('SHPpas', 'LHPpas')) %>%
  .[SHPpas>=LHPpas]


#Launch DCI analysis for all combinations
tic()
#Initialize parallel analysis
cl <- parallel::makeCluster(bigstatsr::nb_cores()) #make cluster based on recommended number of cores
on.exit(stopCluster(cl))
doParallel::registerDoParallel(cl)

DCIsensitivity <- foreach(j=allpass[,.I], #basinList, ## Loop over basins
                          .packages = c("data.table", "Rcpp", "RcppAlgos", "magrittr", 'plyr','fst', 'igraph')) %dopar% {
                            
                            overwrite=FALSE
                            outf <-  file.path(outdir_sensitivity, paste0('DCIsensitivity_', paste(allpass[j],collapse='_'), '.fst'))
                            
                            if (!file.exists(outf) | overwrite==TRUE) {
                              DCIpasscen <- DCI_bybasscen(netdt = NetworkBRAZIL,
                                                          damformatdt = dampassformat(damdt = DamAttributes, 
                                                                                      pasvec = unlist(allpass[j])),
                                                          scencol = 'scenario')
                              write.fst(DCIpasscen, outf)
                            }
                          }
stopCluster(cl)
toc()