#To do: create function to prepare dam passability for a scenario based on year, dam size, and passability

## Packages
require(tictoc)
require(plyr)

# Import network and dams dataset (Mathis folder structure)
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
dcigdb <- file.path(resdir, 'dci.gdb')

#Source functions
source(file.path(rootdir, 'src', 'BrazilDCI_R', 'DCIAnalysis.R'))
DCIfunc <- DCIp_opti3 ## Choose what type of DCI function to run (DCIp or DCIi)

#Import formatted data
DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% setDT
NetworkBRAZIL <- read.fst(file.path(resdir, 'NetworkBRAZIL.fst')) %>% setDT

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

DCI_bybasscen <- function(netdt, damdt, scencol = 'scenario') {
  damdt[,
        list(DCI = DCIfunc(d3=netdt[HYBAS_ID08ext==DAMBAS_ID08ext,
                                    list(id=as.character(SEGID), 
                                         l=Shape_Length)],
                           d2=.SD[, list(id1 = DownSeg, 
                                         id2 = UpSeg, 
                                         pass = passability)],
                           print = F),
             passpar),
        by=.(DAMBAS_ID08ext, get(scencol))]
}


#Generates all combinations of passability for LHP and SHP for SHP<=LHP passability
allpass <- do.call(CJ, rep(list(seq(0, 0.5, 0.1)), 2)) %>%
  setnames(c('V1', 'V2'), c('SHPpas', 'LHPpas')) %>%
  .[SHPpas>=LHPpas]

#Launch DCI analysis for all combinations
DCI_bybasscen(netdt = NetworkBRAZIL,
              damdt = dampassformat(DamAttributes, pasvec = c(0.1, 0.1)),
              scencol = 'scenario')