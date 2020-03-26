

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


DamAttributes[, `:=`(All_current = ifelse(DamAttributes$ESTAGIO_1 == 'Operation', 0.1, 1),
                     All_future = 1,
                     SHP_current = ifelse(DamAttributes$ESTAGIO_1 == 'Operation' & 
                                            DamAttributes$Tipo_1 == 'SHP', 0.1, 1),
                     LHP_current = ifelse(DamAttributes$ESTAGIO_1 == 'Operation' & 
                                            DamAttributes$Tipo_1 == 'LHP', 0.1, 1),
                     SHP_future = ifelse(DamAttributes$Tipo_1 == 'SHP', 0.1, 1),
                     LHP_future = ifelse(DamAttributes$Tipo_1 == 'LHP', 0.1, 1)
)]