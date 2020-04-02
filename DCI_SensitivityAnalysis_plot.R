library(rprojroot)
library(magrittr)
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
  unique

#Identify basins that are currently free of hydropower (1), currently regulated, and those with both SHPs and LHPs in the future
DCIsensitivity[, `:=`(prevfree = +(DAMBAS_ID08ext %in% 
                                     DCIsensitivity[scenario=='Allcurrent' & DCI==100, 
                                                    DAMBAS_ID08ext]),
                      bothtypes = +(DAMBAS_ID08ext %in% 
                                      DCIsensitivity[(scenario=='SHPfuture' & DCI!=100) | (scenario=='LHPfuture' & DCI!=100), 
                                                     DAMBAS_ID08ext]))]

#calculate differences between future and current
DCIsensitivity[, `:=`(scentime=ifelse(grepl('current', scenario), 'current', 'future'), #Separate scenarios into type and time
                      scentype = gsub('current|future', '', scenario))]

DCIloss <- dcast(DCIsensitivity, DAMBAS_ID08ext+passpar+prevfree+bothtypes+scentype~scentime, value.var = 'DCI')  %>% #Compute loss
  .[, `:=`(DCIloss = future-current,
           DCIlossperc = 100*(future-current)/current)]
      
#Compute average brazil-wide DCI percentage loss for all basins, previously free of dams and not previously free of dams
passcols <- c("SHP passability","LHP passability")
DCIlossstats <- DCIloss[, list(meanloss = mean(DCIlossperc)), by=.(passpar, prevfree, scentype)] %>% #Compute DCIloss for previously free vs not previously free
  rbind(DCIloss[, list(meanloss = mean(DCIloss), prevfree = 2), by=.(passpar, scentype)]) %>%  #bind with statistics for all together
  .[, prevfree := factor(prevfree, labels = c("Regulated basins", "Hydro-free basins", "All basins"))] %>%
  .[, (passcols) := tstrsplit(passpar, ", ", fixed=TRUE)] %>% #Separate passability scenarios into two columns
  .[, (passcols) := lapply(.SD, as.numeric), .SDcols = passcols] %>% #Convert to numeric
  setorderv(passcols, order=-1L) #Order in descending order the passability so that low passability shows on top


#Plot
tiff(filename = file.path(resdir, "Figure_sensitivity.tiff"), height = 2000, width = 3200, res = 300, compression = c("lzw"))
ggplot(DCIlossstats, aes(x=`SHP passability`, y=-meanloss, color=`LHP passability`)) + 
  geom_point(size=2, alpha=1/2) + 
  facet_grid(prevfree~scentype, labeller=label_value) + 
  scale_y_sqrt(name='Loss in river connectivity (%)', expand=c(0,0), breaks=c(0,1,5,10,20,30,40)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_distiller(palette='Spectral') +
  coord_cartesian(clip = 'off') +
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(.65, "cm"),
        text = element_text(size=15),
        strip.background = element_rect(fill='white', color='white'),
        strip.text = element_text(face="bold"))
dev.off()