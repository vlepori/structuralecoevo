library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(ggbeeswarm)
library(scatterplot3d)
library(see)
source("R/functions.R")
source("R/figures_funs2.R")

## This is NOT the "modules" pkg that is on CRAN !
# devtools::install_github("klmr/modules")
tb = modules::import("R/toolbox/functions_structural_stability")
tb16 = modules::import("R/toolbox/toolbox_coexistence_2016")
tb22 = modules::import("R/toolbox/toolbox_niche_difference_2022")

## Folder
figfolder <- "tex/figures/"
dir.create(figfolder)

##########
# - Main -
##########

simA = readRDS("out/fig1.rds")
simA = get_metrics(simA) %>% mutate(branching = S-lag(S)) %>%  rename(time = t)

mcA =
  bind_rows(readRDS("out/mc.rds"),
            readRDS("out/mc2.rds"),
            .id = "id")

mcA = get_metrics(mcA) %>%  mutate(time = NA)


#*******************
# ---- FIGURE 1 ----
#*******************

fig1 = figure1(simA)

png(paste0(figfolder,"figure1.png"), width = 10, height = 8, res=300, units = "in")
  fig1
dev.off()

#*******************
# ---- FIGURE 2 ----
#*******************


fig2=
  (figure2(filter(simA, S<6), filter(mcA,S<6))/ figure2(filter(simA, S==6), filter(mcA,S==6)))
  fig2grob = patchwork::patchworkGrob(fig2)

png(paste0(figfolder,"figure3.png"), width = 11, height = 8, res=300, units = "in")
  gridExtra::grid.arrange(fig2grob, 
                        bottom = "Structural niche difference", 
                        left = "Structural fitness difference")
dev.off()


#*******************
# ---- FIGURE 4 ----
#*******************

fig4 = figure4g(mutate(simA,time = time + 2.0),mcA, c(1,2), "sum_n_stars", thin=1, textsize = 4.5)+scale_x_continuous(trans = "sqrt")

ggsave(plot=fig4,paste0(figfolder,"figure4.png"),  width = 11, height = 7)

#*******************
# ---- FIGURE 5 ----
#*******************

fig5 = figure4g_log(filter(simA,S>1),filter(mcA,S>1), c(1,2), "eta", thin=1, labs_y = c(500,220,100),textsize = 4.5)

png(paste0(figfolder,"figure5.png"), width = 11, height = 7, res=300, units = "in")
 fig5
dev.off()


