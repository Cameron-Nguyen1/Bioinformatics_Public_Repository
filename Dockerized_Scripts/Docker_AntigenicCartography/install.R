dep = c("devtools","ggplot2","stringr","reshape2","RColorBrewer","paletteer","optparse")
for (i in dep){install.packages(i,Ncpus=4,INSTALL_opts = '--no-lock')}
library(devtools)
devtools::install_github("acorg/Racmacs")