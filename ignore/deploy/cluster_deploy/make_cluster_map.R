#########################################################################
# Purpose: Make Cluster Map
#
# Author: Nicholas F. Brazeau
#
# Date: April 17 2020
#########################################################################
library(tidyverse)
lvls <- c("Cumulative", "Time-Series")
currday <- c(25, 50, 100)
prior <- c("intermed", "stonger", "strong")
burnin <- c(1e3, 1e4)
rungs <- c(1, 10, 25)
gti <- c(2, 3)

mstrmap <- expand.grid(list(lvls, currday, prior, burnin, rungs, gti), stringsAsFactors = F)
colnames(mstrmap) <- c("lvl", "curr_day", "prior", "burnin", "rungs", "gti")
mstrmap <- mstrmap %>%
  dplyr::mutate(model = paste0("paramset", 1:nrow(mstrmap)))

reps <- 2
mstrmap <- lapply(1:reps, function(x){mstrmap}) %>%
  dplyr::bind_rows(.) %>%
  dplyr::group_by(model) %>%
  dplyr::mutate(iter = 1:n())


#..............................................................
# out
#..............................................................
outnames <- paste0(mstrmap$model, "_iter_", mstrmap$iter)
mstrmap.list <- split(mstrmap, 1:nrow(mstrmap))
mapply(function(x, y){
  saveRDS(object = x, file = paste0("ignore/deploy/cluster_deploy/paramset/", y, ".rds"))
}, x = mstrmap.list, y = outnames)

write.table(x = mstrmap,
            file = "ignore/deploy/cluster_deploy/paramset/mastermap_map.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)


snakemap <- as.data.frame(outnames)
colnames(snakemap) <- c("#paramset")
write.table(x = snakemap,
            file = "ignore/deploy/cluster_deploy/paramset/snake_map.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)
