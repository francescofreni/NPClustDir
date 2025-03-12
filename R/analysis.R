# LOADING PACKAGES #############################################################
library(pdfCluster)
library(mvtnorm)
library(rgl)
library(icosa)
library(plotly)
library(igraph)
library(ks)
library(tidyverse)
library(clue)
source('utils.R')
source('functions.R')


# MESH DIVISION AND DATA #######################################################
ENERGY <- c(alldata_new$ENERGY, bckF_sim$ENERGY)
ID <- c(alldata_new$ID, bckF_sim$ID)
PSF_T <- c(alldata_new$PSF_T, bckF_sim$PSF_T)
sorg <- gal2cart(alldata_new[,c('L', 'B')])
bkg <- gal2cart(bckF_sim[,c('L', 'B')])
all <- rbind(sorg, bkg)
posizioni9 <- locate(grid_9_20, all)
posizioni11 <- locate(grid_11_20, all)
posizioni9.bkg <- locate(grid_9_20, bkg)
posizioni11.bkg <- locate(grid_11_20, bkg)

# Final, tessellation (9,20), elimination of the two connected components
# with empty cells
faces3d(final, col='blue')

# Largest connected component
faces3d(central, col='blue')
pos_central <- which(locate(central, all) != 'F0')
energy_central <- ENERGY[pos_central]
id_central <- ID[pos_central]
psf_t_central <- PSF_T[pos_central]
dati_central <- all[pos_central,]
points3d(dati_central[,1], dati_central[,2], dati_central[,3], col='red', size=3)

# Poles found by excluding the central portion (central) from final
poli <- final[!(row.names(final@faces) %in% row.names(central@faces))]
pos_poli <- which(locate(poli, all) != 'F0')
id_poli <- ID[pos_poli]
dati_poli <- all[pos_poli,]
sorg_poli <- sorg[which(locate(poli, sorg) != 'F0'),]
faces3d(poli, col='blue')
points3d(dati_poli[,1], dati_poli[,2], dati_poli[,3], col='red', size=3)

# Connected components found starting from the largest connected component
# removing empty cells and those with one or two observations
# Tessellation (9,20)
faces3d(final_central, col='blue')
pos_final_central <- which(locate(final_central, all) != 'F0')
energy_final_central <- ENERGY[pos_final_central]
id_final_central <- ID[pos_final_central]
psf_t_final_central <- PSF_T[pos_final_central]
dati_final_central <- all[pos_final_central,]
points3d(dati_final_central[,1], dati_final_central[,2], dati_final_central[,3],
         col='red', size=3)

# Largest component among those remaining, tessellation (9,20)
faces3d(final_component, col='blue')
pos_final_component <- which(locate(final_component, all) != 'F0')
energy_final_component <- ENERGY[pos_final_component]
id_final_component <- ID[pos_final_component]
psf_t_final_component <- PSF_T[pos_final_component]
dati_final_component <- all[pos_final_component,]
points3d(dati_final_component[,1], dati_final_component[,2], dati_final_component[,3],
         col='red', size=3)

# Largest central component, tessellation (11,20). This replaces
# the central, obtained by removing connected components with empty cells
faces3d(central_2, col='blue')
pos_central_2 <- which(locate(central_2, all) != 'F0')
energy_central_2 <- ENERGY[pos_central_2]
id_central_2 <- ID[pos_central_2]
psf_t_central_2 <- PSF_T[pos_central_2]
dati_central_2 <- all[pos_central_2,]
points3d(dati_central_2[,1], dati_central_2[,2], dati_central_2[,3], col='red', size=3)

# Connected components found starting from final_11_20, removing empty cells
# and those with one or two observations
faces3d(final_central_2, col='blue')
pos_final_central_2 <- which(locate(final_central_2, all) != 'F0')
energy_final_central_2 <- ENERGY[pos_final_central_2]
id_final_central_2 <- ID[pos_final_central_2]
psf_t_final_central_2 <- PSF_T[pos_final_central_2]
dati_final_central_2 <- all[pos_final_central_2,]
points3d(dati_final_central_2[,1], dati_final_central_2[,2], dati_final_central_2[,3],
         col='red', size=3)

# Largest connected component in final_central_2
faces3d(final_component_2, col='blue')
pos_final_component_2 <- which(locate(final_component_2, all) != 'F0')
energy_final_component_2 <- ENERGY[pos_final_component_2]
id_final_component_2 <- ID[pos_final_component_2]
psf_t_final_component_2 <- PSF_T[pos_final_component_2]
dati_final_component_2 <- all[pos_final_component_2,]
sorg_final_component_2 <- sorg[which(locate(final_component_2, sorg) != 'F0'),]
points3d(dati_final_component_2[,1], dati_final_component_2[,2], dati_final_component_2[,3],
         col='red', size=3)

# Central region, largest connected component excluded
final_central_3 <- final_central_2[!(row.names(final_central_2@faces) %in% 
                                       row.names(final_component_2@faces))]
faces3d(final_central_3, col='blue')
pos_final_central_3 <- which(locate(final_central_3, all) != 'F0')
id_final_central_3 <- ID[pos_final_central_3]
dati_final_central_3 <- all[pos_final_central_3,]
sorg_final_central_3 <- sorg[which(locate(final_central_3, sorg) != 'F0'),]
sorg_final_central_2 <- sorg[which(locate(final_central_2, sorg) != 'F0'),]

# The regions I consider are:
# - poli
# - final_central_3
# - final_component_2
faces3d(poli)
faces3d(final_central_3, col='blue')
faces3d(final_component_2, col='red')
faces3d(grid_11_20, col='black', lit=F)
points3d(sorg[,1], sorg[,2], sorg[,3], col='orange', size=2)
points3d(bkg[,1], bkg[,2], bkg[,3], col='lightblue', size=2)


# ANALYSIS #####################################################################
# Vector that will contain the cluster ID to which each observation belongs. 
# 0 represents the background.
clus_obs <- rep(0, dim(all)[1])

# First, classify the observations in the poles.
com_poli <- components(poli@graph)
ret_poli <- vector(length=count_components(poli@graph), mode='list')
limit <- unique(com_poli$membership)
for (cc in limit) {
  print(cc)
  nomi <- names(com_poli$membership[com_poli$membership==cc])
  if (length(nomi) <= 2) {
    pos <- posizioni9 %in% nomi
    n <- sum(pos)
    if (n >= 5) {
      clus_obs[pos] <- max(clus_obs) + 1
    }
  }
  if (length(nomi) > 2) {
    pos <- which(posizioni9 %in% nomi)
    dati <- all[pos,]
    ret <- pdfCluster.fermi.2(x=dati, mesh=poli, h=hns(dati), nomi=nomi)
    ret_poli[[cc]] <- ret
    obs_ret <- ret@obs.cluster
    obs_ret[is.na(obs_ret)] = 0
    obs_ret[obs_ret != 0] <- obs_ret[obs_ret != 0] + max(clus_obs)
    clus_obs[pos] <- obs_ret
  }
}

save(clus_obs, ret_poli, file='risultati_poli2.RData')
write.table(table(clus_obs[pos_poli], id_poli), file='tabella_poli2.txt')

ret_poli[[25820]]
plot(ret_poli[[25820]]@tree, center=T)
table(clus_obs[pos_prova_603], id_prova_603)

x_poli <- clus_obs[pos_poli]
x_poli[x_poli == '0'] <- 'bck'
sum(x_poli == id_poli) # This is the correctly classified background
# 48771/73590 == 0.663
sum((x_poli != 'bck') == (id_poli != 'bck')) # Correctly classified observations
sum((x_poli != 'bck') & (id_poli != 'bck')) # Correctly classified sources
sum(id_poli!='bck') # These are the true sources
sum(id_poli=='bck')

# Classify the observations in the central region, excluding the largest component.
com_central <- components(final_central_3@graph)
ret_central <- vector(length=count_components(final_central_3@graph), mode='list')
limit2 <- unique(com_central$membership)
for (cc in limit2) {
  print(cc)
  nomi <- names(com_central$membership[com_central$membership==cc])
  if (length(nomi) <= 2) {
    pos <- posizioni11 %in% nomi
    n <- sum(pos)
    if (n >= 5) clus_obs[pos] <- max(clus_obs) + 1
  }
  if (length(nomi) > 2) {
    pos <- which(posizioni11 %in% nomi)
    dati <- all[pos,]
    ret <- pdfCluster.fermi.2(x=dati, mesh=final_central_3, h=hns(dati), nomi=nomi)
    ret_central[[cc]] <- ret
    obs_ret <- ret@obs.cluster
    obs_ret[is.na(obs_ret)] = 0
    obs_ret[obs_ret != 0] <- obs_ret[obs_ret != 0] + max(clus_obs)
    clus_obs[pos] <- obs_ret
  }
}

save(clus_obs, ret_central, file='risultati_central.RData')
write.table(table(clus_obs[pos_final_central_3], id_final_central_3),
            file='tabella_central.txt')

ret_central[[7165]]
which(com_central$csize==2020)

x <- clus_obs[pos_final_central_3]
x[x == '0'] <- 'bck'
sum(x == id_final_central_3) # This is the correctly classified background
sum((x != 'bck') == (id_final_central_3 != 'bck')) # Correctly classified observations
sum((x != 'bck') & (id_final_central_3 != 'bck')) # Correctly classified sources
sum(id_final_central_3!='bck') # These are the true sources
sum(id_final_central_3=='bck')

# Classify the observations in the largest component.
nomi_comp <- row.names(final_component_2@faces)
ret_comp <- pdfCluster.fermi.2(x=dati_final_component_2,
                                    mesh=final_component_2, 
                                    h=hns(dati_final_component_2), 
                                    nomi=nomi_comp)
obs_ret_comp <- ret_comp1@obs.cluster
obs_ret_comp[is.na(obs_ret_comp)] = 0
clus_obs[pos_final_component_2] <- obs_ret_comp

x_comp <- obs_ret_comp
x_comp[x_comp == '0'] <- 'bck'
sum(x_comp == id_final_component_2)
sum((x_comp != 'bck') == (id_final_component_2 != 'bck')) # Correctly classified observations
sum((x_comp != 'bck') & (id_final_component_2 != 'bck')) # Correctly classified sources
sum(id_final_component_2!='bck') # These are the true sources
sum(id_final_component_2=='bck')

write.table(table(obs_ret_comp, id_final_component_2),
            file='tabella_comp_grande.txt')

# Combining results
obs_ret_comp[obs_ret_comp != 0] <- obs_ret_comp[obs_ret_comp != 0] + max(clus_obs)
clus_obs[pos_final_component_2] <- obs_ret_comp
clus_obs[clus_obs==0] <- 'bck' # Assign 'bck' to observations classified as background

sum((clus_obs != 0) & (ID != 'bck')) # Sources correctly classified
sum(ID != 'bck')
66079/73318 # True positive rate

gr_stim <- ifelse(clus_obs == 'bck', 0, 1) # Assign 0 to background, 1 to sources
table(gr_stim)
gr <- ifelse(ID == 'bck', 0, 1) # True classification: 0 for background, 1 for sources
table(gr)
adj.rand.index(gr_stim, gr) # Adjusted Rand Index to measure clustering similarity

sum((clus_obs != 'bck') & (ID == 'bck')) # Background misclassified as a source
sum((clus_obs != 'bck') & (ID == 'bck')) + sum(clus_obs == ID) # Total correctly classified + misclassified
sum((clus_obs != 'bck') & (ID == 'bck')) / (
  sum((clus_obs != 'bck') & (ID == 'bck')) + sum(clus_obs == ID)
) # Proportion of background misclassified as a source

length(ID[ID=='bck']) / length(ID)

table(gr_stim, gr) # Contingency table comparing predicted vs true classifications


# ATTEMPT TO REMOVE THE MISCLASSIFIED BACKGROUND ###############################
# APPROACH 1 ###################################################################
faces3d(poli)
faces3d(final_central_3)
faces3d(final_component_2)

table(com_poli$csize)

una <- poli[row.names(poli@faces) %in%
              names(com_poli$membership[
                com_poli$membership==which(com_poli$csize==1)[1]])]
una
osservazione <- all[which(locate(una, all) != 'F0'),]
osservazione <- t(as.matrix(osservazione))

dens(osservazione, poli, hns(osservazione), "F100050")

pdf_est <- rep(0, dim(all)[1])

# Compute density at the poles
limit <- unique(com_poli$membership)
for (cc in limit) {
  print(cc)
  nomi <- names(com_poli$membership[com_poli$membership==cc])
  pos <- which(posizioni9 %in% nomi)
  dati <- all[pos,]
  if (is.null(dim(dati))) dati <- t(as.matrix(dati)) # un solo dato
  densita <- as.vector(dens(dati, poli, hns(dati), nomi))
  pdf_est[pos] <- densita
}

# Compute density in the central region, excluding the largest connected component
limit2 <- unique(com_central$membership)
for (cc in limit2) {
  print(cc)
  nomi <- names(com_central$membership[com_central$membership==cc])
  pos <- which(posizioni11 %in% nomi)
  dati <- all[pos,]
  if (is.null(dim(dati))) dati <- t(as.matrix(dati)) # un solo dato
  densita <- as.vector(dens(dati, final_central_3, hns(dati), nomi))
  pdf_est[pos] <- densita
}

# Take the results of the largest connected component
ret_comp1@pdf
obs_comp <- locate(final_component_2, dati_final_component_2)
dens_comp <- sapply(1:dim(dati_final_component_2)[1],
                    function(i) ret_comp1@pdf[which(nomi_comp == obs_comp[i])])
pdf_est[pos_final_component_2] <- dens_comp

# Compute density in the regions that were excluded to find the large component
faces3d(central, col='blue')
points3d(dati_central[,1], dati_central[,2], dati_central[,3], col='red', size=2)

plot3d(grid_11_20, guides=F, arcs=T, sphere=1, col='blue')

# Identify empty regions that do not contain observations
empty <- grid_11_20[!(occupied(grid_11_20, dati_central))@values]
com_empty <- components(empty@graph)
table(com_empty$csize)
faces3d(empty)

# Identify and visualize the largest empty components
to_del_1 <- which.max(com_empty$csize)
to_del_2 <- which.max(com_empty$csize == sort(com_empty$csize, decreasing=T)[2])

empty_c <- grid_11_20[row.names(grid_11_20@faces) %in%
                        names(
                          com_empty$membership[
                            (com_empty$membership==to_del_1) |
                              (com_empty$membership==to_del_2)
                          ]
                        )
]
faces3d(empty_c, col='red')

# Remove small isolated regions from the grid
final_11_20 <- grid_11_20[!(row.names(grid_11_20@faces) %in% row.names(empty_c@faces))]
faces3d(final_11_20, col='blue')

pos_final <- locate(final_11_20, dati_central)
t_final <- table(pos_final)

# Further remove small disconnected components
to_delete_2 <- final_11_20[
  (!(row.names(final_11_20@faces) %in% names(t_final))) |
    #(row.names(final_11_20@faces) %in% names(t_final[t_final==2])) |
    (row.names(final_11_20@faces) %in% names(t_final[t_final==1]))]

count_components(to_delete_2@graph)
to_del_com_2 <- components(to_delete_2@graph)
table(to_del_com_2$csize)

# Remove remaining small components
to_del.1 <- which.max(to_del_com_2$csize)
to_del.2 <- which.max(to_del_com_2$csize == sort(to_del_com_2$csize, decreasing=T)[2])

to_rem <- final_11_20[row.names(final_11_20@faces) %in%
                        names(
                          to_del_com_2$membership[
                            #(to_del_com$membership==to_del.2) |
                            (to_del_com_2$membership==to_del.1)
                          ]
                        )
]
faces3d(to_rem)

faces3d(final_11_20[!(row.names(final_11_20@faces) %in% row.names(to_rem@faces))])

# Identify and visualize the data points that should be removed
dati_to_rem <- all[which(locate(to_rem, all) != 'F0'),]
points3d(dati_to_rem[,1], dati_to_rem[,2], dati_to_rem[,3], col='red', size=2)

dens_rem <- as.vector(dens(dati_to_rem, to_rem, hns(dati_to_rem), row.names(to_rem@faces)))


# APPROACH 2 ###################################################################
gr_stim_copy <- gr_stim
gr_stim_copy[pdf_est <= 3] <- 0
table(gr_stim_copy)

gr_stim_copy[ENERGY <= 15000] <- 0
table(gr_stim_copy)
table(gr_stim_copy, gr)
adj.rand.index(gr_stim_copy, gr)

table(gr_stim[pos_final_component_2], gr[pos_final_component_2])
gr_stim_copy[ret_comp1@pdf <= 1.38] <- 0
table(gr_stim_copy[pos_final_component_2], gr[pos_final_component_2])
adj.rand.index(gr_stim_copy[pos_final_component_2], gr[pos_final_component_2])
adj.rand.index(gr_stim[pos_final_component_2], gr[pos_final_component_2])


# APPROACH 3 ###################################################################
range(clus_obs)
clus_obs_copy <- clus_obs
clus_obs_copy[clus_obs_copy == 'bck'] = 0
clus_obs_copy <- as.numeric(clus_obs_copy)
range(clus_obs_copy)

pdf_est_copy <- pdf_est

# Remove densities below the first quartile within each individual cluster
for (i in 1:max(clus_obs_copy)) {
  print(i)
  pos <- which(clus_obs_copy == i)
  dens <- pdf_est_copy[pos]
  soglia <- quantile(dens, 0.25)
  gr_stim_copy[pos][dens <= soglia/2] <- 0
}
table(gr_stim_copy, gr)
adj.rand.index(gr_stim_copy, gr)

# Further refine by removing densities below the first quartile, 
# considering only clusters inside `pos_poli`
for (i in 1:max(clus_obs_copy)) {
  print(i)
  pos <- which(clus_obs_copy == i)
  if (all(pos %in% pos_poli)) {
    dens <- pdf_est_copy[pos]
    soglia <- quantile(dens, 0.25)
    gr_stim_copy[pos][dens <= soglia] <- 0
  }
}
table(gr_stim_copy, gr)
adj.rand.index(gr_stim_copy, gr)

# Visualize specific data points in red
prova <- ret_poli[[4]]
dat <- prova@x
points3d(dat[,1], dat[,2], dat[,3], col='red')

ret_poli[[4]]@graph

# Remove densities below the first quartile in the central region, 
# considering only clusters inside `pos_final_central_3`
for (i in 1:max(clus_obs_copy)) {
  print(i)
  pos <- which(clus_obs_copy == i)
  if (all(pos %in% pos_final_central_3)) {
    dens <- pdf_est_copy[pos]
    soglia <- quantile(dens, 0.25)
    gr_stim_copy[pos][dens <= soglia] <- 0
  }
}

# Count the number of clusters in `ret_poli`
noc <- sapply(1:length(ret_poli), function(x) if (!is.null(ret_poli[[x]])) ret_poli[[x]]@noc)
table(unlist(noc))
noc

# Remove densities below the first quartile in each region
limit <- unique(com_poli$membership)
for (cc in limit) {
  print(cc)
  if (!is.null(ret_poli[[cc]])) {
    if (ret_poli[[cc]]@noc == 1) {
      nomi <- names(com_poli$membership[com_poli$membership==cc])
      pos <- which(posizioni9 %in% nomi)
      dati <- all[pos,]
      obs <- locate(poli, dati)
      dens <- ret_poli[[cc]]@pdf
      #dens <- sapply(1:dim(dati)[1], function(i) ret_poli[[cc]]@pdf[which(nomi == obs[i])])
      pos.bkg <- which(posizioni9.bkg %in% nomi)
      dati.bkg <- bkg[pos.bkg,]
      dens.bkg <- compute_dens(dati.bkg, poli, hns(dati), nomi)
      dens_diff <- dens - dens.bkg
      dens_diff <- sapply(1:dim(dati)[1], function(i) dens_diff[which(nomi == obs[i])])
      gr_stim_copy[pos][dens_diff <= 0] <- 0
    }
  }
}
table(gr_stim_copy[pos_poli], gr[pos_poli])
adj.rand.index(gr_stim_copy[pos_poli], gr[pos_poli])

# Remove densities below the first quartile in the central region
limit2 <- unique(com_central$membership)
for (cc in limit2) {
  print(cc)
  if (!is.null(ret_central[[cc]])) {
    if (ret_central[[cc]]@noc == 1) {
      nomi <- names(com_central$membership[com_central$membership==cc])
      pos <- which(posizioni11 %in% nomi)
      dati <- all[pos,]
      obs <- locate(final_central_3, dati)
      dens <- ret_central[[cc]]@pdf
      #dens <- sapply(1:dim(dati)[1], function(i) ret_central[[cc]]@pdf[which(nomi == obs[i])])
      pos.bkg <- which(posizioni11.bkg %in% nomi)
      dati.bkg <- bkg[pos.bkg,]
      dens.bkg <- compute_dens(dati.bkg, final_central_3, hns(dati), nomi)
      dens_diff <- dens - dens.bkg
      dens_diff <- sapply(1:dim(dati)[1], function(i) dens_diff[which(nomi == obs[i])])
      gr_stim_copy[pos][dens_diff <= 0] <- 0
    }
  }
}
table(gr_stim_copy[pos_final_central_3], gr[pos_final_central_3])
adj.rand.index(gr_stim_copy[pos_final_central_3], gr[pos_final_central_3])

# Compare classifications, excluding `pos_poli` and `pos_final_component_2`
table(gr_stim_copy[-c(pos_poli, pos_final_component_2)],
      gr[-c(pos_poli, pos_final_component_2)])
adj.rand.index(gr_stim_copy[-c(pos_poli, pos_final_component_2)],
               gr[-c(pos_poli, pos_final_component_2)])

adj.rand.index(gr_stim_copy, gr)
table(gr_stim_copy, gr)

# Compute density differences for the largest connected component
nomi <- row.names(final_component_2@faces)
pos <- which(posizioni11 %in% nomi)
dati <- all[pos,]
obs <- locate(final_component_2, dati)
dens <- ret_comp1@pdf
pos.bkg <- which(posizioni11.bkg %in% nomi)
dati.bkg <- bkg[pos.bkg,]
dens.bkg <- compute_dens(dati.bkg, final_component_2, hns(dati), nomi)
dens_diff <- dens - dens.bkg
dens_diff <- sapply(1:dim(dati)[1], function(i) dens_diff[which(nomi == obs[i])])
gr_stim_copy[pos][dens_diff <= 0] <- 0

table(gr_stim_copy[pos_final_component_2], gr[pos_final_component_2])
adj.rand.index(gr_stim_copy[pos_final_component_2], gr[pos_final_component_2])


# PLOT #########################################################################
# table(gr_stim[pos_poli], gr[pos_poli])
# table(gr_stim[-c(pos_poli, pos_final_component_2)], gr[-c(pos_poli, pos_final_component_2)])
# table(gr_stim[pos_final_central_3], gr[pos_final_central_3])
# table(gr_stim[pos_final_component_2], gr[pos_final_component_2])
# table(gr_stim[-c(pos_poli, pos_final_central_3, pos_final_component_2)],
#       gr[-c(pos_poli, pos_final_central_3, pos_final_component_2)])
# 
# plot3d(grid_11_20, guides=F, arcs=T, sphere=1, col='gray')
# points3d(sorg_final_central_3[,1], sorg_final_central_3[,2], 
#          sorg_final_central_3[,3], col='green', size=3)
# points3d(dati_final_central_3[,1], dati_final_central_3[,2],
#          dati_final_central_3[,3], col=gr_stim[pos_final_central_3]+1, size=3)
# 
# plot3d(grid_11_20, guides=F, arcs=T, sphere=1, col='gray')
# points3d(sorg_final_component_2[,1], sorg_final_component_2[,2], 
#          sorg_final_component_2[,3], col='green', size=3)
# points3d(dati_final_component_2[,1], dati_final_component_2[,2],
#          dati_final_component_2[,3], col=gr_stim[pos_final_component_2]+1, size=3)
# 
# plot3d(grid_9_20, guides=F, arcs=T, sphere=1, col='gray')
# points3d(sorg_poli[,1], sorg_poli[,2], 
#          sorg_poli[,3], col='green', size=3)
# points3d(dati_poli[,1], dati_poli[,2],
#          dati_poli[,3], col=gr_stim_2[pos_poli]+1, size=3)
# 
# faces3d(final_central_2)
# plot3d(grid_11_20, guides=F, arcs=T, sphere=1, col='gray')
# points3d(sorg_final_central_3[,1], sorg_final_central_3[,2], 
#          sorg_final_central_3[,3], col='red', size=3)
# points3d(dati_final_central_3[,1], dati_final_central_3[,2],
#          dati_final_central_3[,3], col='green', size=3)
# sorg_central_2 <- sorg[which(locate(central_2, sorg) != 'F0'),]
# points3d(sorg_final_central_2[,1], sorg_final_central_2[,2], 
#          sorg_final_central_2[,3], col='blue', size=3)
# 
# plot3d(grid_9_20, guides=F, arcs=T, sphere=1, col='gray')
# faces3d(poli, col='blue')
# points3d(sorg_poli[,1], sorg_poli[,2], 
#          sorg_poli[,3], col='green', size=3)
# points3d(dati_poli[,1], dati_poli[,2],
#          dati_poli[,3], col=gr_stim[pos_poli]+1, size=3)


# POTENTIAL SOLUTION ###########################################################
# # I need to evaluate the density in the cells when sources are absent
# bkg
# 
# # I need to consider the same mesh partition and the same smoothing parameter
# faces3d(poli)
# faces3d(final_central_3, col='blue')
# faces3d(final_component_2, col='red')
# 
# # For now, I will try only at the poles
# pdfCluster.data.2.sol <- function(x, mesh, h, nomi, x_bkg,
#                                   n.grid=min(round((5 + sqrt(length(nomi)))*4),
#                                              length(nomi)),
#                                   ... )
# {
#   call <- match.call()
#   x <- data.matrix(x)
#   x_bkg <- data.matrix(x_bkg)
#   if(any(sapply(x[1,],is.factor))) stop("All variables must be numeric")
#   if(any(is.na(x))) 
#   {warning(cat("NA in object",deparse(substitute(x)), "have been omitted", "\n"));
#     x <- na.omit(x)}
#   
#   # compute the elements to be used for both density estimation and connected components
#   l <- comp.centr.wei.ico.2(x, mesh, nomi)
#   t <- l$t
#   graph.nb <- l$adj # this is my graph
#   obs <- l$obs 
#   n <- l$n
#   
#   # density estimation
#   coarse_estimate <- do.call(vmf.bskde.2, list(x = t, h = h/10, n = n))
#   g <- exp(mean(log(coarse_estimate)))
#   h_variable <- h * (1/g * coarse_estimate)^(-0.5)
#   estimate <- do.call(vmf.bskde.2, list(x = t, h = as.vector(h_variable), n = n))
#   
#   # density estimation with background only
#   l_bkg <- comp.centr.wei.ico.2(x_bkg, mesh, nomi)
#   t_bkg <- l_bkg$t
#   graph.nb_bkg <- l_bkg$adj # this is my graph
#   obs_bkg <- l_bkg$obs 
#   n_bkg <- l_bkg$n
#   coarse_estimate_bkg <- do.call(vmf.bskde.2, list(x = t_bkg, h = hns(x_bkg)/10, n = n_bkg))
#   g_bkg <- exp(mean(log(coarse_estimate_bkg)))
#   h_variable <- hns(x_bkg) * (1/g_bkg * coarse_estimate_bkg)^(-0.5)
#   estimate_bkg <- do.call(vmf.bskde.2, list(x = t_bkg, h = as.vector(h_variable), n = n_bkg))
#   
#   # compute the difference
#   diff_est <- estimate - estimate_bkg
#   diff_est[diff_est < 0] = 0
#   
#   # check given arguments
#   N <- nrow(t)
#   if (n.grid > N) {
#     warning("n.grid too large, set equal to N")
#     n.grid <- min(n.grid, nrow(t))
#   }
#   
#   # connected components
#   nc <- pdfCluster:::num.con(t, diff_est, graph.nb,
#                              profile.grid = n.grid-2, correct=TRUE)  
#   struct <- pdfCluster:::con2tree(nc, diff_est)
#   if (struct$bad) {
#     message("No output given")
#   }
#   else {
#     g <- struct$g
#     g[struct$g == 0] <- NA
#     pdf.estim <- diff_est
#     names(l$adj) <- nomi
#     obs.cl <- sapply(1:dim(x)[1], function(i) g[which(names(l$adj) == obs[i])])
#     out <- new("pdfCluster2", call = call, x = data.matrix(x), 
#                pdf = pdf.estim,  nc = nc, graph = graph.nb, cluster.cores = g,
#                tree = struct$tree, noc = struct$noc, obs.cluster = obs.cl)
#     out
#   }
# }
# 
# # clus_obs.2 <- rep(0, dim(all)[1])
# # com_poli <- components(poli@graph)
# # ret_poli.2 <- vector(length=count_components(poli@graph), mode='list')
# # limit <- unique(com_poli$membership)
# # for (cc in limit) {
# #   print(cc)
# #   nomi <- names(com_poli$membership[com_poli$membership==cc])
# #   if (length(nomi) <= 2) {
# #     pos <- posizioni9 %in% nomi
# #     n <- sum(pos)
# #     if (n >= 5) {
# #       clus_obs.2[pos] <- max(clus_obs.2) + 1
# #     }
# #   }
# #   if (length(nomi) > 2) {
# #     pos <- which(posizioni9 %in% nomi)
# #     pos.bkg <- which(posizioni9.bkg %in% nomi)
# #     dati <- all[pos,]
# #     dati.bkg <- bkg[pos.bkg,]
# #     ret <- pdfCluster.data.2.sol(x=dati, x_bkg=dati.bkg, mesh=poli, h=hns(dati), nomi=nomi)
# #     ret_poli.2[[cc]] <- ret
# #     obs_ret <- ret@obs.cluster
# #     obs_ret[is.na(obs_ret)] = 0
# #     obs_ret[obs_ret != 0] <- obs_ret[obs_ret != 0] + max(clus_obs.2)
# #     clus_obs.2[pos] <- obs_ret
# #   }
# # }
# 
# primo <- ret_poli[[25820]]@x
# primo.mesh <- poli[row.names(poli@faces) %in%
#                      names(com_poli$membership[com_poli$membership==
#                                                  which(com_poli$csize==603)])]
# faces3d(primo.mesh)
# points3d(primo[,1], primo[,2], primo[,3], col='red')
# 
# pos.bkg <- which(posizioni9.bkg %in% row.names(primo.mesh@faces))
# 
# ret.sol <- pdfCluster.data.2.sol(x=primo, x_bkg=bkg[pos.bkg,], 
#                                  mesh=poli, h=hns(primo), 
#                                  nomi=row.names(primo.mesh@faces))
# ret.sol@pdf
# ret.sol@noc


# TPR AND FPR ##################################################################
# table(clus_obs[pos_prova_603], id_prova_603)
# 
# tab <- table(clus_obs[pos_prova_603], id_prova_603)
# tab <- tab[-nrow(tab), -ncol(tab)]
# tab <- t(tab)
# sol <- solve_LSAP(tab, maximum = T)
# mat <- matrix(0, nrow = nrow(tab), ncol = ncol(tab))
# mat[cbind(seq_along(sol), sol)] <- tab[cbind(seq_along(sol), sol)]
# mat
# # fpr
# sum(colSums(mat)==0)/ncol(mat)
# # tpr 
# sum(rowSums(mat)!=0)/nrow(mat)

tpr_fpr(clus_obs[pos_poli], id_poli)
tpr_fpr(clus_obs[pos_final_central_3], id_final_central_3)

sum(table(clus_obs_copy)<=4)
clus_min_4 <- names(table(clus_obs_copy)[table(clus_obs_copy) <= 4])
to_set_bck <- sapply(1:length(clus_obs_copy), 
                     function(i) clus_obs_copy[i] %in% clus_min_4)
clus_obs_copy[to_set_bck] <- 0
length(unique(clus_obs_copy))

tpr_fpr(clus_obs_copy[pos_poli], id_poli)
tpr_fpr(clus_obs_copy[pos_final_central_2], id_final_central_2)

tpr_fpr(clus_obs_copy, ID) # 0.88, 0.79

gr_stim_2 <- rep(0, length(clus_obs_copy))
gr_stim_2[clus_obs_copy > 0] = 1 
table(gr_stim_2)
adj.rand.index(gr_stim_2, gr)
table(gr_stim_2, gr)


# RESULTS ######################################################################
# Total
tpr_fpr(clus_obs, ID)
adj.rand.index(gr_stim, gr)
table(gr_stim, gr)
med(gr_stim, gr)

# poli
tpr_fpr(clus_obs[pos_poli], id_poli)
adj.rand.index(gr_stim[pos_poli], gr[pos_poli])
table(gr_stim[pos_poli], gr[pos_poli])
med(gr_stim[pos_poli], gr[pos_poli])

# Central
tpr_fpr(clus_obs[pos_final_central_3],
        ID[pos_final_central_3])
adj.rand.index(gr_stim[-c(pos_poli, pos_final_component_2)],
               gr[-c(pos_poli, pos_final_component_2)])
table(gr_stim[-c(pos_poli, pos_final_component_2)], gr[-c(pos_poli, pos_final_component_2)])
med(gr_stim[-c(pos_poli, pos_final_component_2)], gr[-c(pos_poli, pos_final_component_2)])

# Large connected component
tpr_fpr(clus_obs[pos_final_component_2],
        ID[pos_final_component_2], trasp=T)
adj.rand.index(gr_stim[pos_final_component_2], gr[pos_final_component_2])
table(gr_stim[pos_final_component_2], gr[pos_final_component_2])
med(gr_stim[pos_final_component_2], gr[pos_final_component_2])

# Remove small clusters
clus_obs_2 <- clus_obs
clus_obs_2[clus_obs_2 == 'bck'] = 0
clus_obs_2 <- as.numeric(clus_obs_2)
sum(table(clus_obs)<4)
clus_min_4 <- names(table(clus_obs_2)[table(clus_obs_2) < 4])
to_set_bck <- sapply(1:length(clus_obs_2), 
                     function(i) clus_obs_2[i] %in% clus_min_4)
clus_obs_2[to_set_bck] <- 0
length(unique(clus_obs_2))
gr_stim_2 <- rep(0, length(clus_obs_2))
gr_stim_2[clus_obs_2 > 0] = 1 

# Total
tpr_fpr(clus_obs_2, ID)
adj.rand.index(gr_stim_2, gr)
table(gr_stim_2, gr)
med(gr_stim_2, gr)

# poli
tpr_fpr(clus_obs_2[pos_poli], id_poli)
adj.rand.index(gr_stim_2[pos_poli], gr[pos_poli])
table(gr_stim_2[pos_poli], gr[pos_poli])
med(gr_stim_2[pos_poli], gr[pos_poli])

# Central
tpr_fpr(clus_obs_2[pos_final_central_3],
        ID[pos_final_central_3])
adj.rand.index(gr_stim_2[-c(pos_poli, pos_final_component_2)],
               gr[-c(pos_poli, pos_final_component_2)])
table(gr_stim_2[-c(pos_poli, pos_final_component_2)], gr[-c(pos_poli, pos_final_component_2)])
med(gr_stim_2[-c(pos_poli, pos_final_component_2)], gr[-c(pos_poli, pos_final_component_2)])

# Large component
tpr_fpr(clus_obs_2[pos_final_component_2],
        ID[pos_final_component_2], trasp=T)
adj.rand.index(gr_stim_2[pos_final_component_2], gr[pos_final_component_2])
table(gr_stim_2[pos_final_component_2], gr[pos_final_component_2])
med(gr_stim_2[pos_final_component_2], gr[pos_final_component_2])


# MED ##########################################################################
med(clus_obs, ID)
med(gr_stim, gr)
med(clus_obs_2, ID)
med(clus_obs[pos_poli], id_poli)
med(gr_stim[pos_poli], gr[pos_poli])
med(clus_obs[-c(pos_poli, pos_final_component_2)], ID[-c(pos_poli, pos_final_component_2)])
med(clus_obs[pos_final_component_2], ID[pos_final_component_2])


# POSSIBLE IMPROVEMENT #########################################################
bkg30 <- bckF_sim %>% arrange(TIME)
bkg70 <- bkg30[(round(nrow(bckF_sim)*0.3)+1):nrow(bckF_sim),]
bkg30 <- bkg30[1:round(nrow(bckF_sim)*0.3),]
bkg70 <- bkg70 %>% select(L, B) %>% gal2cart()
bkg30 <- bkg30 %>% select(L, B) %>% gal2cart()

all70 <- rbind(sorg, bkg70)

ID_70 <- c(alldata_new$ID,
           bckF_sim %>% arrange(TIME) %>% 
             slice((round(nrow(bckF_sim)*0.3)+1):nrow(bckF_sim)) %>% pull(ID))
pos9_70 <- locate(grid_9_20, all70)
pos11_70 <- locate(grid_11_20, all70)
pos9_70_bkg <- locate(grid_9_20, bkg70)
pos11_70_bkg <- locate(grid_11_20, bkg70)
pos9_30_bkg <- locate(grid_9_20, bkg30)
pos11_30_bkg <- locate(grid_11_20, bkg30)

# Mesh and partition
plot3d(grid_9_20, guides=F, arcs=T, sphere=1, col='grey')
points3d(all70[,1], all70[,2], all70[,3])

t_pos_9 <- table(pos9_70)
celle_vuote <- grid_9_20[!(row.names(grid_9_20@faces) %in% names(t_pos_9))]
faces3d(grid_9_20[row.names(grid_9_20@faces) %in% names(t_pos_9)], col='blue')
faces3d(celle_vuote, col='red')
is_connected(celle_vuote@graph)
count_components(celle_vuote@graph)
com_celle_vuote <- components(celle_vuote@graph)
table(com_celle_vuote$csize)
to_del_1 <- which.max(com_celle_vuote$csize)
to_del_2 <- which.max(com_celle_vuote$csize == sort(com_celle_vuote$csize, decreasing=T)[2])
empty_cells <- grid_9_20[row.names(grid_9_20@faces) %in%
                           names(
                             com_celle_vuote$membership[
                               (com_celle_vuote$membership==to_del_1) |
                                 (com_celle_vuote$membership==to_del_2)
                             ]
                           )
]
faces3d(empty_cells, col='red')
final <- grid_9_20[!(row.names(grid_9_20@faces) %in% row.names(empty_cells@faces))]
plot3d(final, guides=F, arcs=T, sphere=1)
faces3d(final, col='blue')
points3d(sorg[,1], sorg[,2], sorg[,3], col='red', size=3)
points3d(bkg70[,1], bkg70[,2], bkg70[,3], col='green', size=3)

com_final <- components(final@graph)
table(com_final$csize)
which.com <- which.max(com_final$csize)
central <- final[row.names(final@faces) %in%
                   names(
                     com_final$membership[
                       com_final$membership==which.com
                     ]
                   )
]
faces3d(central)

dati_central <- all70[which(locate(central, all70) != 'F0'),]
points3d(dati_central[,1], dati_central[,2], dati_central[,3], col='violet', size=1)

pos_central <- locate(central, dati_central)
t_central <- table(pos_central)

to_delete <- central[
  (!(row.names(central@faces) %in% names(t_central))) |
    (row.names(central@faces) %in% names(t_central[t_central==1])) |
    (row.names(central@faces) %in% names(t_central[t_central==2]))]
faces3d(to_delete, col='red')
faces3d(central[!(row.names(central@faces) %in% row.names(to_delete@faces))])

count_components(to_delete@graph)
to_del_com <- components(to_delete@graph)
table(to_del_com$csize)

to_del_com_1 <- which.max(to_del_com$csize)

to_del <- central[row.names(central@faces) %in%
                    names(
                      to_del_com$membership[
                        (to_del_com$membership==to_del_com_1)
                      ]
                    )
]
faces3d(to_del)
final_central <- central[!(row.names(central@faces) %in% row.names(to_del@faces))]

faces3d(final_central)
count_components(final_central@graph)
points3d(all70[which(locate(final_central, all70) != 'F0'),][,1], 
         all70[which(locate(final_central, all70) != 'F0'),][,2],
         all70[which(locate(final_central, all70) != 'F0'),][,3], col='red', size=3)

com_final_central <- components(final_central@graph)
table(com_final_central$csize)

final_component <- final_central[row.names(final_central@faces) %in% 
                                   names(com_final_central$membership[
                                     com_final_central$membership==
                                       which.max(com_final_central$csize)
                                   ])]
faces3d(final_component, col='blue')


faces3d(central, col='blue')
points3d(dati_central[,1], dati_central[,2], dati_central[,3], col='red', size=2)

plot3d(grid_11_20, guides=F, arcs=T, sphere=1, col='blue')

empty <- grid_11_20[!(occupied(grid_11_20, dati_central))@values]
com_empty <- components(empty@graph)
table(com_empty$csize)
faces3d(empty)

to_del_1 <- which.max(com_empty$csize)
to_del_2 <- which.max(com_empty$csize == sort(com_empty$csize, decreasing=T)[2])

empty_c <- grid_11_20[row.names(grid_11_20@faces) %in%
                        names(
                          com_empty$membership[
                            (com_empty$membership==to_del_1)
                          ]
                        )
]
faces3d(empty_c, col='red')

final_11_20 <- grid_11_20[!(row.names(grid_11_20@faces) %in% row.names(empty_c@faces))]
faces3d(final_11_20, col='blue')

pos_final <- locate(final_11_20, dati_central)
t_final <- table(pos_final)

to_delete_2 <- final_11_20[
  (!(row.names(final_11_20@faces) %in% names(t_final))) |
    (row.names(final_11_20@faces) %in% names(t_final[t_final==2])) |
    (row.names(final_11_20@faces) %in% names(t_final[t_final==1]))]

count_components(to_delete_2@graph)
to_del_com_2 <- components(to_delete_2@graph)
table(to_del_com_2$csize)

to_del.1 <- which.max(to_del_com_2$csize)
to_del.2 <- which.max(to_del_com_2$csize == sort(to_del_com_2$csize, decreasing=T)[2])

to_ren <- final_11_20[row.names(final_11_20@faces) %in%
                        names(
                          to_del_com_2$membership[
                            #(to_del_com$membership==to_del.2) |
                            (to_del_com_2$membership==to_del.1)
                          ]
                        )
]
faces3d(to_ren)
final_central_2 <- final_11_20[!(row.names(final_11_20@faces) %in% row.names(to_ren@faces))]
faces3d(final_central_2)
faces3d(final_component, col='blue')
com_final_central_2 <- components(final_central_2@graph)
table(com_final_central_2$csize)
final_component_2_70 <- final_central_2[row.names(final_central_2@faces) %in% 
                                          names(com_final_central_2$membership[
                                            com_final_central_2$membership==
                                              which.max(com_final_central_2$csize)
                                          ])]
faces3d(final_component_2_70, col='blue')
poli70 <- final[!(row.names(final@faces) %in% row.names(central@faces))]
faces3d(poli70)
comp.grande70 <- final_component_2_70
centr70 <- final_central_2[!(row.names(final_central_2@faces) %in% 
                               row.names(comp.grande70@faces))]
faces3d(centr70)

# Final regions
faces3d(poli70)
faces3d(centr70, col='blue')
faces3d(comp.grande70, col='red')

# Poles
pos_poli70 <- which(locate(poli70, all70) != 'F0')
id_poli70 <- ID_70[pos_poli70]
dati_poli70 <- all70[pos_poli70,]
sorg_poli70 <- sorg[which(locate(poli70, sorg) != 'F0'),]

# Central region excluding the largest connected component
pos_centr70 <- which(locate(centr70, all70) != 'F0')
id_centr70 <- ID_70[pos_centr70]
dati_centr70 <- all70[pos_centr70,]
sorg_centr70 <- sorg[which(locate(centr70, sorg) != 'F0'),]

# Largest connected component
pos_comp.grande70 <- which(locate(comp.grande70, all70) != 'F0')
id_comp.grande70 <- ID_70[pos_comp.grande70]
dati_comp.grande70 <- all70[pos_comp.grande70,]
sorg_comp.grande70 <- sorg[which(locate(comp.grande70, sorg) != 'F0'),]

# Clustering
clus_obs <- rep(0, dim(all70)[1])

# First, I classify the observations in the poles
com_poli70 <- components(poli70@graph)
ret_poli70 <- vector(length=count_components(poli70@graph), mode='list')
limit <- unique(com_poli70$membership)
for (cc in limit) {
  print(cc)
  nomi <- names(com_poli70$membership[com_poli70$membership==cc])
  if (length(nomi) <= 2) {
    pos <- pos9_70 %in% nomi
    n <- sum(pos)
    if (n >= 5) {
      clus_obs[pos] <- max(clus_obs) + 1
    }
  }
  if (length(nomi) > 2) {
    pos <- which(pos9_70 %in% nomi)
    dati <- all70[pos,]
    ret <- pdfCluster.data.2.new.2(x=dati, mesh=poli70, h=hns(dati)*1.5, nomi=nomi)
    ret_poli70[[cc]] <- ret
    obs_ret <- ret@obs.cluster
    obs_ret[is.na(obs_ret)] = 0
    obs_ret[obs_ret != 0] <- obs_ret[obs_ret != 0] + max(clus_obs)
    clus_obs[pos] <- obs_ret
  }
}

# I classify the observations in the central area, excluding the largest component
com_central70 <- components(centr70@graph)
ret_central70 <- vector(length=count_components(centr70@graph), mode='list')
limit2 <- unique(com_central70$membership)
for (cc in limit2) {
  print(cc)
  nomi <- names(com_central70$membership[com_central70$membership==cc])
  if (length(nomi) <= 2) {
    pos <- pos11_70 %in% nomi
    n <- sum(pos)
    if (n >= 5) clus_obs[pos] <- max(clus_obs) + 1
  }
  if (length(nomi) > 2) {
    pos <- which(pos11_70 %in% nomi)
    dati <- all70[pos,]
    ret <- pdfCluster.data.2.new.2(x=dati, mesh=centr70, h=hns(dati)*1.5, nomi=nomi)
    ret_central70[[cc]] <- ret
    obs_ret <- ret@obs.cluster
    obs_ret[is.na(obs_ret)] = 0
    obs_ret[obs_ret != 0] <- obs_ret[obs_ret != 0] + max(clus_obs)
    clus_obs[pos] <- obs_ret
  }
}

# I classify the observations in the largest component
nomi_comp <- row.names(comp.grande70@faces)
ret_comp <- pdfCluster.data.2.new.2(x=dati_comp.grande70,
                                    mesh=comp.grande70, 
                                    h=hns(dati_comp.grande70), 
                                    nomi=nomi_comp)
obs_ret_comp <- ret_comp@obs.cluster
obs_ret_comp[is.na(obs_ret_comp)] = 0
clus_obs[pos_comp.grande70] <- obs_ret_comp

# Results
obs_ret_comp[obs_ret_comp != 0] <- obs_ret_comp[obs_ret_comp != 0] + max(clus_obs)
clus_obs[pos_comp.grande70] <- obs_ret_comp
clus_obs[clus_obs==0] <- 'bck'

gr_stim <- ifelse(clus_obs == 'bck', 0, 1)
table(gr_stim)
gr <- ifelse(ID_70 == 'bck', 0, 1)
table(gr)
adj.rand.index(gr_stim, gr)
table(gr_stim, gr)

table(gr_stim[pos_poli70], gr[pos_poli70])
table(gr_stim[-c(pos_poli70, pos_comp.grande70)], gr[-c(pos_poli70, pos_comp.grande70)])
table(gr_stim[pos_centr70], gr[pos_centr70])
table(gr_stim[pos_comp.grande70], gr[pos_comp.grande70])
table(gr_stim[-c(pos_poli70, pos_centr70, pos_comp.grande70)],
      gr[-c(pos_poli70, pos_centr70, pos_comp.grande70)])

# plot3d(grid_11_20, guides=F, arcs=T, sphere=1, col='gray')
# faces3d(centr70, col='blue')
# points3d(sorg_centr70[,1], sorg_centr70[,2], 
#          sorg_centr70[,3], col='green', size=3)
# points3d(dati_centr70[,1], dati_centr70[,2],
#          dati_centr70[,3], col=gr_stim_2[pos_centr70]+1, size=3)
# 
# plot3d(grid_11_20, guides=F, arcs=T, sphere=1, col='gray')
# points3d(sorg_comp.grande70[,1], sorg_comp.grande70[,2], 
#          sorg_comp.grande70[,3], col='green', size=3)
# points3d(dati_comp.grande70[,1], dati_comp.grande70[,2],
#          dati_comp.grande70[,3], col=gr_stim_2[pos_comp.grande70]+1, size=3)
# 
# plot3d(grid_9_20, guides=F, arcs=T, sphere=1, col='gray')
# faces3d(poli70, col='blue')
# points3d(sorg_poli70[,1], sorg_poli70[,2], 
#          sorg_poli70[,3], col='green', size=3)
# points3d(dati_poli70[,1], dati_poli70[,2],
#          dati_poli70[,3], col=gr_stim_2[pos_poli70]+1, size=3)

# Correction to improve the results - 1
clus_obs_copy <- clus_obs
clus_obs_copy[clus_obs_copy == 'bck'] = 0
clus_obs_copy <- as.numeric(clus_obs_copy)
range(clus_obs_copy)

# Poles
gr_stim_copy <- gr_stim
limit <- unique(com_poli70$membership)
for (cc in limit) {
  print(cc)
  if (!is.null(ret_poli70[[cc]])) {
    #if (ret_poli70[[cc]]@noc == 1) {
    nomi <- names(com_poli70$membership[com_poli70$membership==cc])
    pos <- which(pos9_70 %in% nomi)
    dati <- all70[pos,]
    obs <- locate(poli70, dati)
    dens <- ret_poli70[[cc]]@pdf
    pos.bkg <- which(pos9_30_bkg %in% nomi)
    dati.bkg <- bkg30[pos.bkg,]
    if (length(dati.bkg) == 0) next
    if (length(dati.bkg) == 3) dati.bkg <- t(as.matrix(dati.bkg))
    dens.bkg <- compute_dens(dati.bkg, poli70, hns(dati)*1.5, nomi)
    dens_diff <- dens - dens.bkg
    dens_diff <- sapply(1:dim(dati)[1], function(i) dens_diff[which(nomi == obs[i])])
    gr_stim_copy[pos][dens_diff <= 0] <- 0
    clus_obs_copy[pos][dens_diff <= 0] <- 0
    #}
  }
}
table(gr_stim_copy[pos_poli70], gr[pos_poli70])
adj.rand.index(gr_stim_copy[pos_poli70], gr[pos_poli70])
table(gr_stim[pos_poli70], gr[pos_poli70])
adj.rand.index(gr_stim[pos_poli70], gr[pos_poli70])

# Central region excluding the largest connected component
limit2 <- unique(com_central70$membership)
for (cc in limit2) {
  print(cc)
  if (!is.null(ret_central70[[cc]])) {
    #if (ret_central70[[cc]]@noc == 1) {
    nomi <- names(com_central70$membership[com_central70$membership==cc])
    pos <- which(pos11_70 %in% nomi)
    dati <- all70[pos,]
    obs <- locate(centr70, dati)
    dens <- ret_central70[[cc]]@pdf
    pos.bkg <- which(pos11_30_bkg %in% nomi)
    dati.bkg <- bkg30[pos.bkg,]
    if (length(dati.bkg) == 0) next
    if (length(dati.bkg) == 3) dati.bkg <- t(as.matrix(dati.bkg))
    dens.bkg <- compute_dens(dati.bkg, centr70, hns(dati)*1.5, nomi)
    dens_diff <- dens - dens.bkg
    dens_diff <- sapply(1:dim(dati)[1], function(i) dens_diff[which(nomi == obs[i])])
    gr_stim_copy[pos][dens_diff <= 0] <- 0
    clus_obs_copy[pos][dens_diff <= 0] <- 0
    #}
  }
}
table(gr_stim_copy[pos_centr70], gr[pos_centr70])
adj.rand.index(gr_stim_copy[pos_centr70], gr[pos_centr70])
table(gr_stim_copy[-c(pos_poli70, pos_comp.grande70)],
      gr[-c(pos_poli70, pos_comp.grande70)])
adj.rand.index(gr_stim_copy[-c(pos_poli70, pos_comp.grande70)],
               gr[-c(pos_poli70, pos_comp.grande70)])

table(gr_stim[pos_centr70], gr[pos_centr70])
adj.rand.index(gr_stim[pos_centr70], gr[pos_centr70])
table(gr_stim[-c(pos_poli70, pos_comp.grande70)],
      gr[-c(pos_poli70, pos_comp.grande70)])
adj.rand.index(gr_stim[-c(pos_poli70, pos_comp.grande70)],
               gr[-c(pos_poli70, pos_comp.grande70)])

# Largest connected component
nomi <- row.names(comp.grande70@faces)
pos <- which(pos11_70 %in% nomi)
dati <- all70[pos,]
obs <- locate(comp.grande70, dati)
dens <- ret_comp@pdf
pos.bkg <- which(pos11_30_bkg %in% nomi)
dati.bkg <- bkg[pos.bkg,]
dens.bkg <- compute_dens(dati.bkg, comp.grande70, hns(dati), nomi)
dens_diff <- dens - dens.bkg
dens_diff <- sapply(1:dim(dati)[1], function(i) dens_diff[which(nomi == obs[i])])
gr_stim_copy[pos][dens_diff <= 0] <- 0

table(gr_stim_copy[pos_comp.grande70], gr[pos_comp.grande70])
adj.rand.index(gr_stim_copy[pos_comp.grande70], gr[pos_comp.grande70])

# TPR and FPR
tpr_fpr(clus_obs[pos_poli70], id_poli70)
tpr_fpr(clus_obs[pos_centr70], id_centr70)

clus_obs_copy_2 <- clus_obs
clus_obs_copy_2[clus_obs_copy_2 == 'bck'] = 0
clus_obs_copy_2 <- as.numeric(clus_obs_copy_2)
sum(table(clus_obs_copy_2)<=4)
clus_min_4 <- names(table(clus_obs_copy_2)[table(clus_obs_copy_2) <= 4])
to_set_bck <- sapply(1:length(clus_obs_copy_2), 
                     function(i) clus_obs_copy_2[i] %in% clus_min_4)
clus_obs_copy_2[to_set_bck] <- 0
length(unique(clus_obs_copy_2))

tpr_fpr(clus_obs_copy_2[pos_poli70], id_poli70)
adj.rand.index(clus_obs_copy_2[pos_poli70], id_poli70)
tpr_fpr(clus_obs_copy_2[pos_centr70], id_centr70)

obs_ret_comp <- ret_comp@obs.cluster
obs_ret_comp[is.na(obs_ret_comp)] = 0
obs_ret_comp[obs_ret_comp != 0] = obs_ret_comp[obs_ret_comp != 0] + max(clus_obs_copy_2)
clus_obs_copy_2[pos_comp.grande70] <- obs_ret_comp
clus_obs[pos_comp.grande70] <- obs_ret_comp

tpr_fpr(id_comp.grande70, clus_obs_copy_2[pos_comp.grande70])

gr_stim_2 <- rep(0, length(clus_obs_copy_2))
gr_stim_2[clus_obs_copy_2 > 0] = 1 
adj.rand.index(gr_stim_2[pos_poli70], gr[pos_poli70])
table(gr_stim_2[pos_poli70], gr[pos_poli70])
adj.rand.index(gr_stim_2[-c(pos_poli70, pos_comp.grande70)],
               gr[-c(pos_poli70, pos_comp.grande70)])
table(gr_stim_2[-c(pos_poli70, pos_comp.grande70)],
      gr[-c(pos_poli70, pos_comp.grande70)])
adj.rand.index(gr_stim_2[pos_comp.grande70], gr[pos_comp.grande70])
table(gr_stim_2[pos_comp.grande70], gr[pos_comp.grande70])

tpr_fpr(clus_obs_copy_2, ID_70)
table(gr_stim_2, gr)
adj.rand.index(gr_stim_2, gr)

# Correction to improve the results - 2
# Poli
gr_stim_copy_2 <- gr_stim_2
clus_obs_copy_3 <- clus_obs_copy_2
limit <- unique(com_poli70$membership)
for (cc in limit) {
  print(cc)
  if (!is.null(ret_poli70[[cc]])) {
    #if (ret_poli70[[cc]]@noc == 1) {
    nomi <- names(com_poli70$membership[com_poli70$membership==cc])
    pos <- which(pos9_70 %in% nomi)
    dati <- all70[pos,]
    obs <- locate(poli70, dati)
    dens <- ret_poli70[[cc]]@pdf
    pos.bkg <- which(pos9_30_bkg %in% nomi)
    dati.bkg <- bkg30[pos.bkg,]
    if (length(dati.bkg) == 0) next
    if (length(dati.bkg) == 3) dati.bkg <- t(as.matrix(dati.bkg))
    dens.bkg <- compute_dens(dati.bkg, poli70, hns(dati)*1.5, nomi)
    dens_diff <- dens - dens.bkg
    dens_diff <- sapply(1:dim(dati)[1], function(i) dens_diff[which(nomi == obs[i])])
    gr_stim_copy_2[pos][dens_diff <= 0] <- 0
    clus_obs_copy_3[pos][dens_diff <= 0] <- 0
    #}
  }
}
table(gr_stim_copy_2[pos_poli70], gr[pos_poli70])
adj.rand.index(gr_stim_copy_2[pos_poli70], gr[pos_poli70])
table(gr_stim[pos_poli70], gr[pos_poli70])
adj.rand.index(gr_stim[pos_poli70], gr[pos_poli70])
tpr_fpr(clus_obs_copy_3[pos_poli70], id_poli70)


# Central region excluding the largest connected component
limit2 <- unique(com_central70$membership)
for (cc in limit2) {
  print(cc)
  if (!is.null(ret_central70[[cc]])) {
    #if (ret_central70[[cc]]@noc == 1) {
    nomi <- names(com_central70$membership[com_central70$membership==cc])
    pos <- which(pos11_70 %in% nomi)
    dati <- all70[pos,]
    obs <- locate(centr70, dati)
    dens <- ret_central70[[cc]]@pdf
    pos.bkg <- which(pos11_30_bkg %in% nomi)
    dati.bkg <- bkg30[pos.bkg,]
    if (length(dati.bkg) == 0) next
    if (length(dati.bkg) == 3) dati.bkg <- t(as.matrix(dati.bkg))
    dens.bkg <- compute_dens(dati.bkg, centr70, hns(dati)*1.5, nomi)
    dens_diff <- dens - dens.bkg
    dens_diff <- sapply(1:dim(dati)[1], function(i) dens_diff[which(nomi == obs[i])])
    gr_stim_copy_2[pos][dens_diff <= 0] <- 0
    clus_obs_copy_3[pos][dens_diff <= 0] <- 0
    #}
  }
}
table(gr_stim_copy_2[pos_centr70], gr[pos_centr70])
adj.rand.index(gr_stim_copy_2[pos_centr70], gr[pos_centr70])
table(gr_stim_copy_2[-c(pos_poli70, pos_comp.grande70)],
      gr[-c(pos_poli70, pos_comp.grande70)])
adj.rand.index(gr_stim_copy_2[-c(pos_poli70, pos_comp.grande70)],
               gr[-c(pos_poli70, pos_comp.grande70)])
tpr_fpr(clus_obs_copy_3[pos_centr70], id_centr70)

table(gr_stim[pos_centr70], gr[pos_centr70])
adj.rand.index(gr_stim[pos_centr70], gr[pos_centr70])
table(gr_stim[-c(pos_poli70, pos_comp.grande70)],
      gr[-c(pos_poli70, pos_comp.grande70)])
adj.rand.index(gr_stim[-c(pos_poli70, pos_comp.grande70)],
               gr[-c(pos_poli70, pos_comp.grande70)])

# Largest connected component
nomi <- row.names(comp.grande70@faces)
pos <- which(pos11_70 %in% nomi)
gr_stim_copy_2[pos][dens_diff <= 0] <- 0
clus_obs_copy_3[pos][dens_diff <= 0] <- 0

tpr_fpr(id_comp.grande70, clus_obs_copy_3[pos_comp.grande70])
adj.rand.index(gr_stim_copy_2[pos_comp.grande70], gr[pos_comp.grande70])
table(gr_stim_copy_2[pos_comp.grande70], gr[pos_comp.grande70])

# Final results
clus_obs[clus_obs == 'bck'] <- 0
clus_obs <- as.numeric(clus_obs)
clus_obs[pos_comp.grande70] <- clus_obs_copy_2[pos_comp.grande70]
clus_obs_copy[pos_comp.grande70] <- clus_obs_copy_3[pos_comp.grande70]
gr_stim[pos_comp.grande70] <- gr_stim_2[pos_comp.grande70]
gr_stim_copy[pos_comp.grande70] <- gr_stim_copy_2[pos_comp.grande70]

# Total
tpr_fpr(clus_obs, ID_70)
adj.rand.index(gr_stim, gr)
table(gr_stim, gr)
med(gr_stim, gr)

# Poles
tpr_fpr(clus_obs[pos_poli70], id_poli70)
adj.rand.index(gr_stim[pos_poli70], gr[pos_poli70])
table(gr_stim[pos_poli70], gr[pos_poli70])

# Central
tpr_fpr(clus_obs[pos_centr70], ID_70[pos_centr70])
adj.rand.index(gr_stim[-c(pos_poli70, pos_comp.grande70)],
               gr[-c(pos_poli70, pos_comp.grande70)])
table(gr_stim[-c(pos_poli70, pos_comp.grande70)], gr[-c(pos_poli70, pos_comp.grande70)])

# Large component
tpr_fpr(id_comp.grande70, clus_obs[pos_comp.grande70], trasp=T)
adj.rand.index(gr_stim[pos_comp.grande70], gr[pos_comp.grande70])
table(gr_stim[pos_comp.grande70], gr[pos_comp.grande70])

# Total corrected
tpr_fpr(clus_obs_copy, ID_70)
adj.rand.index(gr_stim_copy, gr)
table(gr_stim_copy, gr)

# Poles corrected
tpr_fpr(clus_obs_copy[pos_poli70], id_poli70)
adj.rand.index(gr_stim_copy[pos_poli70], gr[pos_poli70])
table(gr_stim_copy[pos_poli70], gr[pos_poli70])

# Central corrected
tpr_fpr(clus_obs_copy[pos_centr70], ID_70[pos_centr70])
adj.rand.index(gr_stim_copy[-c(pos_poli70, pos_comp.grande70)],
               gr[-c(pos_poli70, pos_comp.grande70)])
table(gr_stim_copy[-c(pos_poli70, pos_comp.grande70)], gr[-c(pos_poli70, pos_comp.grande70)])

# Large component corrected
tpr_fpr(id_comp.grande70, clus_obs_copy[pos_comp.grande70], trasp=T)
adj.rand.index(gr_stim_copy[pos_comp.grande70], gr[pos_comp.grande70])
table(gr_stim_copy[pos_comp.grande70], gr[pos_comp.grande70])

# Remove the small clusters
clus_obs_copy_v2 <- clus_obs
clus_obs_copy_v2[clus_obs_copy_v2 == 'bck'] = 0
clus_obs_copy_v2 <- as.numeric(clus_obs_copy_v2)
sum(table(clus_obs_copy_v2)<4)
clus_min_4 <- names(table(clus_obs_copy_v2)[table(clus_obs_copy_v2) < 4])
to_set_bck <- sapply(1:length(clus_obs_copy_v2), 
                     function(i) clus_obs_copy_v2[i] %in% clus_min_4)
clus_obs_copy_v2[to_set_bck] <- 0
length(unique(clus_obs_copy_v2))
gr_stim_v2 <- rep(0, length(clus_obs_copy_v2))
gr_stim_v2[clus_obs_copy_v2 > 0] = 1 

# Total minus small clusters
tpr_fpr(clus_obs_copy_v2, ID_70)
adj.rand.index(gr_stim_v2, gr)
table(gr_stim_v2, gr)
med(gr_stim_v2, gr)

# Poles minus small clusters
tpr_fpr(clus_obs_copy_v2[pos_poli70], id_poli70)
adj.rand.index(gr_stim_v2[pos_poli70], gr[pos_poli70])
table(gr_stim_v2[pos_poli70], gr[pos_poli70])
med(gr_stim_v2[pos_poli70], gr[pos_poli70])

# Central minus small clusters
tpr_fpr(clus_obs_copy_v2[pos_centr70], ID_70[pos_centr70])
adj.rand.index(gr_stim_v2[-c(pos_poli70, pos_comp.grande70)],
               gr[-c(pos_poli70, pos_comp.grande70)])
table(gr_stim_v2[-c(pos_poli70, pos_comp.grande70)], gr[-c(pos_poli70, pos_comp.grande70)])
med(gr_stim_v2[-c(pos_poli70, pos_comp.grande70)], gr[-c(pos_poli70, pos_comp.grande70)])

# Large component minus small clusters
tpr_fpr(id_comp.grande70, clus_obs_copy_v2[pos_comp.grande70], trasp=T)
adj.rand.index(gr_stim_v2[pos_comp.grande70], gr[pos_comp.grande70])
table(gr_stim_v2[pos_comp.grande70], gr[pos_comp.grande70])
med(gr_stim_v2[pos_comp.grande70], gr[pos_comp.grande70])

# Correction to improve the results - 3
gr_stim_copy_v2 <- gr_stim_v2
gr_stim_copy_v2[pos_comp.grande70] <- gr_stim_copy_2[pos_comp.grande70]
clus_obs_copy_v3 <- clus_obs_copy_v2
clus_obs_copy_v3[pos_comp.grande70] <- clus_obs_copy_3[pos_comp.grande70]
limit <- unique(com_poli70$membership)
for (cc in limit) {
  print(cc)
  if (!is.null(ret_poli70[[cc]])) {
    #if (ret_poli70[[cc]]@noc == 1) {
    nomi <- names(com_poli70$membership[com_poli70$membership==cc])
    pos <- which(pos9_70 %in% nomi)
    dati <- all70[pos,]
    obs <- locate(poli70, dati)
    dens <- ret_poli70[[cc]]@pdf
    pos.bkg <- which(pos9_30_bkg %in% nomi)
    dati.bkg <- bkg30[pos.bkg,]
    if (length(dati.bkg) == 0) next
    if (length(dati.bkg) == 3) dati.bkg <- t(as.matrix(dati.bkg))
    dens.bkg <- compute_dens(dati.bkg, poli70, hns(dati)*1.5, nomi)
    dens_diff <- dens - dens.bkg
    dens_diff <- sapply(1:dim(dati)[1], function(i) dens_diff[which(nomi == obs[i])])
    gr_stim_copy_v2[pos][dens_diff <= 0] <- 0
    clus_obs_copy_v3[pos][dens_diff <= 0] <- 0
    #}
  }
}

# Central region excluding the largest connected component
limit2 <- unique(com_central70$membership)
for (cc in limit2) {
  print(cc)
  if (!is.null(ret_central70[[cc]])) {
    #if (ret_central70[[cc]]@noc == 1) {
    nomi <- names(com_central70$membership[com_central70$membership==cc])
    pos <- which(pos11_70 %in% nomi)
    dati <- all70[pos,]
    obs <- locate(centr70, dati)
    dens <- ret_central70[[cc]]@pdf
    pos.bkg <- which(pos11_30_bkg %in% nomi)
    dati.bkg <- bkg30[pos.bkg,]
    if (length(dati.bkg) == 0) next
    if (length(dati.bkg) == 3) dati.bkg <- t(as.matrix(dati.bkg))
    dens.bkg <- compute_dens(dati.bkg, centr70, hns(dati)*1.5, nomi)
    dens_diff <- dens - dens.bkg
    dens_diff <- sapply(1:dim(dati)[1], function(i) dens_diff[which(nomi == obs[i])])
    gr_stim_copy_v2[pos][dens_diff <= 0] <- 0
    clus_obs_copy_v3[pos][dens_diff <= 0] <- 0
    #}
  }
}

# Total corrected minus small clusters
tpr_fpr(clus_obs_copy_v3, ID_70)
adj.rand.index(gr_stim_copy_v2, gr)
table(gr_stim_copy_v2, gr)
med(gr_stim_copy_v2, gr)

# Poles corrected minus small clusters
tpr_fpr(clus_obs_copy_v3[pos_poli70], id_poli70)
adj.rand.index(gr_stim_copy_v2[pos_poli70], gr[pos_poli70])
table(gr_stim_copy_v2[pos_poli70], gr[pos_poli70])
med(gr_stim_copy_v2[pos_poli70], gr[pos_poli70])

# Central corrected minus small clusters
tpr_fpr(clus_obs_copy_v3[pos_centr70], ID_70[pos_centr70])
adj.rand.index(gr_stim_copy_v2[-c(pos_poli70, pos_comp.grande70)],
               gr[-c(pos_poli70, pos_comp.grande70)])
table(gr_stim_copy_v2[-c(pos_poli70, pos_comp.grande70)], gr[-c(pos_poli70, pos_comp.grande70)])
med(gr_stim_copy_v2[-c(pos_poli70, pos_comp.grande70)], gr[-c(pos_poli70, pos_comp.grande70)])

# Large component corrected minus small clusters
tpr_fpr(id_comp.grande70, clus_obs_copy_v3[pos_comp.grande70], trasp=T)
adj.rand.index(gr_stim_copy_v2[pos_comp.grande70], gr[pos_comp.grande70])
table(gr_stim_copy_v2[pos_comp.grande70], gr[pos_comp.grande70])
med(gr_stim_copy_v2[pos_comp.grande70], gr[pos_comp.grande70])

