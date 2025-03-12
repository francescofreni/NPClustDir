# LOADING LIBRARIES ############################################################
library(rcosmo)
library(pdfCluster)
library(mvtnorm)
library(icosa)
library(rgl)
library(igraph)
library(tidyverse)
library(clue)
library(ks)
library(Directional)
source('functions.R')


# FERMI DATA ###################################################################
sorg <- gal2cart(alldata_new[,c('L', 'B')])
bkg <- gal2cart(bckF_sim[,c('L', 'B')])
all <- rbind(sorg, bkg)
id_sorg <- alldata_new[,'ID']
id <- c(id_sorg, rep(0, nrow(bckF_sim)))

# I try to apply the method only to a portion of the sphere  
# window <- win %>% data.frame %>% pol2cart %>% cart2gal  

# I select both the coordinates and the IDs, which will be useful for
# evaluating the results  
sorg_gal <- cart2gal(sorg)
pos_s <- which(sorg_gal[,1] <= 30 & sorg_gal[,1] >= -10 & 
             sorg_gal[,2] <= 55 & sorg_gal[,2] >= 40)
subsec_s <-  sorg %>% cart2gal %>% as_tibble %>% 
  filter(l <= 30 &
           l >= -10 &
           b <= 55 &
           b >= 40) %>% 
  as.data.frame %>% gal2cart
id_subsec_s <- id_sorg[pos_s]

bkg_gal <- cart2gal(bkg)
pos_b <- which(bkg_gal[,1] <= 30 & bkg_gal[,1] >= -10 & 
               bkg_gal[,2] <= 55 & bkg_gal[,2] >= 40)
subsec_b <-  bkg %>% cart2gal %>% as_tibble %>% 
  filter(l <= 30 &
           l >= -10 &
           b <= 55 &
           b >= 40) %>% 
  as.data.frame %>% gal2cart
id_subsec_b <- rep(0, nrow(subsec_b))

subsec <- rbind(subsec_s, subsec_b)
id_subsec <- c(id_subsec_s, id_subsec_b)
points3d(subsec_s[,1], subsec_s[,2], subsec_s[,3], col=2)
points3d(subsec_b[,1], subsec_b[,2], subsec_b[,3])
#points3d(all[,1], all[,2], all[,3], col=2)
#displayPixelBoundaries(128)

# I want to use an nside of 128 for now. The problem is how to section  
# the sphere to obtain only the cells of interest.  
# There are two ways to do this: using coordinates or using the graph.  
# For now, I prefer to use the graph, which is the slower method but the one  
# I will use later.  

# In this case, the grid size used for the procedure is the same as the one  
# used for partitioning the grid. I could also perform the partitioning  
# with this size and then group the regions. For now, I want to try  
# a coarser grid to see if it is sufficient near the poles. 

ns2 <- 128
j2 <- 7
cmbdf4 <- CMBDataFrame(nside = ns2, ordering = "nested", coords = 'cartesian')
cmbdf4$pix <- 1:nrow(cmbdf4)
g.128 <- build.graph.v2(subsec, cmbdf4, j2)

# Compute the number of observations per cell 
loc <- localize(subsec, j2)
intensity <- rep(0, 12*ns2^2)
tab <- table(loc)
intensity[as.numeric(names(tab))] <- tab
intensity <- intensity[cmbdf4$pix]
cmbdf4$I <- intensity

# Insert them into the graph, considering the correct position, because  
# the pixels in the graph are not necessarily ordered.  
V(g.128)$intensity <- intensity[as.numeric(names(V(g.128)))]

empty_cells <- delete_vertices(g.128, V(g.128)$intensity != 0)
count_components(empty_cells)
comp <- components(empty_cells)
to_del <- names(comp$membership[comp$membership==1])
occupied <- delete_vertices(g.128, to_del)
cmbdf.win4 <- cmbdf4 %>% filter(pix %in% V(occupied)$name)

points3d(subsec_s[,1], subsec_s[,2], subsec_s[,3], col=2)
points3d(subsec_b[,1], subsec_b[,2], subsec_b[,3])
points3d(cmbdf.win4$x, cmbdf.win4$y, cmbdf.win4$z, col = 'purple', cex = 2)
# Maybe I could also try with nside 64 to subdivide the sphere and then  
# perform clustering with nside 128  

# Clustering test 
test.128 <- pdfCluster.healpix(subsec, cmbdf.win4, j2, 0.1)
plot(test.128@tree, center=T)
cl.128 <- test.64@obs.cluster
cl.128[is.na(cl.128)] = 0
table(cl.128, id_subsec)

# This gives an error because the detected region consists of multiple  
# connected components, so I should apply the clustering algorithm  
# to each connected component separately. To prevent multiple connected  
# components from forming in this example, I calculate the connected  
# components with nside = 64, remove the empty cells, and perform  
# clustering with a tessellation of nside = 128.  
cmbdf5 <- CMBDataFrame(nside = 64, ordering = "nested", coords = 'cartesian')
cmbdf5$pix <- 1:nrow(cmbdf5)
g.64 <- build.graph.v2(subsec, cmbdf5, 6)

loc <- localize(subsec, 6)
intensity <- rep(0, 12*64^2)
tab <- table(loc)
intensity[as.numeric(names(tab))] <- tab
intensity <- intensity[cmbdf5$pix]
cmbdf5$I <- intensity

V(g.64)$intensity <- intensity[as.numeric(names(V(g.64)))]
empty_cells <- delete_vertices(g.64, V(g.64)$intensity != 0)
count_components(empty_cells)
comp <- components(empty_cells)
to_del <- names(comp$membership[comp$membership==1])
occupied <- delete_vertices(g.64, to_del)

kids <- children(as.numeric(V(occupied)$name))
cmbdf.win5 <- cmbdf4 %>% filter(pix %in% kids)

points3d(subsec_s[,1], subsec_s[,2], subsec_s[,3], col=2)
points3d(subsec_b[,1], subsec_b[,2], subsec_b[,3])
points3d(cmbdf.win5$x, cmbdf.win5$y, cmbdf.win5$z, col = 'purple', cex = 2)

# Clustering 
test.128 <- pdfCluster.healpix(subsec, cmbdf.win5, j2, 0.02)
plot(test.128@tree, center=T)
cl.128 <- test.128@obs.cluster
cl.128[is.na(cl.128)] = 0
table(cl.128, id_subsec)

# Graphical evaluation  
points3d(subsec[,1], subsec[,2], subsec[,3], col=cl.128+1)
rgl.open()
rgl.bg(color='white')

colori <- match(id_subsec, sort(unique(id_subsec)))
points3d(subsec[,1], subsec[,2], subsec[,3], col=colori)

# Results evaluation
tpr_fpr(cl.128, id_subsec)
adj.rand.index(cl.128, id_subsec)
med(cl.128, id_subsec)

# With 0 e 1
grstim.128 <- ifelse(cl.128 > 0, 1, 0)
gr.subsec <- ifelse(colori > 1, 1, 0)
adj.rand.index(grstim.128, gr.subsec)
med(grstim.128, gr.subsec)
table(grstim.128, gr.subsec)

# Bandwidth
hns(subsec) # h_ns

k <- Directional::vmf.mle(subsec)$kappa
n <- nrow(subsec)
h <- ( 8 * sinh(k)^2 /
         ( k * n * ( (1 + 4 * k^2) * sinh(2 * k) -  2 * k * cosh(2 * k) ) ) )^(1/6)
h # garcia portugues

hns(subsec)/2
# h_ns in this case is twice the smoothing parameter found using the rule 
# of thumb proposed by Garcia-Portugues

# Clustering 2
test.128.v2 <- pdfCluster.healpix(subsec, cmbdf.win5, j2, hns(subsec)/2)
plot(test.128.v2@tree, center=T)
cl.128.v2 <- test.128.v2@obs.cluster
cl.128.v2[is.na(cl.128.v2)] = 0
table(cl.128.v2, id_subsec)

tpr_fpr(cl.128.v2, id_subsec)
adj.rand.index(cl.128.v2, id_subsec)
med(cl.128.v2, id_subsec)
grstim.128.v2 <- ifelse(cl.128.v2 > 0, 1, 0)
med(grstim.128.v2, gr.subsec)
table(grstim.128.v2, gr.subsec)


# SPHERE SUBDIVISION ###########################################################
# First, I check if a tessellation with nside 64 is enough to identify 
# connected components of reasonable size
cmbdf5 <- CMBDataFrame(nside = 64, ordering = "nested", coords = 'cartesian')
cmbdf5$pix <- 1:nrow(cmbdf5)

cmbdf6 <- CMBDataFrame(nside = 128, ordering = "nested", coords = 'cartesian')
cmbdf6$pix <- 1:nrow(cmbdf6)

loc.64 <- localize(all, 6)
loc.128 <- localize(all, 7)

# With nside 64, disconnected regions cannot be found for grouping
intensity <- rep(0, 12*128^2)
tab <- table(loc.128)
intensity[as.numeric(names(tab))] <- tab
intensity <- intensity[cmbdf6$pix]
cmbdf6$I <- intensity

V(g.128)$intensity <- intensity[as.numeric(names(V(g.128)))]
empty_cells <- delete_vertices(g.128, V(g.128)$intensity != 0)
count_components(empty_cells)
comp <- components(empty_cells)

to_del_1 <- names(comp$membership[comp$membership==which.max(comp$csize)])
to_del_2 <- names(comp$membership[
  comp$membership==which(comp$csize == sort(comp$csize, decreasing=T)[2])])
occupied <- delete_vertices(g.128, c(to_del_1, to_del_2))
count_components(occupied)
components(occupied)
table(components(occupied)$csize)

cmbdf.win6 <- cmbdf6 %>% filter(pix %in% V(occupied)$name)
points3d(cmbdf.win6$x, cmbdf.win6$y, cmbdf.win6$z, col = 'purple', cex = 2)


# CLUSTERING + CORRECTION ######################################################
pos_s2 <- which(sorg_gal[,1] <= 180 & sorg_gal[,1] >= 120 & 
             sorg_gal[,2] <= 90 & sorg_gal[,2] >= 45)
subsec_s2 <-  sorg %>% cart2gal %>% as_tibble %>% 
  filter(l <= 180 &
           l >= 120 &
           b <= 90 &
           b >= 45) %>% 
  as.data.frame %>% gal2cart
id_subsec_s2 <- id_sorg[pos_s2]

pos_b2 <- which(bkg_gal[,1] <= 180 & bkg_gal[,1] >= 120 & 
               bkg_gal[,2] <= 90 & bkg_gal[,2] >= 45)
subsec_b2 <-  bkg %>% cart2gal %>% as_tibble %>% 
  filter(l <= 180 &
           l >= 120 &
           b <= 90 &
           b >= 45) %>% 
  as.data.frame %>% gal2cart
id_subsec_b2 <- rep(0, nrow(subsec_b2))

subsec2 <- rbind(subsec_s2, subsec_b2)
id_subsec2 <- c(id_subsec_s2, id_subsec_b2)
points3d(subsec_s2[,1], subsec_s2[,2], subsec_s2[,3], col=2)
points3d(subsec_b2[,1], subsec_b2[,2], subsec_b2[,3])

# Grid preparation
# cmbdf7 <- CMBDataFrame(nside = 64, ordering = "nested", coords = 'cartesian')
# cmbdf7$pix <- 1:nrow(cmbdf7)
# g.64.2 <- build.graph.v2(subsec2, cmbdf7, 6)
# 
# loc2 <- localize(subsec2, 6)
# intensity2 <- rep(0, 12*64^2)
# tab2 <- table(loc2)
# intensity2[as.numeric(names(tab2))] <- tab2
# intensity2 <- intensity2[cmbdf7$pix]
# cmbdf7$I <- intensity2
# 
# V(g.64.2)$intensity <- intensity2[as.numeric(names(V(g.64.2)))]
# empty_cells2 <- delete_vertices(g.64.2, V(g.64.2)$intensity != 0)
# count_components(empty_cells2)
# comp2 <- components(empty_cells2)
# to_del2 <- names(comp2$membership[comp2$membership==1])
# occupied2 <- delete_vertices(g.64.2, to_del2)
# 
# kids2 <- children(as.numeric(V(occupied2)$name))
# cmbdf.win7 <- cmbdf4 %>% filter(pix %in% kids2)

cmbdf4.gal <- cart2gal(as.matrix(cmbdf4[,1:3]))
cmbdf.win7 <- cmbdf4 %>% filter(cmbdf4.gal[,1] <= 182 & 
                                  cmbdf4.gal[,1] >= 119 &
                                  cmbdf4.gal[,2] <= 91 &
                                  cmbdf4.gal[,2] >= 44)

points3d(subsec_s2[,1], subsec_s2[,2], subsec_s2[,3], col=2)
points3d(subsec_b2[,1], subsec_b2[,2], subsec_b2[,3])
points3d(cmbdf.win7$x, cmbdf.win7$y, cmbdf.win7$z, col = 'purple', cex = 2)

# Clustering
test.128.v3 <- pdfCluster.healpix(subsec2, cmbdf.win7, j2, 0.02, n.grid = 3000)
plot(test.128.v3@tree, center=T)
cl.128.v3 <- test.128.v3@obs.cluster
cl.128.v3[is.na(cl.128.v3)] = 0
table(cl.128.v3, id_subsec2)

# Graphical evaluation
points3d(subsec2[,1], subsec2[,2], subsec2[,3], col=cl.128.v3+1)
rgl.open()
rgl.bg(color='white')

colori.v3 <- match(id_subsec2, sort(unique(id_subsec2)))
points3d(subsec2[,1], subsec2[,2], subsec2[,3], col=colori.v3)

# TPR, FPR, ARI, MED
tpr_fpr(cl.128.v3, id_subsec2)
adj.rand.index(cl.128.v3, id_subsec2)
med(cl.128.v3, id_subsec2)

# With 0 e 1
grstim.128.v3 <- ifelse(cl.128.v3 > 0, 1, 0)
gr.subsec.v3 <- ifelse(colori.v3 > 1, 1, 0)
adj.rand.index(grstim.128.v3, gr.subsec.v3)
med(grstim.128.v3, gr.subsec.v3)
table(grstim.128.v3, gr.subsec.v3)

# ..............................................................................
# I try to implement the correction procedure.
# First, I need to divide the dataset into the first 30% and the second 70%.
bkg30 <- bckF_sim %>% arrange(TIME)
bkg70 <- bkg30[(round(nrow(bckF_sim)*0.3)+1):nrow(bckF_sim),]
bkg30 <- bkg30[1:round(nrow(bckF_sim)*0.3),]
bkg70 <- bkg70 %>% select(L, B) %>% gal2cart()
bkg30 <- bkg30 %>% select(L, B) %>% gal2cart()

all70 <- rbind(sorg, bkg70)

# The sources remain the same
pos_s70 <- pos_s2
subsec_s70 <- subsec_s2
id_subsec_s70 <- id_subsec_s2
# The background observations change
# For the second 70%
bkg70_gal <- cart2gal(bkg70)
pos_b70 <- which(bkg70_gal[,1] <= 180 & bkg70_gal[,1] >= 120 & 
                   bkg70_gal[,2] <= 90 & bkg70_gal[,2] >= 45)
subsec_b70 <-  bkg70 %>% cart2gal %>% as_tibble %>% 
  filter(l <= 180 &
           l >= 120 &
           b <= 90 &
           b >= 45) %>% 
  as.data.frame %>% gal2cart
id_subsec_b70 <- rep(0, nrow(subsec_b70))
subsec70 <- rbind(subsec_s70, subsec_b70)
id_subsec70 <- c(id_subsec_s70, id_subsec_b70)
# First 30%
subsec_b30 <-  bkg30 %>% cart2gal %>% as_tibble %>% 
  filter(l <= 180 &
           l >= 120 &
           b <= 90 &
           b >= 45) %>% 
  as.data.frame %>% gal2cart
id_subsec_b30 <- rep(0, nrow(subsec_b30))

points3d(subsec_s70[,1], subsec_s70[,2], subsec_s70[,3], col=2)
points3d(subsec_b70[,1], subsec_b70[,2], subsec_b70[,3])
points3d(subsec_b30[,1], subsec_b30[,2], subsec_b30[,3], col='purple3')

# Perform clustering in the considered region
test.128.70 <- pdfCluster.healpix(subsec70, cmbdf.win7, j2, 0.0175, n.grid = 3000)
cl.128.70 <- test.128.70@obs.cluster
cl.128.70[is.na(cl.128.70)] = 0
table(cl.128.70, id_subsec70)
tpr_fpr(cl.128.70, id_subsec70)
adj.rand.index(cl.128.70, id_subsec70)
med(cl.128.70, id_subsec70)

points3d(subsec70[,1], subsec70[,2], subsec70[,3], col=cl.128.70+1)
rgl.open()
rgl.bg(color='white')

colori.70 <- match(id_subsec70, sort(unique(id_subsec70)))
points3d(subsec70[,1], subsec70[,2], subsec70[,3], col=colori.70)

grstim.128.70 <- ifelse(cl.128.70 > 0, 1, 0)
gr.subsec.70 <- ifelse(colori.70 > 1, 1, 0)
adj.rand.index(grstim.128.70, gr.subsec.70)
med(grstim.128.70, gr.subsec.70)
table(grstim.128.70, gr.subsec.70)

# Perform the correction
grstim.128.70.copy <- grstim.128.70
cl.128.70.copy <- cl.128.70
obs <- localize(subsec70, j2)
dens.bkg <- compute_dens(subsec_b30, cmbdf.win7, j2, 0.0175)
dens_diff <- test.128.70@pdf - dens.bkg
dens_diff <- sapply(1:nrow(subsec70), 
                    function(i) dens_diff[which(cmbdf.win7$pix == obs[i])])
grstim.128.70.copy[dens_diff <= 0] <- 0
adj.rand.index(grstim.128.70.copy, gr.subsec.70)
med(grstim.128.70.copy, gr.subsec.70)
table(grstim.128.70.copy, gr.subsec.70)
cl.128.70.copy[dens_diff <= 0] <- 0
tpr_fpr(cl.128.70.copy, id_subsec70)

