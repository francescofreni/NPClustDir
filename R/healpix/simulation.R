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
source('utils.R')
source('functions.R')


# DATA GENERATION  #############################################################
xbg <- cbind(runif(3000,22.5,50), c(runif(300,-2,2), 
                                    runif(300,-4,4), 
                                    runif(300,-6,6), 
                                    runif(300,-8,8), 
                                    runif(300,-10,10), 
                                    runif(300,-12,12), 
                                    runif(300,-14,14),
                                    runif(300,-16,16),
                                    runif(300,-18,18),
                                    runif(300,-20,20)))
xsg <- rbind(rmvnorm(400, c(45,-2), diag(1,2)),
             rmvnorm(400, c(27.5,5), diag(.5,2)),
             rmvnorm(400, c(31.8,3), diag(1.1,2)))
xg <- rbind(xbg, xsg) 
xc <- PolToCar(xg, radius = 1) 
lab <- rep(c(0,1,2,3), c(3000,400,400,400))

ns <- 64
j <- 6
win <- CMBWindow(theta = c(pi*0.36,pi*0.36,pi*0.62,pi*0.62),
                 phi = c(pi/3.25,pi/10,pi/10,pi/3.25))
plot(win)
points3d(xc[,1], xc[,2], xc[,3], col=lab+1)
displayPixelBoundaries(nside = ns)


# ANALYSIS #####################################################################
loc <- localize(xc, j)
intensity <- rep(0, 12*ns^2)
tab <- table(loc)
intensity[as.numeric(names(tab))] <- tab
cmbdf <- CMBDataFrame(nside = ns, I = intensity, ordering = "nested",
                      coords = 'cartesian') # trovo anche le coordinate dei centri
cmbdf$pix <- 1:nrow(cmbdf)
cmbdf.win <- window(cmbdf, new.window = win)

# Plot
plot(win)
points3d(xc[,1], xc[,2], xc[,3], col=lab+1)
points3d(cmbdf.win$x, cmbdf.win$y, cmbdf.win$z, col = 'purple', cex = 2)
#text3d(cmbdf.win$x, cmbdf.win$y, cmbdf.win$z, text = cmbdf.win$pix)
displayPixelBoundaries(nside = 64)

# First test with nside=64
cmbdf <- CMBDataFrame(nside = ns, ordering = "nested", coords = 'cartesian')
cmbdf$pix <- 1:nrow(cmbdf)
cmbdf.win <- window(cmbdf, new.window = win)
test.64 <- pdfCluster.healpix(xc, cmbdf.win, j, 0.07)
plot(test.64@tree, center=T)
cl.64 <- test.64@obs.cluster
cl.64[is.na(cl.64)] = 0
table(cl.64, lab)

plot(win)
points3d(xc[,1], xc[,2], xc[,3], col=cl.64+1)
points3d(cmbdf.win$x, cmbdf.win$y, cmbdf.win$z, col = 'purple', cex = 2)


# GRAPH ########################################################################
# The next step is constructing the graph from the HEALPix structure.
# Two nodes (cells) are connected if they are adjacent.
# To avoid creating very small regions, we use a smaller tessellation (lower nside).
# Once we identify the connected components with this structure, we move to a larger nside.
# This way, we are more precise and avoid very small regions.
ns1 <- 32
j1 <- 5

# Test cases
cmbdf1 <- CMBDataFrame(nside = 2, ordering = "nested", coords = 'cartesian')
cmbdf1$pix <- 1:nrow(cmbdf1)
g <- build.graph(xc, cmbdf1, 1)

cmbdf.win1 <- window(cmbdf1, new.window = win)
g1 <- build.graph(xc, cmbdf.win1, 1)

# When removing a vertex from an igraph object, its edges are also removed.
# So, the plan is to build the graph for the entire sphere, identify empty cells, 
# and remove them along with their edges.

delete_vertices(g, 48)  # Removes the index
delete_vertices(g, '48') # Removes the node with this name

# Test with simulated data
cmbdf2 <- CMBDataFrame(nside = ns1, ordering = "nested", coords = 'cartesian')
cmbdf2$pix <- 1:nrow(cmbdf2)
g <- build.graph(xc, cmbdf2, j1)
empty_cells <- delete_vertices(g, V(g)$intensity != 0)
count_components(empty_cells)
comp <- components(empty_cells)

plot(win)
points3d(xc[,1], xc[,2], xc[,3], col=lab+1)
points3d(cmbdf.win$x, cmbdf.win$y, cmbdf.win$z, col = 'purple', cex = 2)
displayPixelBoundaries(nside = ns1)

to_del <- names(comp$membership[comp$membership==1])
occupied <- delete_vertices(g, to_del)

# This is the window created using nside 32
cmbdf.win2 <- cmbdf2 %>% filter(pix %in% V(occupied)$name)

plot(win)
points3d(xc[,1], xc[,2], xc[,3], col=lab+1)
points3d(cmbdf.win2$x, cmbdf.win2$y, cmbdf.win2$z, col = 'purple', cex = 2)
displayPixelBoundaries(nside = ns1)

# This is the grid window with nside 64, identified from connected components at nside 32.
# If a cell at nside 32 is occupied, not all of its child cells at nside 64 must be occupied.
kids <- children(as.numeric(V(occupied)$name))
cmbdf.win3 <- cmbdf %>% filter(pix %in% kids)

plot(win)
points3d(xc[,1], xc[,2], xc[,3], col=lab+1)
points3d(cmbdf.win3$x, cmbdf.win3$y, cmbdf.win3$z, col = 'purple', cex = 2)

test.graph <- pdfCluster.healpix(xc, cmbdf.win3, j, 0.057)
plot(test.graph@tree, center=T)
cl.graph <- test.graph@obs.cluster
cl.graph[is.na(cl.graph)] = 0
table(cl.graph, lab)

points3d(xc[,1], xc[,2], xc[,3], col=cl.graph+1)
