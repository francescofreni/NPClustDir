# LOADING PACKAGES #############################################################
library(pdfCluster)
library(mvtnorm)
library(rgl)
library(icosa)
library(plotly)
library(igraph)
library(mvmesh)
source('utils.R')
source('functions.R')


# DATA GENERATION ##############################################################
xbg3 <- cbind(runif(3000,22.5,50), c(runif(300,-2,2), 
                                     runif(300,-4,4), 
                                     runif(300,-6,6), 
                                     runif(300,-8,8), 
                                     runif(300,-10,10), 
                                     runif(300,-12,12), 
                                     runif(300,-14,14),
                                     runif(300,-16,16),
                                     runif(300,-18,18),
                                     runif(300,-20,20)))
plot(xbg3, pch=20)
xsg3 <- rbind(rmvnorm(400, c(45,-2), diag(1,2)),
              rmvnorm(400, c(27.5,5), diag(.5,2)),
              rmvnorm(400, c(31.8,3), diag(1.1,2)))
xg3 <- rbind(xbg3, xsg3) # Galactic coordinates
xc3 <- gal2cart(xg3) # Cartesian coordinates
lab3 <- rep(c(0,1,2,3), c(3000,400,400,400))
plot(xg3, pch=20, col=lab3+1)

# Clusters closer together
xbg4 <- cbind(runif(3000,22.5,50), c(runif(300,-2,2), 
                                     runif(300,-4,4), 
                                     runif(300,-6,6), 
                                     runif(300,-8,8), 
                                     runif(300,-10,10), 
                                     runif(300,-12,12), 
                                     runif(300,-14,14),
                                     runif(300,-16,16),
                                     runif(300,-18,18),
                                     runif(300,-20,20)))
plot(xbg4, pch=20)
xsg4 <- rbind(rmvnorm(400, c(45,-2), diag(1,2)),
              rmvnorm(400, c(27.5,5), diag(.5,2)),
              rmvnorm(400, c(30.5,3), diag(c(1,2))))
xg4 <- rbind(xbg4, xsg4) # Galactic coordinates
xc4 <- gal2cart(xg4) # Cartesian coordinates
lab4 <- rep(c(0,1,2,3), c(3000,400,400,400))
plot(xg4, pch=20, col=lab4+1)


# SIMULATIONS ##################################################################
# Fixed bandwidth, tessellation = c(3,20)
grid_3_20 <- trigrid(tessellation=c(3, 20), radius=1)
subgrid_3_20 <- grid_3_20[c(lomin=20, lomax=52.5, lamin=-22.5,lamax=22.5)]
plot3d(subgrid_3_20, guides=F, arcs=T, sphere=1)
points3d(xc3[,1],xc3[,2],xc3[,3],col=lab3+1,size=2)
test.ico_3_20 <- pdfCluster.data.2(xc3, subgrid_3_20, 0.03)
plot(test.ico_3_20@tree, center=T)
cl_3_20 <- test.ico_3_20@obs.cluster
cl_3_20[is.na(cl_3_20)] = 0
table(cl_3_20, lab3)

# Variable bandwidth, tessellation = c(5,20)
grid_5_20 <- trigrid(tessellation=c(5, 20), radius=1)
subgrid_5_20 <- grid_5_20[c(lomin=20, lomax=52.5, lamin=-22.5,lamax=22.5)]
plot3d(subgrid_5_20, guides=F, arcs=T, sphere=1)
points3d(xc3[,1],xc3[,2],xc3[,3],col=lab3+1,size=2)
test.ico_5_20.new <- pdfCluster.fermi.2(xc3, subgrid_5_20, 0.08)
plot(test.ico_5_20.new@tree, center=T)
cl_5_20.new <- test.ico_5_20.new@obs.cluster
cl_5_20.new[is.na(cl_5_20.new)] = 0
table(cl_5_20.new, lab3)
