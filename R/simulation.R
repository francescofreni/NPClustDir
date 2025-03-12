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
xbg <- cbind(runif(1000,0,45), c(runif(100,-2,2), 
                                 runif(100,-4,4), 
                                 runif(100,-6,6), 
                                 runif(100,-8,8), 
                                 runif(100,-10,10), 
                                 runif(100,-12,12), 
                                 runif(100,-14,14),
                                 runif(100,-16,16),
                                 runif(100,-18,18),
                                 runif(100,-20,20)))
#plot(xbg, pch=20)
library(mvtnorm)
xsg<- rbind(rmvnorm(100, c(25,2), diag(1,2)),
            rmvnorm(100, c(10,5), diag(.5,2)),
            rmvnorm(100, c(7.8,3), diag(1.1,2)))
xg <- rbind(xbg, xsg) # galactic coordinates
xc <- gal2cart(xg) # cartesian coordinates
lab <- rep(c(0,1,2,3), c(1000,100,100,100))
#points(xsg, col=2, pch=20)
#plot(xg, pch=20, col=lab+1)
#points3d(xc[,1],xc[,2],xc[,3],col=lab+1,size=2)

# Save the mesh
mesh6 <- readRDS('mesh6.rds')
mesh5 <- readRDS('mesh5')

# Generating sub-mesh
cm <- cutmesh(mesh5, 0, 45,-20, 20)
#plot.cutmesh(cm)
#points3d(xc[,1],xc[,2],xc[,3],col=lab+1, pch=20)
cm6 <- cutmesh(mesh6, 0, 45,-20, 20)

# Generating data in a different way
xbg2 <- cbind(runif(1000,22.5,67.5), c(runif(100,-2,2), 
                                       runif(100,-4,4), 
                                       runif(100,-6,6), 
                                       runif(100,-8,8), 
                                       runif(100,-10,10), 
                                       runif(100,-12,12), 
                                       runif(100,-14,14),
                                       runif(100,-16,16),
                                       runif(100,-18,18),
                                       runif(100,-20,20)))
xsg2 <- rbind(rmvnorm(100, c(55,2), diag(1,2)),
              rmvnorm(100, c(35,5), diag(.5,2)),
              rmvnorm(100, c(38.8,3), diag(1.1,2)))
xg2 <- rbind(xbg2, xsg2)
xc2 <- gal2cart(xg2)
lab2 <- rep(c(0,1,2,3), c(1000,100,100,100))
#plot(xg2, pch=20, col=lab2+1)
cm2 <- cutmesh(mesh5, 22.5, 67.5,-20, 20)
#plot.cutmesh(cm3)
#points3d(xc3[,1],xc3[,2],xc3[,3],col=lab+1, pch=20)


# MESH AND SIMULATION ##########################################################
# Creation of the grid and subgrid
grid <- trigrid(tessellation=c(20), radius=1)
subgrid <- grid[c(lomin=-10, lomax=50, lamin=-25,lamax=25)] # inverto lat e long
plot3d(subgrid, guides=F, arcs=T, sphere=1)
points3d(xc[,1],xc[,2],xc[,3],col=lab+1,size=2)

test.ico <- pdfCluster.fermi(xc, subgrid, 0.033)
test.ico@pdf
plot(test.ico@tree, center=T)
cl <- test.ico@obs.cluster
cl[is.na(cl)] = 0
table(cl, lab)

system.time(pdfCluster.fermi(xc, subgrid, 0.033))

# Creation of a finer mesh
grid_5_20 <- readRDS('grid_5_20.rds') # trigrid(tessellation=c(5, 20), radius=1)
subgrid.5 <- grid_5_20[c(lomin=-10, lomax=50, lamin=-25,lamax=25)]
plot3d(grid_5_20, guides=F, arcs=T, sphere=1)
plot3d(subgrid.5, guides=F, arcs=T, sphere=1)
points3d(xc[,1],xc[,2],xc[,3],col=lab+1,size=2)

# Test with finer mesh
grid_3_20 <- trigrid(tessellation=c(3, 20), radius=1)
subgrid.3 <- grid_3_20[c(lomin=-10, lomax=50, lamin=-25,lamax=25)]
plot3d(subgrid.3, guides=F, arcs=T, sphere=1)
points3d(xc[,1],xc[,2],xc[,3],col=lab+1,size=2)
test.ico.3 <- pdfCluster.fermi(xc, subgrid.3, 0.04)
plot(test.ico.3@tree, center=T)
cl.3 <- test.ico.3@obs.cluster
cl.3[is.na(cl.3)] = 0
table(cl.3, lab)

# Test with data where the three clusters are more separated
subgrid.2 <- grid_3_20[c(lomin=20, lomax=70, lamin=-25,lamax=25)]
plot3d(subgrid.2, guides=F, arcs=T, sphere=1)
points3d(xc2[,1],xc2[,2],xc2[,3],col=lab2+1,size=2)
test.ico.2 <- pdfCluster.fermi(xc2, subgrid.2, 0.025)
plot(test.ico.2@tree, center=T)
cl.2 <- test.ico.2@obs.cluster
cl.2[is.na(cl.2)] = 0
table(cl.2, lab2)

l.3.20 <- comp.centr.wei.ico(xc2, subgrid.2)
tg2.3.20 <- cart2gal(l.3.20$t)
n2.3.20 <- l.3.20$n
ret.3.20 <- vmf.bskde.contour(tg2.3.20, n2.3.20, 0.025, den.ret=T)
fig.3.20 <- plot_ly(x = ret.3.20$long, y = ret.3.20$lat, z = ret.3.20$den) %>% add_surface()
hide_colorbar(fig.3.20)


# PLOT #########################################################################
grid <- trigrid(tessellation=c(20), radius=1)
plot3d(grid, guides=F, arcs=T, sphere=1)

# Plot grid sections
plot3d(grid[c(lomin=-180, lomax=-90, lamin=-90,lamax=0)], guides=F, arcs=F,
       sphere=1, col='blue')
lines3d(grid[c(lomin=-180, lomax=-90, lamin=0,lamax=90)], guides=F, arcs=F,
        sphere=1, col='red')

lines3d(grid[c(lomin=-90, lomax=0, lamin=-90,lamax=0)], guides=F, arcs=F,
        sphere=1, col='green')
lines3d(grid[c(lomin=-90, lomax=0, lamin=0,lamax=90)], guides=F, arcs=F,
        sphere=1, col='yellow')

lines3d(grid[c(lomin=0, lomax=90, lamin=-90,lamax=0)], guides=F, arcs=F,
        sphere=1, col='orange')
lines3d(grid[c(lomin=0, lomax=90, lamin=0,lamax=90)], guides=F, arcs=F,
        sphere=1, col='purple')

lines3d(grid[c(lomin=90, lomax=180, lamin=-90,lamax=0)], guides=F, arcs=F,
        sphere=1, col='brown')
lines3d(grid[c(lomin=90, lomax=180, lamin=0,lamax=90)], guides=F, arcs=F,
        sphere=1, col='violet')

# Simulated data + source data
data <- gal2cart(alldata_new[,c('L', 'B')])
sg <- gal2cart(sorgenti[,c('l', 'b')])
data_bkg <- gal2cart(bckF_sim[,c('L', 'B')])
plot3d(grid, guides=F, arcs=T, sphere=1)
points3d(data[,1], data[,2], data[,3])
points3d(data_bkg[,1], data_bkg[,2], data_bkg[,3])
points3d(sg[,1], sg[,2], sg[,3], color=2)

# Real data
data_ft <- gal2cart(ft1[,c('L', 'B')])
plot3d(grid, guides=F, arcs=T, sphere=1)
points3d(data_ft[,1], data_ft[,2], data_ft[,3])

# Portions can be taken in this way
faces3d(grid[grid@belts == 1], col='blue')
faces3d(grid[grid@belts == 2], col='red')


# GENERATION OF DENSE DATA #####################################################
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


# SIMULATIONS WITH DENSE DATA ##################################################
# tessellation = c(3,20)
grid_3_20 <- readRDS('grid_3_20.rds') # trigrid(tessellation=c(3, 20), radius=1)
subgrid_3_20 <- grid_3_20[c(lomin=20, lomax=52.5, lamin=-22.5,lamax=22.5)]
plot3d(subgrid_3_20, guides=F, arcs=T, sphere=1)
points3d(xc3[,1],xc3[,2],xc3[,3],col=lab3+1,size=2)
test.ico_3_20 <- pdfCluster.fermi(xc3, subgrid_3_20, 0.02)
plot(test.ico_3_20@tree, center=T)
cl_3_20 <- test.ico_3_20@obs.cluster
cl_3_20[is.na(cl_3_20)] = 0
saveRDS(table(cl_3_20, lab3), 'tab_3_20_v2.rds')
readRDS('tab_3_20_v2.rds') # with bandwidth equal to 0.02
readRDS('tab_3_20.rds') # with bandwidth equal to 0.0275

# tessellation = c(5,20)
subgrid_5_20 <- grid_5_20[c(lomin=20, lomax=52.5, lamin=-22.5,lamax=22.5)]
plot3d(subgrid_5_20, guides=F, arcs=T, sphere=1)
points3d(xc3[,1],xc3[,2],xc3[,3],col=lab3+1,size=2)
test.ico_5_20 <- pdfCluster.fermi(xc3, subgrid_5_20, 0.02)
plot(test.ico_5_20@tree, center=T)
cl_5_20 <- test.ico_5_20@obs.cluster
cl_5_20[is.na(cl_5_20)] = 0
saveRDS(table(cl_5_20, lab3), 'tab_5_20_v2.rds')
readRDS('tab_5_20.rds') # with bandwidth equal to 0.03
readRDS('tab_5_20_v2.rds') # with bandwidth equal to 0.02

# tessellation = c(7,20)
grid_7_20 <- readRDS('grid_7_20.rds') # trigrid(tessellation=c(7, 20), radius=1)
subgrid_7_20 <- grid_7_20[c(lomin=20, lomax=52.5, lamin=-22.5,lamax=22.5)]
plot3d(subgrid_7_20, guides=F, arcs=T, sphere=1)
points3d(xc3[,1],xc3[,2],xc3[,3],col=lab3+1,size=2)
test.ico_7_20 <- pdfCluster.fermi(xc3, subgrid_7_20, 0.02)
plot(test.ico_7_20@tree, center=T)
cl_7_20 <- test.ico_7_20@obs.cluster
cl_7_20[is.na(cl_7_20)] = 0
saveRDS(table(cl_7_20, lab3), 'tab_7_20_v2.rds')
readRDS('tab_7_20.rds') # with bandwidth equal to 0.03
readRDS('tab_7_20_v2.rds') # with bandwidth equal to 0.02

# Plot
l_7_20 <- comp.centr.wei.ico(xc3, subgrid_7_20)
tg3_7_20 <- cart2gal(l_7_20$t)
n3_7_20 <- l_7_20$n
ret_7_20 <- vmf.bskde.contour(tg3_7_20, n3_7_20, 0.02, den.ret=T)
fig_7_20 <- plot_ly(x = ret_7_20$long, y = ret_7_20$lat, z = ret_7_20$den) %>% add_surface()
hide_colorbar(fig_7_20)

l_5_20 <- comp.centr.wei.ico(xc3, subgrid_5_20)
tg3_5_20 <- cart2gal(l_5_20$t)
n3_5_20 <- l_5_20$n
ret_5_20 <- vmf.bskde.contour(tg3_5_20, n3_5_20, 0.02, den.ret=T)
fig_5_20 <- plot_ly(x = ret_5_20$long, y = ret_5_20$lat, z = ret_5_20$den) %>% add_surface()
hide_colorbar(fig_5_20)

# Creating a finer grid
grid_9_20 <- readRDS('grid_9_20.rds') # trigrid(tessellation=c(9, 20), radius=1)
grid_11_20 <- readRDS('grid_11_20.rds') # trigrid(tessellation=c(11, 20), radius=1)

# Trying with data where two clusters are very close
plot3d(subgrid_7_20, guides=F, arcs=T, sphere=1)
points3d(xc4[,1],xc4[,2],xc4[,3],col=lab4+1,size=2)
test.ico_7_20_xc4 <- pdfCluster.fermi(xc4, subgrid_7_20, 0.0175)
plot(test.ico_7_20_xc4@tree, center=T)
cl_7_20_xc4 <- test.ico_7_20_xc4@obs.cluster
cl_7_20_xc4[is.na(cl_7_20_xc4)] = 0
saveRDS(table(cl_7_20_xc4, lab4), 'tab_7_20_xc4.rds')
readRDS('tab_7_20_xc4.rds') # with bandwidth equal to 0.0175

test.ico_7_20_n.grid1000 <- pdfCluster.fermi(xc4, subgrid_7_20, 0.01, n.grid=1000)
# No output given

l_7_20_xc4 <- comp.centr.wei.ico(xc4, subgrid_7_20)
tg4_7_20 <- cart2gal(l_7_20_xc4$t)
n4_7_20 <- l_7_20_xc4$n
ret_7_20_xc4 <- vmf.bskde.contour(tg4_7_20, n4_7_20, 0.01, den.ret=T)
fig_7_20_xc4 <- plot_ly(x = ret_7_20_xc4$long, y = ret_7_20_xc4$lat,
                        z = ret_7_20_xc4$den) %>% add_surface()
hide_colorbar(fig_7_20_xc4)
# In this case, it would not recognize the two groups because they merge,
# but if I set the smoothing parameter to 0.01 the function does not return output,
# probably due to the tree or similar operations, because there would be
# too many small peaks.

# Plot of simulated data and evaluation of optimal grid
plot3d(grid_9_20, guides=F, arcs=T, sphere=1)
points3d(data[,1], data[,2], data[,3])
points3d(data_bkg[,1], data_bkg[,2], data_bkg[,3])

# Variable mesh
faces3d(grid_5_20[grid_5_20@belts %in% 1:(300/3)], col='blue')
faces3d(grid_5_20[grid_5_20@belts %in% 201:300], col='blue')
faces3d(grid_7_20[grid_7_20@belts %in% 141:280], col='red')


# VARIABLE BANDWIDTH ###########################################################
# tessellation = c(5,20)
test.ico_5_20.new <- pdfCluster.fermi.2(xc3, subgrid_5_20, 0.15)
plot(test.ico_5_20.new@tree, center=T)
cl_5_20.new <- test.ico_5_20.new@obs.cluster
cl_5_20.new[is.na(cl_5_20.new)] = 0
saveRDS(table(cl_5_20.new, lab3), 'tab_5_20_new.rds')
readRDS('tab_5_20_new.rds') # with initial bandwidth equal to 0.05

# tessellation = c(7,20)
plot3d(subgrid_7_20, guides=F, arcs=T, sphere=1)
points3d(xc3[,1],xc3[,2],xc3[,3],col=lab3+1,size=5)
test.ico_7_20.new <- pdfCluster.fermi.2(xc3, subgrid_7_20, 0.055)
plot(test.ico_7_20.new@tree, center=T)
cl_7_20.new <- test.ico_7_20.new@obs.cluster
cl_7_20.new[is.na(cl_7_20.new)] = 0
saveRDS(table(cl_7_20.new, lab3), 'tab_7_20_new_v2.rds')
readRDS('tab_7_20_new.rds') # with initial bandwidth equal to 0.05
readRDS('tab_7_20_new_v2.rds') # with initial bandwidth equal to 0.055

# tessellation = c(7,20) + dati fitti
plot3d(subgrid_7_20, guides=F, arcs=T, sphere=1)
points3d(xc4[,1],xc4[,2],xc4[,3],col=lab4+1,size=2)
test.ico_7_20_xc4.new <- pdfCluster.fermi.2(xc4, subgrid_7_20, 0.055)
plot(test.ico_7_20_xc4.new@tree, center=T)
cl_7_20_xc4.new <- test.ico_7_20_xc4.new@obs.cluster
cl_7_20_xc4.new[is.na(cl_7_20_xc4.new)] = 0
saveRDS(table(cl_7_20_xc4.new, lab4), 'tab_7_20_xc4_new_v3.rds')
readRDS('tab_7_20_xc4_new.rds') # with initial bandwidth equal to 0.03
readRDS('tab_7_20_xc4_new_v2.rds') # with initial bandwidth equal to 0.05
readRDS('tab_7_20_xc4_new_v3.rds') # with initial bandwidth equal to 0.055
