ang.dist.cart <- function(x, cart =T){
  r2d=180/pi
  d2r=pi/180
  n <- nrow(x)
  if(cart) xgal <- cart2gal(xcart)*d2r else xgal <- xgal*d2r
  out <- matrix(0, n,n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      l1 <- xgal[i,1]
      l2 <- xgal[j,1]
      b1 <- xgal[i,2]
      b2 <- xgal[j,2]
      delta_l <- abs(l1-l2)
      radq <- sqrt((cos(b2)*sin(delta_l))^2+(cos(b1)*sin(b2)-sin(b1)*cos(b2)*cos(delta_l))^2)
      d <- sin(b1)*sin(b2)+cos(b1)*cos(b2)*cos(delta_l)
      tg <- atan2(radq,d)
      tg[tg<0.] <- tg+2*pi
      out[i,j] <- out[j,i]  <- tg*r2d
    }
  }
  out
}


cart2gal <- function(xcart)
{
  r2d=180/pi
  l<-atan2(xcart[,2],xcart[,1])*r2d
  l[l<0.]<-l[l<0.]+360 # if the value are negative sum 2pi
  b <-atan2(xcart[,3], sqrt(xcart[,1]^2+xcart[,2]^2))*r2d
  return (cbind(l,b))
}


cartesian2gal <- function(x,y,z)
{
  r2d=180/pi
  d2r=pi/180
  l<-atan2(y,x)*r2d
  l[l<0.]<-l[l<0.]+360 # if the value are negative sum 2pi
  b<-atan2(z, sqrt(x*x+y*y))*r2d
  return (cbind(l,b))
}


cutmesh <- function(mesh, La, Lb, Ba, Bb){
  data <- cart2gal(mesh$V)
  indsel <- which(data[,1]>=La & data[,1]<Lb & data[,2]>Ba & data[,2]<=Bb)
  allgal <- data[indsel,]
  allcart <- gal2cart(allgal)
  seltriangoli <- which(mesh$SVI[1,]%in%indsel)
  seltriangoli <- c(seltriangoli, which(mesh$SVI[2,]%in%indsel))
  seltriangoli <- c(seltriangoli, which(mesh$SVI[3,]%in%indsel))
  seltriangoli <- unique(seltriangoli)
  indselwithb <- mesh$SVI[,seltriangoli]
  allgalwithb <- data[unique(as.vector(indselwithb)),]
  allcartwithb <- gal2cart(allgalwithb)
  S <- mesh$S[,,seltriangoli]
  SVI <- mesh$SVI[,seltriangoli]
  list(meshgal = allgal, meshcart = allcart, meshgalwithb=allgalwithb, meshcartwithb =allcartwithb, S=S, SVI =SVI)
}


plot.cutmesh<- function(cutmesh, gal=FALSE, ...){
  if(!gal){
    for (k in 1:dim(cutmesh$S)[3]){ 
      for (i in 1:(nrow(cutmesh$S[,,k])-1)) {
        for (j in (i+1):(nrow(cutmesh$S[,,k])) ) {
          lines3d(cutmesh$S[,,k][c(i, j),1],cutmesh$S[,,k][c(i, j),2],cutmesh$S[,,k][c(i, j),3], ...) 
        }  
      }
      rgl.triangles(cutmesh$S[,,k]) # needed?
    }# else
    #{
    #  
  }
}


estrai.real <- function(data= alldata, sorgenti= sorgenti,La, Lb, Ba, Bb){
  all <- data[data[,4]>=La & data[,4]<Lb & data[,5]>Ba & data[,5]<=Bb,]
  cooso<-sorgenti[sorgenti[,2]>=La & sorgenti[,2]<Lb & sorgenti[,3]>Ba & sorgenti[,3]<=Bb,2:3]
  rownames(cooso)<- 1:nrow(cooso)
  list(all=all, cooso=cooso)
}


estrai0 <- function(data= alldata, sorgenti=sorgenti, estremi){
  which.sel <- match(data$ID, factor(sorgenti$name))
  data <- data[!is.na(which.sel),]
  lab <- sorgenti$id[which.sel[!is.na(which.sel)]]
  data <- data.frame(id=lab, data)
  La <- estremi[1]
  Lb <- estremi[2]
  Ba <- estremi[3]
  Bb <- estremi[4] 
  datatemp <- data[,5:6]
  if(La<0|Lb<0){
    ind.neg <- which(data[,5]> 180)
    datatemp[ind.neg, 1] <- data[ind.neg, 5] -360    
  }
  if(Ba>90|Bb>90){
    ind.pos <- which(data[,6]< 0)
    datatemp[ind.pos, 2] <- 180 + data[ind.pos, 6]     
  }
  all <- data[datatemp[,1]>=La & datatemp[,1]<=Lb & datatemp[,2]>=Ba & datatemp[,2]<=Bb,]
  which.sel <- match(all$id, sorgenti$id)
  cooso <-sorgenti[unique(which.sel),3:4]
  rownames(cooso)<- 1:nrow(cooso)
  list(all=all, cooso=cooso)
}


gal2cart <- function(xgal){
  xrad <- xgal*pi/180
  l <- xrad[,1]
  b <- xrad[,2]
  x<-cos(l)*cos(b)
  y<-sin(l)*cos(b)
  z<-sin(b)
  cbind(x,y,z)
}


gal2cartesian <- function(file= NULL, L=NULL, B=NULL){
  r2d=180/pi
  d2r=pi/180
  
  if(!is.null(file)) {
    l<-file$L*r2d
    b<-file$B*r2d} else {
      l <- L*r2d
      b <- B*r2d}
  x<-cos(l)*cos(b)
  y<-sin(l)*cos(b)
  z<-sin(b)
  return (list(coord=cbind(x,y,z)))
  #return (cbind(x,y,z))
  #return(file$x, file$y, file$z)
}


gal2cartesian2 <- function(file= NULL, L=NULL, B=NULL){
  r2d=180/pi
  d2r=pi/180
  
  if(!is.null(file)) {
    l<-file$L*d2r
    b<-file$B*d2r} else {
      l <- L*d2r
      b <- B*d2r}
  x<-cos(l)*cos(b)
  y<-sin(l)*cos(b)
  z<-sin(b)
  return (list(coord=cbind(x,y,z), ft_ph=file$PH, ft_id=file$ID))
  #return(file$x, file$y, file$z)
  
}


regioni <- function(dati, mesh) {
  occ <- occupied(mesh, dati)@values
  m <- mesh[occ]
  com <- components(m@graph)
  lista_reg <- vector(length=count_components(m@graph), mode='list')
  lista_obs <- vector(length=count_components(m@graph), mode='list')
  limit <- unique(com$membership)
  for (c in limit) {
    print(c)
    lista_reg[c] <- m[row.names(m@faces) %in% names(com$membership[com$membership==c])]
    lista_obs[[c]] <- dati[which(locate(lista_reg[[c]], dati) != 'F0'),]
  }
  return(list(regioni=lista_reg, dati=lista_obs))
}


risultati <- function() {
  # vector that will contain the cluster membership IDs. 0 represents the background
  clus_obs <- rep(0, dim(all)[1])
  
  # first, classify the observations in the poles
  com_poli <- components(poli@graph)
  ret_poli <- vector(length=count_components(poli@graph), mode='list')
  limit <- unique(com_poli$membership)
  for (c in limit) {
    #print(c)
    nomi <- names(com_poli$membership[com_poli$membership==c])
    if (length(nomi) > 2) {
      pos <- which(posizioni9 %in% nomi)
      dati <- all[pos,]
      ret <- pdfCluster.data.2.new.2(x=dati, mesh=poli, h=h, nomi=nomi)
      ret_poli[[c]] <- ret
      obs_ret <- ret@obs.cluster
      obs_ret[is.na(obs_ret)] = 0
      obs_ret[obs_ret != 0] <- obs_ret[obs_ret != 0] + max(clus_obs) + 1
      clus_obs[pos] <- obs_ret
    }
  }
  
  # now classify the observations in the central portion
  com_central <- components(final_central_2)
  ret_central <- vector(length=count_components(final_central_2@graph), mode='list')
  limit2 <- unique(com_central$membership)
  for (c in limit2) {
    print(c)
    nomi <- names(com_central$membership[com_central$membership==c])
    if (length(nomi) > 2) {
      pos <- which(posizioni11 %in% nomi)
      dati <- all[pos,]
      ret <- pdfCluster.data.2.new.2(x=dati, mesh=final_central_2, h=h, nomi=nomi)
      ret_central[[c]] <- ret
      obs_ret <- ret@obs.cluster
      obs_ret[is.na(obs_ret)] = 0
      obs_ret[obs_ret != 0] <- obs_ret[obs_ret != 0] + max(clus_obs) + 1
      clus_obs[pos] <- obs_ret
    }
  }
  
  return(list(cluster=clus_obs, ret_poli=ret_poli, ret_central=ret_central))
}


subregions <- function(mesh, points) {
  cells <- locate(mesh, points)
  tCell <- table(cells)
  
  l <- as.list(rep(0, dim(mesh@faces)[1]))
  names(l) <- row.names(mesh@faces)
  l[row.names(tCell)] <- tCell
  
  soglia <- round(mean(tCell)) + 1
  reg <- list() # each element will contain the names of the faces belonging to
                # the subregion
  visit <- c() # visited and discarded cells
  k = 1 # subregion number
  
  end = FALSE
  
  while (length(l) > 0) {
    #print(length(l))
    print(max(unlist(l)))
    if (max(unlist(l)) < soglia) {
      end = TRUE
      break
    }
    
    maxcell <- names(l)[which.max(unlist(l))]
    i = 1
    lastbelt <- vicinity(mesh, maxcell, order=i, self=F)
    
    # it can happen that a neighbor has already been visited. In that case, 
    # I start from a single cell
    if ( any(sapply(1:length(lastbelt), function(j) lastbelt[j] %in% visit)) )
      lastbelt <- maxcell
    
    # this is the inner region that I will discard to keep only the belt
    lastreg <- lastbelt
    
    while (!all(unlist(l[lastbelt]) < soglia)) {
      i = i + 1
      lastbelt <- vicinity(mesh, maxcell, order=i, self=F)
      
      discard <- sapply(1:length(lastbelt),
                        function(j) lastbelt[j] %in% lastreg)
      lastreg <- lastbelt
      lastbelt <- lastbelt[!discard] # I only keep the belt
      
      # remove those that have already been visited
      if (length(lastbelt) > 0) {
        discard2 <- sapply(1:length(lastbelt),
                           function(j) lastbelt[j] %in% visit)
        lastbelt <- lastbelt[!discard2]
      }
    }
    
    # reconstruct the region found
    region <- vicinity(mesh, maxcell, order=i)
    discard <- sapply(1:length(region), function(j) region[j] %in% visit)
    region <- region[!discard]
    reg[[k]] <- region
    
    # add the found region to the visited cells
    visit <- c(visit, region)
    
    # update the size of the list
    pos <- sapply(1:length(l), function(j) names(l)[j] %in% reg[[k]])
    l <- l[!pos]
    k = k + 1
  }
  if (end) {
    reg[[k]] <- setdiff(row.names(mesh@faces), visit)
  }
  return(reg)
}


vmf.bskde.contour <- function(u, n, h, ngrid = 100, den.ret=T, full=F) {
  N <- dim(u)[1]
  x <- euclid(u)
  
  if (full) {
    x1 <- seq( 0, 180, length = ngrid )  
    x2 <- seq( 0, 360, length = ngrid )  
  } else {
    x1 <- seq( min(u[, 1]) - 5, max(u[, 1]) + 5, length = ngrid ) 
    x2 <- seq( min(u[, 2]) - 5, max(u[, 2]) + 5, length = ngrid ) 
  }
  
  mat <- matrix(nrow = ngrid, ncol = ngrid)
  
  for (i in 1:ngrid) {
    for (j in 1:ngrid) {
      y <- Directional::euclid( c(x1[i], x2[j]) )
      a <- as.vector( tcrossprod(x, y / h^2) )
      adj = max(a)
      argsum = exp(t(a - adj)) %*% n
      can = exp(log(argsum) + adj - 2*log(h) - log(2*pi) - 1/h^2)/ngrid
      if (abs(can) < Inf)   mat[i, j] <- can
    }
  }
  
  if (den.ret) {
    return(list(lat = x1, long = x2, h = h, den = mat))
  } else {
    par(fg = NA)
    filled.contour(x1, x2, mat,
                   nlevels = 200,
                   color.palette = colorRampPalette( c( "blue",
                                                        "cyan",
                                                        "yellow",
                                                        "red") ),
                   plot.axes = { axis(1, col = "black", cex.axis = 1.2);
                     axis(2, col = "black", cex.axis = 1.2);
                     contour(x1, x2, mat,
                             col="black",
                             nlevels=10,
                             labcex = 0.8,
                             lwd = 1.5,
                             add = TRUE) },
                   key.axes = {axis(4, col = "black", cex.axis = 1.2)},
                   xlab = "Latitude",
                   ylab = "Longitude",
                   cex.lab = 1.2)
  }
}