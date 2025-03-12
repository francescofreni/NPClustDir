cart2gal <- function(xcart)
{
  r2d=180/pi
  l<-atan2(xcart[,2],xcart[,1])*r2d
  l[l<0.]<-l[l<0.]+360 # if the value are negative sum 2pi
  b <-atan2(xcart[,3], sqrt(xcart[,1]^2+xcart[,2]^2))*r2d
  return (cbind(l,b))
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


pol2cart <- function(xpol) {
  theta <- xpol[,1]
  phi <- xpol[,2]
  x <- cos(theta)
  y <- sin(theta) * cos(phi)
  z <- sin(theta) * sin(phi)
  return(cbind(x,y,z))
}


# @param data data set di punti che si vogliono localizzare
# @param j indica la dimensione della griglia (nside = 2^j)
# @return vettore di lunghezza pari al numero di punti e con 
localize <- function(data, j) {
  loc <- sapply(1:nrow(data), function(i) {
    cell <- nestSearch(target = data[i,], nside = 2^(j+1))
    return(rcosmo::parent(cell$pix))
  })
  return(loc)
}


# Dato l'indice di un pixel e j (dimensione griglia), trova i pixel adiacenti
# av.pix (available pixels) indica i pixel contenuti nella griglia, perché 
# si potrebbe avere una sottogriglia piuttosto che una griglia
nb <- function(p, j, av.pix) {
  neig <- as.integer(neighbours(p, j))
  if ( length(neig) == 9 ) {
    poss.neig <- neig[c(2,4,6,8)] # possibili vicini
    poss.neig <- poss.neig[poss.neig %in% av.pix]
    return(poss.neig)
  } else {
    nside <- 2^j
    sequence <- c(0, seq(nside^2, 12*nside^2, nside^2))
    pos <- min(which(sequence >= p))
    
    if ( p < (sequence[pos-1] + sequence[pos])/2 ) {
      poss.neig <- neig[c(2,3,5,7)]
      poss.neig <- poss.neig[poss.neig %in% av.pix]
      return(poss.neig)
    }
    poss.neig <- neig[c(2,4,6,7)]
    poss.neig <- poss.neig[poss.neig %in% av.pix]
    return(poss.neig)
  }
}


# Data una finestra in input (volendo potrebbe essere anche l'intera sfera) che
# contiene anche l'informazione relativa all'indice dei pixel, le osservazioni
# e la dimensione della griglia, ritorna: i centri della griglia considera;
# la lista contentente, per ogni cella, i suoi vicini; per ogni osservazione 
# la cella (indice del pixel) in cui essa è contenuta; un vettore contenente,
# per ogni cella, il numero di osservazioni che cadono in quella cella
comp.centr.wei <- function(data, win, j) {
  loc <- localize(data, j)
  ns <- 2^j
  intensity <- rep(0, 12*ns^2)
  tab <- table(loc)
  intensity[as.numeric(names(tab))] <- tab
  intensity <- intensity[win$pix]
  win$I <- intensity
  
  l <- lapply(1:nrow(win), function(i) nb(win$pix[i], j, win$pix))
  adj <- lapply(1:length(l), function(i) {
    sapply(1:length(l[[i]]), 
           function(x) which(win$pix == l[[i]][x]))
  })
  
  return(list(tt = win[,1:3], adj = adj, obs = loc, n = win$I))
}


# stima kernel binned con parametro di lisciamento variabile
vmf.bskde.2 = function (x, h, n, y=x){
  N <- sum(n)
  d = tcrossprod(y, x)*(1/h^2)
  argsum = (exp(d - 2*log(h) - log(2*pi) - 1/h^2)) %*% n
  f = exp(log(argsum))/N
  f
}


# Modifico opportunamente la funzione pdfCluster
setClass("pdfCluster.healpix", 
         representation(call="call", x="matrix", pdf="vector", 
                        nc="list",  graph= "list", cluster.cores="ANY", 
                        tree="dendrogram", noc="numeric", obs.cluster='vector'))

pdfCluster.healpix <- function(x, win, j, h, 
                               n.grid=min(round((5 + sqrt(nrow(win)))*4), nrow(win)),
                               ...)
{
  call <- match.call()
  x <- data.matrix(x)
  if(any(sapply(x[1,],is.factor))) stop("All variables must be numeric")
  if(any(apply(x,2, pdfCluster:::check.discrete)<10)) 
    warning("One or more variables look to be discrete: pdfCluster is designed for continuous data")
  if(any(is.na(x))) 
  {warning(cat("NA in object",deparse(substitute(x)), "have been omitted", "\n"));  x <- na.omit(x)}
  
  # compute the elements to be used for both density estimation and connected components
  l <- comp.centr.wei(x, win, j)
  tt <- as.matrix(l$tt)
  graph.nb <- l$adj # this is my graph
  obs <- l$obs
  n <- l$n
  
  # density estimation
  coarse_estimate <- do.call(vmf.bskde.2, list(x = tt, h = h, n = n))
  print(range(coarse_estimate))
  #h_variable <- h/sqrt(coarse_estimate)
  g <- exp(mean(log(coarse_estimate)))
  h_variable <- h * (1/g * coarse_estimate)^(-0.5)
  print(range(h_variable))
  estimate <- do.call(vmf.bskde.2, list(x = tt, h = as.vector(h_variable), n = n))
  print(range(estimate))
  
  # check given arguments
  N <- nrow(tt)
  if (n.grid > N) {
    warning("n.grid too large, set equal to N")
    n.grid <- min(n.grid, nrow(tt))
  }
  
  # connected components
  nc <- pdfCluster:::num.con(tt, estimate, graph.nb,
                             profile.grid = n.grid-2, correct=TRUE)  
  struct <- pdfCluster:::con2tree(nc, estimate)
  if (struct$bad) {
    message("No output given")
  }
  else {
    g <- struct$g
    g[struct$g == 0] <- NA
    pdf.estim <- estimate
    names(l$adj) <- win$pix
    obs.cl <- sapply(1:dim(x)[1], function(i) g[which(as.numeric(names(l$adj)) == obs[i])])
    out <- new("pdfCluster.healpix", call = call, x = data.matrix(x), 
               pdf = pdf.estim,  nc = nc, graph = graph.nb, cluster.cores = g,
               tree = struct$tree, noc = struct$noc, obs.cluster = obs.cl)
    out
  }
}


# Calcola la densita' per poter effettuare la correzione
compute_dens <- function(x, win, j, h) {
  l <- comp.centr.wei(x, win, j)
  tt <- as.matrix(l$tt)
  obs <- l$obs 
  n <- l$n
  
  coarse_estimate <- do.call(vmf.bskde.2, list(x = tt, h = h, n = n))
  g <- exp(mean(log(coarse_estimate)))
  h_variable <- h * (1/g * coarse_estimate)^(-0.5)
  estimate <- do.call(vmf.bskde.2, list(x = tt, h = as.vector(h_variable), n = n))
  
  # Calcolo la densità di ogni osservazione (assegno la densita' della cella a
  # cui appartiene).
  #est_obs <- sapply(1:dim(x)[1], function(i) estimate[which(nomi == obs[i])]) 
  
  return(estimate)
}

# This method constructs the graph, but with a larger nside, it might be inefficient.
# However, this operation is done once and for all.
# The function takes the window and nside as input.
# The window must contain coordinate information, intensity (number of observations per cell), 
# and pixel indices.
build.graph <- function(data, win, j) {
  # Construct the edges of the graph
  # To improve efficiency, we preallocate a data frame of the correct size.
  # By the handshaking lemma, the sum of all vertex degrees equals twice the number of edges.
  # Since each vertex connects to four nodes, we estimate the matrix size accordingly.
  dim <- (nrow(win) * 4)/2
  mat <- matrix(NA, ncol=2, nrow=dim)
  
  count = 1
  for (i in 1:nrow(win)) {
    names <- nb(win$pix[i], j, win$pix)
    idxs <- unlist(sapply(1:length(names), function(x) which(win$pix == names[x])))
    idxs <- idxs[!(idxs %in% mat[,1])] # to avoid duplicates
    if (length(idxs) > 0) {
      mat[count:(count+length(idxs)-1),1] <- rep(i, length(idxs))
      mat[count:(count+length(idxs)-1),2] <- idxs
      count <- count + length(idxs)
    }
  }
  
  mat <- mat[!rowSums(is.na(mat)),]
  g <- graph_from_data_frame(data.frame(mat), directed = F)
  
  if (all(is.na(win$I))) {
    loc <- localize(data, j)
    intensity <- rep(0, 12*ns^2)
    tab <- table(loc)
    intensity[as.numeric(names(tab))] <- tab
    intensity <- intensity[win$pix]
    win$I <- intensity
  }
  
  pos <- as.numeric(names(V(g)))
  V(g)$name <- win$pix[pos]
  V(g)$intensity <- intensity[pos]
  return(g)
}

# I separate the graph construction from the calculation of the  
# number of observations per cell.
build.graph.v2 <- function(data, win, j) {
  dim <- (nrow(win) * 4)/2
  mat <- matrix(NA, ncol=2, nrow=dim)
  
  count = 1
  for (i in 1:nrow(win)) {
    names <- nb(win$pix[i], j, win$pix)
    idxs <- unlist(sapply(1:length(names), function(x) which(win$pix == names[x])))
    idxs <- idxs[!(idxs %in% mat[,1])] # to avoid duplicates
    if (length(idxs) > 0) {
      mat[count:(count+length(idxs)-1),1] <- rep(i, length(idxs))
      mat[count:(count+length(idxs)-1),2] <- idxs
      count <- count + length(idxs)
    }
  }
  
  mat <- mat[!rowSums(is.na(mat)),]
  g <- graph_from_data_frame(data.frame(mat), directed = F)
  
  pos <- as.numeric(names(V(g)))
  V(g)$name <- win$pix[pos]
  return(g)
}

tpr_fpr <- function(labels1, labels2){
  n <- length(labels1)
  N <- table(labels1, labels2)
  N <- N[-1,-1]
  r <- nrow(N)
  s <- ncol(N)
  if (r>s) {
    N <- t(N)
    r <- nrow(N)
    s <- ncol(N)
    optimal.permutation <- solve_LSAP(N, maximum = T)
    mat <- matrix(0, nrow = r, ncol = s)
    mat[cbind(seq_along(optimal.permutation), optimal.permutation)] <- 
      N[cbind(seq_along(optimal.permutation), optimal.permutation)]
    tpr <- sum(rowSums(mat)!=0)/nrow(mat)
    fpr <- sum(colSums(mat)==0)/ncol(mat)
  } else {
    N <- rbind(N, matrix(0, nrow = s-r, ncol = s))
    optimal.permutation <- solve_LSAP(N, maximum = T)
    mat <- matrix(0, nrow = s, ncol = s)
    mat[cbind(seq_along(optimal.permutation), optimal.permutation)] <- 
      N[cbind(seq_along(optimal.permutation), optimal.permutation)]
    if (s>r) {mat <- mat[-((r+1):s),]}
    tpr <- sum(colSums(mat)!=0)/ncol(mat)
    fpr <- sum(rowSums(mat)==0)/nrow(mat)
  }
  return(list(TPR = tpr, FPR = fpr))
}


med <- function(labels1, labels2){
  n <- length(labels1)
  N <- table(labels1, labels2)
  r <- nrow(N)
  s <- ncol(N)
  if (r>s) N <- t(N); r <- nrow(N); s <- ncol(N)
  if (r<s) N <- rbind(N, matrix(0, nrow = s-r, ncol = s))
  M <- matrix(rowSums(N), nrow = s, ncol = s) +
    matrix(colSums(N), nrow = s, ncol = s, byrow = TRUE) - 2 * N
  optimal.permutation <- solve_LSAP(M)
  result <- sum(M[cbind(seq_along(optimal.permutation),
                        optimal.permutation)]) / (2 * n)
  return(result)
}