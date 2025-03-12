# MAIN #########################################################################
vmf.bskde <- function (x, h, n, y=x){
  N <- sum(n)
  d = tcrossprod(y, x)/h^2
  adj = max(d)
  argsum = exp(d - adj) %*% n
  f = exp(log(argsum) + adj - 2*log(h) - log(2*pi) - 1/h^2)/N
  f
}


vmf.bskde.2 <- function (x, h, n, y=x){
  N <- sum(n)
  d = tcrossprod(y, x)*(1/h^2)
  argsum = (exp(d - 2*log(h) - log(2*pi) - 1/h^2)) %*% n
  f = exp(log(argsum))/N
  f
}


comp.centr.wei.ico <- function(x, mesh) {
  l <- vicinity(mesh, row.names(mesh@faces), output='list', order=1, self=F)
  adj <- lapply(1:length(l), function(i) {
    sapply(1:length(l[[i]]), 
           function(x) which(row.names(mesh@faces) == l[[i]][x]))
  })
  obs <- locate(mesh, x)
  n <- sapply(1:length(adj), function(i) length(which(obs == row.names(mesh@faces)[i])))
  return(list(t=mesh@faceCenters, adj=adj, obs=obs, n=n))
}


setClass("pdfCluster2", representation(call="call", x="matrix", pdf="vector", 
                                       nc="list",  graph= "list", cluster.cores="ANY", 
                                       tree="dendrogram", noc="numeric", obs.cluster='vector'))


# original, variable bandwidth
pdfCluster.fermi <- function(
    x, mesh, h, 
    n.grid=min(round((5 + sqrt(dim(mesh@faces)[1]))*4), dim(mesh@faces)[1]), 
    ...
    )
  # n.grid defines the length of the grid on which evaluating the connected
  # components of the density level sets.
{
  call <- match.call()
  x <- data.matrix(x)
  if(any(sapply(x[1,],is.factor))) stop("All variables must be numeric")
  #if(any(apply(x,2, pdfCluster:::check.discrete)<10)) 
  #  warning("One or more variables look to be discrete: pdfCluster is designed for continuous data")
  if(any(is.na(x))) 
  {warning(cat("NA in object",deparse(substitute(x)), "have been omitted", "\n"));  x <- na.omit(x)}
  
  # compute the elements to be used for both density estimation and connected components
  l <- comp.centr.wei.ico(x, mesh)
  t <- l$t
  graph.nb <- l$adj # this is my graph
  obs <- l$obs 
  n <- l$n
  
  # density estimation
  coarse_estimate <- do.call(vmf.bskde.2, list(x = t, h = hns(x)/10, n = n))
  print(range(coarse_estimate))
  #h_variable <- h/sqrt(coarse_estimate)
  g <- exp(mean(log(coarse_estimate)))
  h_variable <- h * (1/g * coarse_estimate)^(-0.5)
  #h_variable <- h * (coarse_estimate)^(-0.5)
  print(range(h_variable))
  estimate <- do.call(vmf.bskde.2, list(x = t, h = as.vector(h_variable), n = n))
  print(range(estimate))
  
  # check given arguments
  N <- nrow(t)
  if (n.grid > N) {
    warning("n.grid too large, set equal to N")
    n.grid <- min(n.grid, nrow(t))
  }
  
  # connected components
  nc <- pdfCluster:::num.con(t, estimate, graph.nb,
                             profile.grid = n.grid-2, correct=TRUE)  
  struct <- pdfCluster:::con2tree(nc, estimate)
  if (struct$bad) {
    message("No output given")
  }
  else {
    g <- struct$g
    g[struct$g == 0] <- NA
    pdf.estim <- estimate
    names(l$adj) <- row.names(mesh@faces)
    obs.cl <- sapply(1:dim(x)[1], function(i) g[which(names(l$adj) == obs[i])])
    out <- new("pdfCluster2", call = call, x = data.matrix(x), 
               pdf = pdf.estim,  nc = nc, graph = graph.nb, cluster.cores = g,
               tree = struct$tree, noc = struct$noc, obs.cluster = obs.cl)
    out
  }
}


comp.centr.wei.ico.2 <- function(x, mesh, nomi) {
  l <- vicinity(mesh, nomi, output='list', order=1, self=F)
  l <- sapply(1:length(l), function(i) l[[i]] <- l[[i]][l[[i]] %in% nomi])
  adj <- lapply(1:length(l), function(i) {
    sapply(1:length(l[[i]]), 
           function(x) which(nomi == l[[i]][x]))
  })
  obs <- locate(mesh, x)
  n <- sapply(1:length(adj), function(i) length(which(obs == nomi[i])))
  pos <- row.names(mesh@faceCenters) %in% nomi
  return(list(t=mesh@faceCenters[pos,], adj=adj, obs=obs, n=n))
}


pdfCluster.fermi.2 <- function(
    x, mesh, h, nomi,
    n.grid=min(round((5 + sqrt(length(nomi)))*4), length(nomi)),
    ... 
    )
  # n.grid defines the length of the grid on which evaluating the connected
  # components of the density level sets.
{
  call <- match.call()
  x <- data.matrix(x)
  if(any(sapply(x[1,],is.factor))) stop("All variables must be numeric")
  #if(any(apply(x,2, pdfCluster:::check.discrete)<10)) 
  #  warning("One or more variables look to be discrete: pdfCluster is designed for continuous data")
  if(any(is.na(x))) 
  {warning(cat("NA in object",deparse(substitute(x)), "have been omitted", "\n"));  x <- na.omit(x)}
  
  # compute the elements to be used for both density estimation and connected components
  l <- comp.centr.wei.ico.2(x, mesh, nomi)
  t <- l$t
  graph.nb <- l$adj # this is my graph
  obs <- l$obs 
  n <- l$n
  
  # density estimation
  coarse_estimate <- do.call(vmf.bskde.2, list(x = t, h = h/10, n = n))
  #print(range(coarse_estimate))
  #h_variable <- h/sqrt(coarse_estimate)
  g <- exp(mean(log(coarse_estimate)))
  h_variable <- h * (1/g * coarse_estimate)^(-0.5)
  #h_variable <- h * (coarse_estimate)^(-0.5)
  #print(range(h_variable))
  estimate <- do.call(vmf.bskde.2, list(x = t, h = as.vector(h_variable), n = n))
  #print(range(estimate))
  
  # check given arguments
  N <- nrow(t)
  if (n.grid > N) {
    warning("n.grid too large, set equal to N")
    n.grid <- min(n.grid, nrow(t))
  }
  
  # connected components
  nc <- pdfCluster:::num.con(t, estimate, graph.nb,
                             profile.grid = n.grid-2, correct=TRUE)  
  struct <- pdfCluster:::con2tree(nc, estimate)
  if (struct$bad) {
    message("No output given")
  }
  else {
    g <- struct$g
    g[struct$g == 0] <- NA
    pdf.estim <- estimate
    names(l$adj) <- nomi
    obs.cl <- sapply(1:dim(x)[1], function(i) g[which(names(l$adj) == obs[i])])
    out <- new("pdfCluster2", call = call, x = data.matrix(x), 
               pdf = pdf.estim,  nc = nc, graph = graph.nb, cluster.cores = g,
               tree = struct$tree, noc = struct$noc, obs.cluster = obs.cl)
    out
  }
}


compute_dens <- function(x, mesh, h, nomi) {
  l <- comp.centr.wei.ico.2(x, mesh, nomi)
  t <- l$t
  obs <- l$obs 
  n <- l$n
  
  if (is.null(dim(t))) t <- t(as.matrix(t))
  
  coarse_estimate <- do.call(vmf.bskde.2, list(x = t, h = h/10, n = n))
  g <- exp(mean(log(coarse_estimate)))
  h_variable <- h * (1/g * coarse_estimate)^(-0.5)
  estimate <- do.call(vmf.bskde.2, list(x = t, h = as.vector(h_variable), n = n))
  
  # I calculate the density at each observation 
  # (assigning the density of the cell to which it belongs).
  #est_obs <- sapply(1:dim(x)[1], function(i) estimate[which(nomi == obs[i])]) 
  
  return(estimate)
}


dens <- function(x, mesh, h, nomi) {
  l <- comp.centr.wei.ico.2(x, mesh, nomi)
  t <- l$t
  obs <- l$obs 
  n <- l$n
  
  if (is.null(dim(t))) t <- t(as.matrix(t)) # una sola cella
  
  coarse_estimate <- do.call(vmf.bskde.2, list(x = t, h = h/10, n = n))
  g <- exp(mean(log(coarse_estimate)))
  h_variable <- h * (1/g * coarse_estimate)^(-0.5)
  estimate <- do.call(vmf.bskde.2, list(x = t, h = as.vector(h_variable), n = n))
  
  # Calcolo la densità di ogni osservazione (assegno la densita' della cella a
  # cui appartiene).
  est_obs <- sapply(1:dim(x)[1], function(i) estimate[which(nomi == obs[i])]) 
  
  return(est_obs)
}


tpr_fpr <- function(labels1, labels2){
  n <- length(labels1)
  N <- table(labels1, labels2)
  N <- N[-nrow(N),-ncol(N)]
  r <- nrow(N)
  s <- ncol(N)
  if (r>s) {N <- t(N); r <- nrow(N); s <- ncol(N)}
  else {N <- rbind(N, matrix(0, nrow = s-r, ncol = s))}
  optimal.permutation <- solve_LSAP(N, maximum = T)
  mat <- matrix(0, nrow = r, ncol = s)
  mat[cbind(seq_along(optimal.permutation), optimal.permutation)] <- 
    N[cbind(seq_along(optimal.permutation), optimal.permutation)]
  tpr <- sum(rowSums(mat)!=0)/nrow(mat)
  fpr <- sum(colSums(mat)==0)/ncol(mat)
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


# OTHERS #######################################################################
# Rule of Thumb - Badwidth selection
#k <- Directional::vmf.mle(all)$kappa
#n <- dim(all)[1]
#h <- ( 8 * sinh(k)^2 /
#         ( k * n * ( (1 + 4 * k^2) * sinh(2 * k) -  2 * k * cosh(2 * k) ) ) )^(1/6)
#h

SEi <- function(cut, E, PSF_T){
  c0 <- NULL
  for(i in 1:dim(cut)[1]){
    if(PSF_T[i] == "0") c0[i] <- 1.56*10^-1
    if(PSF_T[i] == "1") c0[i] <- 9.64*10^-2
    if(PSF_T[i] == "2") c0[i] <- 7.02*10^-2
    if(PSF_T[i] == "3") c0[i] <- 4.97*10^-2
  }
  c1 <- NULL
  for(i in 1:dim(cut)[1]){
    if(PSF_T[i] == "0") c1[i] <- 5.70*10^-3
    if(PSF_T[i] == "1") c1[i] <- 1.78*10^-3
    if(PSF_T[i] == "2") c1[i] <- 1.07*10^-3
    if(PSF_T[i] == "3") c1[i] <- 6.13*10^-4
  }
  return(sqrt((c0*(E/100)^(-0.8))^2 + c1^2))
}


# Modified to account for bandwidth selection based on the point spread function
pdfCluster.fermi.psf <- function(
    x, mesh, energy, psf_t, nomi,
    n.grid=min(round((5 + sqrt(length(nomi)))*4), length(nomi)),
    ... 
    )
{
  call <- match.call()
  x <- data.matrix(x)
  if(any(sapply(x[1,],is.factor))) stop("All variables must be numeric")
  #if(any(apply(x,2, pdfCluster:::check.discrete)<10)) 
  #  warning("One or more variables look to be discrete: pdfCluster is designed for continuous data")
  if(any(is.na(x))) 
  {warning(cat("NA in object",deparse(substitute(x)), "have been omitted", "\n"));  x <- na.omit(x)}
  
  # compute the elements to be used for both density estimation and connected components
  l <- comp.centr.wei.ico.2(x, mesh, nomi)
  t <- l$t
  graph.nb <- l$adj # this is my graph
  obs <- l$obs 
  n <- l$n
  
  # density estimation
  h_oss <- SEi(x, energy, psf_t)
  h_variable <- sapply(1:length(nomi), 
                       function(i) {
                         if (nomi[i] %in% obs) mean(h_oss[which(obs == nomi[i])])
                       })
  h_variable <- sapply(1:length(h_variable),
                       function(i) {
                         if (is.null(h_variable[[i]])) {
                           max(unlist(h_variable)) + max(unlist(h_variable))/100 
                         } else {
                           h_variable[[i]]
                         }
                       })
  h_variable <- unlist(h_variable)
  print(h_variable)
  estimate <- do.call(vmf.bskde.2, list(x = t, h = h_variable, n = n))
  
  # check given arguments
  N <- nrow(t)
  if (n.grid > N) {
    warning("n.grid too large, set equal to N")
    n.grid <- min(n.grid, nrow(t))
  }
  
  # connected components
  nc <- pdfCluster:::num.con(t, estimate, graph.nb,
                             profile.grid = n.grid-2, correct=TRUE)  
  struct <- pdfCluster:::con2tree(nc, estimate)
  if (struct$bad) {
    message("No output given")
  }
  else {
    g <- struct$g
    g[struct$g == 0] <- NA
    pdf.estim <- estimate
    names(l$adj) <- nomi
    obs.cl <- sapply(1:dim(x)[1], function(i) g[which(names(l$adj) == obs[i])])
    out <- new("pdfCluster2", call = call, x = data.matrix(x), 
               pdf = pdf.estim,  nc = nc, graph = graph.nb, cluster.cores = g,
               tree = struct$tree, noc = struct$noc, obs.cluster = obs.cl)
    out
  }
}

# Correction
vmf.bskde.2.test = function (x, h, n, y=x){
  N <- sum(n)
  d = tcrossprod(y, x)*(1/h^2)
  argsum = (exp(d - 2*log(h) - log(2*pi) - 1/h^2)) %*% n[n>0]
  f = exp(log(argsum))/N
  f
}


pdfCluster.fermi.2.test <- function(
    x, mesh, h, nomi,
    n.grid=min(round((5 + sqrt(length(nomi)))*10), length(nomi)),
    ...
    )
  # n.grid defines the length of the grid on which evaluating the connected
  # components of the density level sets.
{
  call <- match.call()
  x <- data.matrix(x)
  if(any(sapply(x[1,],is.factor))) stop("All variables must be numeric")
  #if(any(apply(x,2, pdfCluster:::check.discrete)<10)) 
  #  warning("One or more variables look to be discrete: pdfCluster is designed for continuous data")
  if(any(is.na(x))) 
  {warning(cat("NA in object",deparse(substitute(x)), "have been omitted", "\n"));  x <- na.omit(x)}
  
  # compute the elements to be used for both density estimation and connected components
  l <- comp.centr.wei.ico.2(x, mesh, nomi)
  t <- l$t
  graph.nb <- l$adj # this is my graph
  obs <- l$obs 
  n <- l$n
  
  # density estimation
  coarse_estimate <- do.call(vmf.bskde.2.test, list(x = t[n>0,], h = h/10, n = n, y = t))
  #print(range(coarse_estimate))
  #h_variable <- h/sqrt(coarse_estimate)
  g <- exp(mean(log(coarse_estimate)))
  h_variable <- h * (1/g * coarse_estimate)^(-0.5)
  #print(range(h_variable))
  estimate <- do.call(vmf.bskde.2.test, list(x = t[n>0,], h = as.vector(h_variable)[n>0], n = n, y = t))
  #print(range(estimate))
  
  # check given arguments
  N <- nrow(t)
  if (n.grid > N) {
    warning("n.grid too large, set equal to N")
    n.grid <- min(n.grid, nrow(t))
  }
  
  # connected components
  nc <- pdfCluster:::num.con(t, estimate, graph.nb,
                             profile.grid = n.grid-2, correct=TRUE)  
  struct <- pdfCluster:::con2tree(nc, estimate)
  if (struct$bad) {
    message("No output given")
  }
  else {
    g <- struct$g
    g[struct$g == 0] <- NA
    pdf.estim <- estimate
    names(l$adj) <- nomi
    obs.cl <- sapply(1:dim(x)[1], function(i) g[which(names(l$adj) == obs[i])])
    out <- new("pdfCluster2", call = call, x = data.matrix(x), 
               pdf = pdf.estim,  nc = nc, graph = graph.nb, cluster.cores = g,
               tree = struct$tree, noc = struct$noc, obs.cluster = obs.cl)
    out
  }
}


pdfCluster.fermi.psf.test <- function(
    x, mesh, energy, psf_t, nomi,
    n.grid=min(round((5 + sqrt(length(nomi)))*4), length(nomi)),
    ...
    )
  # n.grid defines the length of the grid on which evaluating the connected
  # components of the density level sets.
{
  call <- match.call()
  x <- data.matrix(x)
  if(any(sapply(x[1,],is.factor))) stop("All variables must be numeric")
  #if(any(apply(x,2, pdfCluster:::check.discrete)<10)) 
  #  warning("One or more variables look to be discrete: pdfCluster is designed for continuous data")
  if(any(is.na(x))) 
  {warning(cat("NA in object",deparse(substitute(x)), "have been omitted", "\n"));  x <- na.omit(x)}
  
  # compute the elements to be used for both density estimation and connected components
  l <- comp.centr.wei.ico.2(x, mesh, nomi)
  t <- l$t
  graph.nb <- l$adj # this is my graph
  obs <- l$obs 
  n <- l$n
  
  # density estimation
  h_oss <- SEi(x, energy, psf_t)
  h_variable <- sapply(1:length(nomi), 
                       function(i) {
                         if (nomi[i] %in% obs){
                           mean(h_oss[which(obs == nomi[i])])
                         }
                       })
  h_variable <- unlist(h_variable)
  estimate <- do.call(vmf.bskde.2.test, list(x = t[n>0,], h = as.vector(h_variable), n = n, y = t))
  #print(range(estimate))
  
  # check given arguments
  N <- nrow(t)
  if (n.grid > N) {
    warning("n.grid too large, set equal to N")
    n.grid <- min(n.grid, nrow(t))
  }
  
  # connected components
  nc <- pdfCluster:::num.con(t, estimate, graph.nb,
                             profile.grid = n.grid-2, correct=TRUE)  
  struct <- pdfCluster:::con2tree(nc, estimate)
  if (struct$bad) {
    message("No output given")
  }
  else {
    g <- struct$g
    g[struct$g == 0] <- NA
    pdf.estim <- estimate
    names(l$adj) <- nomi
    obs.cl <- sapply(1:dim(x)[1], function(i) g[which(names(l$adj) == obs[i])])
    out <- new("pdfCluster2", call = call, x = data.matrix(x), 
               pdf = pdf.estim,  nc = nc, graph = graph.nb, cluster.cores = g,
               tree = struct$tree, noc = struct$noc, obs.cluster = obs.cl)
    out
  }
}