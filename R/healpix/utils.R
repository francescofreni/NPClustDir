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
