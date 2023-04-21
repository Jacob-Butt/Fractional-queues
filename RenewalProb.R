#install.packages("MittagLeffleR")

##########################################
# This code is used to numerically compute P(X_alpha > X_beta) over a two dimensional matrix of values for (alpha, beta) in (0,1).
# This calculation is performed by directly integrating the relevant pdf and survival funciton
###########################################

n = 100 # matrix size
nsteps = 100000 #number of steps for numerical integral
ntick = 4

#################
# parameters that set the range of alpha, beta values to consider, x corresponding to alpha and y to beta
###################
min.x = 0.1 
max.x = 1
min.y = 0.1
max.y = 1

#######################################
# These can be altered to change the viewing range of the plot, should you use it
##########################################

#max.z=0.503
#min.z=0.497 

#####################################
# matrix initialisation
###################################
M <- matrix(,nrow = n, ncol = n) # matrix of probability values
X <- seq(min.x,max.x, length.out=n) #vector of alpha values
Y <- seq(min.y,max.y, length.out=n) #vector of beta values


alpha1 = 0
alpha2 = 0 #dummy variables that will be changed late

lambda = 1 # value of lambda for X_alpha
mu = 1 # value of mu for X_beta

################################
# This section fills out the matrix by first defining the relevant function to be integrated and performing this through all values 
# of the matrix M
################################

integrand <- function(x){
  return(dml(x, alpha, lambda) * pml(x, beta, mu))
}

k = 1

while(k <= n){
  alpha = X[k]
  l=1
  while(l <= n){
    beta = Y[l]
    x = integrate(integrand, 0, Inf)
    M[k,l] <- x$value
    if(x$value > max.z){
    }
    if(x$value < min.z){
    }
    l = l + 1
  }
  k = k + 1
}

write.table(M, file="mlProbMatrix.txt", row.names=FALSE, col.names=FALSE) # returns the values of the probability matrix to a text file for later use


###############################################
# This last section of code can be used to plot the matrix M, though one might prefer to use the above text file to plot in another language.
# If using, remember to comment back in max.z and min.z from earlier.
##############################################


# Zaxis <- c(min.z, max.z)
# Xaxis <- seq(min.x,max.x, length.out=ntick) 
# Yaxis <- seq(min.y,max.y, length.out=ntick)
# 
# 
# par(mar=c(5,6,4,1)+.2)
# pmat <- persp3D( x=X, y=Y, z=M,zlim = c(min.z, max.z), ticktype="detailed", theta=40, phi = 30, xlab="alpha", ylab="beta", zlab="P", box = TRUE, axes=FALSE)
# 
# lines(trans3d(Xaxis, min.y, min.z, pmat) , col="black")
# lines(trans3d(max.x, Yaxis, min.z, pmat) , col="black")
# lines(trans3d(min.x, min.y, Zaxis, pmat) , col="black")
# 
# 
# tick.start <- trans3d(Xaxis, min.y, min.z, pmat)
# tick.end <- trans3d(Xaxis, (min.y - 0.03), min.z, pmat)
# segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)
# 
# 
# tick.start <- trans3d(max.x, Yaxis, min.z, pmat)
# tick.end <- trans3d(max.x + 0.03, Yaxis, min.z, pmat)
# segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)
# 
# tick.start <- trans3d(min.x, min.y, Zaxis, pmat)
# tick.end <- trans3d(min.x, (min.y - 0.03), Zaxis, pmat)
# segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)
# 
# labels <- seq(min.x, max.x, length.out = ntick)
# label.pos <- trans3d(Xaxis, (min.y - 0.08), min.z, pmat)
# text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), srt=270, cex=1.7)
# 
# labels <- seq(max.y, min.y, length.out = ntick)
# label.pos <- trans3d((max.x + 0.08), Yaxis, min.z, pmat)
# text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), srt=270, cex=1.7)
# 
# x.axisLabel<-expression(alpha)
# xLabel.pos <- trans3d(0.55, (min.y - 0.25), min.z, pmat)
# text(xLabel.pos$x, xLabel.pos$y, labels=x.axisLabel, srt=320, adj=c(1, NA), cex=2.1)
# 
# y.axisLabel<-expression(beta)
# yLabel.pos <- trans3d((max.x + 0.25),0.55, min.z, pmat)
# text(yLabel.pos$x, yLabel.pos$y, labels=y.axisLabel, srt=320, adj=c(1, NA), cex=2.1)
# 
# mid.z = (max.z - min.z)/2
# z.axisLabel<-expression(plain(P)(X[alpha] < X[beta]))
# zLabel.pos <- trans3d((min.x - 0.25),min.y, 0.5 , pmat)
# text(zLabel.pos$x, zLabel.pos$y, labels=z.axisLabel, srt=100, adj=c(1, NA), cex=2.1)
# 
# h = findobj('tag', 'Colorbar');
# pos = get(h, 'position')
# 
# 
# set(h, 'position', pos - 3)
