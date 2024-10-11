m1d<-rep(NA,dim(THETA)[1])
for(i in 1:dim(THETA)[1]){m1d[i]<-cov_distr[[i]]$m_d[1]}
m1nd<-rep(NA,dim(THETA)[1])
for(i in 1:dim(THETA)[1]){m1nd[i]<-cov_distr[[i]]$m_nd[1]}
m1ed<-rep(NA,dim(THETA)[1])
for(i in 1:dim(THETA)[1]){m1ed[i]<-cov_distr[[i]]$m_ed[1]}
m1ld<-rep(NA,dim(THETA)[1])
for(i in 1:dim(THETA)[1]){m1ld[i]<-cov_distr[[i]]$m_ld[1]}

m2d<-rep(NA,dim(THETA)[1])
for(i in 1:dim(THETA)[1]){m2d[i]<-cov_distr[[i]]$m_d[2]}
m2nd<-rep(NA,dim(THETA)[1])
for(i in 1:dim(THETA)[1]){m2nd[i]<-cov_distr[[i]]$m_nd[2]}
m2ed<-rep(NA,dim(THETA)[1])
for(i in 1:dim(THETA)[1]){m2ed[i]<-cov_distr[[i]]$m_ed[2]}
m2ld<-rep(NA,dim(THETA)[1])
for(i in 1:dim(THETA)[1]){m2ld[i]<-cov_distr[[i]]$m_ld[2]}

m5d<-rep(NA,dim(THETA)[1])
for(i in 1:dim(THETA)[1]){m5d[i]<-cov_distr[[i]]$m_d[3]}
m5nd<-rep(NA,dim(THETA)[1])
for(i in 1:dim(THETA)[1]){m5nd[i]<-cov_distr[[i]]$m_nd[3]}
m5ed<-rep(NA,dim(THETA)[1])
for(i in 1:dim(THETA)[1]){m5ed[i]<-cov_distr[[i]]$m_ed[3]}
m5ld<-rep(NA,dim(THETA)[1])
for(i in 1:dim(THETA)[1]){m5ld[i]<-cov_distr[[i]]$m_ld[3]}

par(mfrow=c(1,1))
jpeg("x1.jpeg",width=4500,height = 2864,units="px",res=500)
hist(m1nd, breaks=50,
     main="",
     xlim=c(-0.5,1), 
     ylim=c(0,15), 
     xlab=expression(X[1]), 
     ylab="",
     col="blue4", 
     density=20, 
     freq = F)
hist(m1ld, breaks=50,
     add=TRUE, 
     col="blue", 
     density=20, 
     freq = F)
hist(m1ed, 
     breaks=50,
     add=TRUE, 
     col="lightblue", 
     density=20, 
     freq = F)
legend("topright", 
       lty=c(1,1,1),
       lwd=2,
       cex=1,
       col= c("blue4", "blue","lightblue"), 
       legend=c("ND","late D", "early D"))
dev.off()

jpeg("x2.jpeg",width=4500,height = 2864,units="px",res=500)
hist(m2nd, 
     breaks=50,
     main="",
     xlab=expression(X[2]),
     ylab="",
     xlim=c(0.1,0.5), 
     ylim=c(0,30), 
     col="red4", 
     density=50, freq = F)
hist(m2ld, breaks=50,
     add=TRUE,  col="red", density=20, freq = F)
hist(m2ed, breaks=50,
     add=TRUE, col="orange", density=20, freq = F)
legend("topright", 
       cex=1,
       lty=c(1,1,1),lwd=2,col= c("red4", "red","orange"), 
       legend=c("ND","late D", "early D"))
dev.off()

jpeg("x3.jpeg",width=4500,height = 2864,units="px",res=500)
hist(m5nd, 
     breaks=50,
     main="",
     xlab=expression(X[3]),
     ylab="",
     xlim=c(0,0.45), 
     ylim=c(0,30), 
     col="green4", 
     density=20, 
     freq = F)
hist(m5ld, 
     breaks=50,
     add=TRUE,  
     col="green",
     density=20, 
     freq = F)
hist(m5ed, 
     breaks=50,
     add=TRUE, 
     col="yellow", 
     density=20, 
     freq = F)
legend("topright", 
       lty=c(1,1),
       lwd=2,
       cex=1,
       col= c("green4", "green","yellow"), 
       legend=c("ND","late D", "early D"))
dev.off()