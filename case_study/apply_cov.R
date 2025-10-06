apply_cov <- function(x){
  data<-x$data
  ID<-x$ID
  D<-x$D
  Y<-data$Y
  RY<-data$RY
  Z<-data$Z
  xs<-as.matrix(observed.data[,grepl("x",names(observed.data))])

  n_cov<-dim(xs)[2]
  
  m_d<-NULL
  for(cov in 1:n_cov){
    m_d[cov]<-mean(xs[ID==0,cov])
  }

  m_nd<-NULL
  for(cov in 1:n_cov){
    m_nd[cov]<-mean(xs[ID==1,cov])
  }
  
  m_ed<-NULL
  for(cov in 1:n_cov){
    m_ed[cov]<-mean(xs[ID==0&D<median(D[ID==0]),cov])
  }
  
  m_ld<-NULL
  for(cov in 1:n_cov){
    m_ld[cov]<-mean(xs[ID==0&D>=median(D[ID==0]),cov])
  }
  
  s1_nd<-median(Y[RY==1&ID==1&Z==1])
  s1_ed<-median(Y[RY==1&ID==0&D<median(D[ID==0])&Z==1])
  s1_ld<-median(Y[RY==1&ID==0&D>=median(D[ID==0])&Z==1])
  s0_nd<-median(Y[RY==1&ID==1&Z==0])
  s0_ed<-median(Y[RY==1&ID==0&D<median(D[ID==0])&Z==0])
  s0_ld<-median(Y[RY==1&ID==0&D>=median(D[ID==0])&Z==0])
  
  list(m_d=m_d, 
       m_nd=m_nd, 
       m_ed=m_ed,
       m_ld=m_ld,
       s1_nd=s1_nd,
       s1_ed=s1_ed,
       s1_ld=s1_ld,
       s0_nd=s0_nd,
       s0_ed=s0_ed,
       s0_ld=s0_ld)
}