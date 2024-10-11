apply_cov <- function(x){
  ID<-x$ID
  D<-x$D
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
  
  list(m_d=m_d, 
       m_nd=m_nd, 
       m_ed=m_ed,
       m_ld=m_ld)
}