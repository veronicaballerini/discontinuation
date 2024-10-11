ace_d<-data.frame(values=rbind(as.matrix(colMeans(ace.d))),
                  CI_lower=as.matrix(hdi(ace.d)[1,]),
                  CI_upper=as.matrix(hdi(ace.d)[2,]),
                  Time=seq(1,10))

jpeg("ACE_d.jpeg",width = 5000, height = 2500, units="px", res=1200)
ggplot(ace_d, aes(x=Time, y=values))+ 
  geom_line(linewidth=.8,) +
  theme_light()+
  labs(x = "Discontinuation time d (months)", y = "") +
  geom_ribbon(aes(ymin=CI_lower, ymax=CI_upper) ,fill="lightblue4", alpha=0.2)+
  theme(legend.position = "none") +
  scale_linetype_manual(values=c("solid",
                                 "dotted"))+
  coord_cartesian(ylim=c(-5, 8),xlim=c(1,4))
dev.off()

mdcend<-apply(dce.nd,2,mean)
hdidcend<-as.matrix(apply(dce.nd,2,hdi95))
dce_nd<-data.frame(values=as.matrix(mdcend),
                   CI_lower=as.matrix(hdidcend[1,]),
                   CI_upper=as.matrix(hdidcend[2,]),
                   Time=seq(1,30))

jpeg("DCE_ND.jpeg",width=5000,height = 2500,units="px",res=1200)
ggplot(dce_nd, aes(x=Time, y=values))+ 
  theme_light()+
  geom_line(size=0.8) +
  labs(x = "Months", y = "") +
  geom_ribbon(aes(ymin=CI_lower, ymax=CI_upper) ,fill="lightblue4", alpha=0.2)+
  theme(legend.position = "none") +
  scale_linetype_manual(values=c("solid",
                                 "dotted"))
dev.off()

dce_d1<-data.frame(values=as.matrix(colMeans(dce.d_1)),
                   CI_lower=as.matrix(hdi(dce.d_1)[1,]),
                   CI_upper=as.matrix(hdi(dce.d_1)[2,]),
                   Time=seq(1,30))

dce_d2<-data.frame(values=as.matrix(colMeans(dce.d_2)),
                   CI_lower=as.matrix(hdi(dce.d_2)[1,]),
                   CI_upper=as.matrix(hdi(dce.d_2)[2,]),
                   Time=seq(1,30))

dce_d3<-data.frame(values=as.matrix(colMeans(dce.d_3)),
                   CI_lower=as.matrix(hdi(dce.d_3)[1,]),
                   CI_upper=as.matrix(hdi(dce.d_3)[2,]),
                   Time=seq(1,30))

dce_d4<-data.frame(values=as.matrix(colMeans(dce.d_4)),
                   CI_lower=as.matrix(hdi(dce.d_4)[1,]),
                   CI_upper=as.matrix(hdi(dce.d_4)[2,]),
                   Time=seq(1,30))

dce_d5<-data.frame(values=as.matrix(colMeans(dce.d_5)),
                   CI_lower=as.matrix(hdi(dce.d_5)[1,]),
                   CI_upper=as.matrix(hdi(dce.d_5)[2,]),
                   Time=seq(1,30))

dce_d<-data.frame(DCE_d=rbind(as.matrix(colMeans(dce.d_1)),
                              as.matrix(colMeans(dce.d_2)),
                              as.matrix(colMeans(dce.d_3)),
                              as.matrix(colMeans(dce.d_4))),
                  CI_lower=rbind(as.matrix(apply(dce.d_1,2,hdi)[1,]),
                                 as.matrix(apply(dce.d_2,2,hdi)[1,]),
                                 as.matrix(apply(dce.d_3,2,hdi)[1,]),
                                 as.matrix(apply(dce.d_4,2,hdi)[1,])),
                  CI_upper=rbind(as.matrix(apply(dce.d_1,2,hdi)[2,]),
                                 as.matrix(apply(dce.d_2,2,hdi)[2,]),
                                 as.matrix(apply(dce.d_3,2,hdi)[2,]),
                                 as.matrix(apply(dce.d_4,2,hdi)[2,])),
                  Discontinuation=rep(c("d=1",
                              "d=2",
                              "d=3",
                              "d=4"),
                            each=30),
                  Time=seq(1,30))

jpeg("DCE_d.jpeg",width=5000,height = 2500,units="px",res=1200)
ggplot(dce_d, aes(x=Time, y=DCE_d, col=Discontinuation))+ 
  theme_light()+
  geom_line(size=1) +
  labs(x = "Time", y = "") +
  geom_ribbon(aes(ymin=CI_lower, ymax=CI_upper) ,fill="gray", alpha=0.2)
dev.off()
