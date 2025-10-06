## TABLE 1 
tab_Z <- table(observed.data$Z)
tab_RD_Z <- table(observed.data$RD.obs,observed.data$Z)
tab_D <- summary(observed.data$Dobs[observed.data$RD.obs==1])
sd_D <- sd(observed.data$Dobs[observed.data$RD.obs==1])
tab_RY <- table(observed.data$RY)
tab_Y <- summary(observed.data$Y[observed.data$RY==1])
sd_Y <- sd(observed.data$Y[observed.data$RY==1])
tab_x1 <- summary((observed.data$x1-mean(observed.data$x1))/sd(observed.data$x1))
sd_x1 <- sd((observed.data$x1-mean(observed.data$x1))/sd(observed.data$x1))
tab_x2 <- table(observed.data$x2)
tab_x3 <- table(observed.data$x3)

Table1<-rbind(
  cbind(paste("Z = 1"), paste(round((tab_Z[2]/sum(tab_Z))*100,2),"% (",tab_Z[2],"/",sum(tab_Z),")",sep=""), paste("--"), paste("--"), paste("--"), paste("--"), paste("--"), paste("--")),
  cbind(paste("I(D < C)"), paste(round((tab_RD_Z[2,2]/tab_Z[2])*100,2),"% (",tab_RD_Z[2,2],"/",tab_Z[2],")",sep=""), paste("--"), paste("--"), paste("--"), paste("--"), paste("--"), paste("--")),
  cbind(paste("D*"), paste(round(tab_D[4],2)), paste(round(sd_D,2)), paste(round(tab_D[1],2)), paste(round(tab_D[2],2)), paste(round(tab_D[3],2)), paste(round(tab_D[5],2)), paste(round(tab_D[6],2))),
  cbind(paste("I(Y < C)"), paste(round((tab_RY[2]/sum(tab_Z))*100,2),"% (",tab_RY[2],"/",sum(tab_Z),")",sep=""), paste("--"), paste("--"), paste("--"), paste("--"), paste("--"), paste("--")),
  cbind(paste("Y*"), paste(round(tab_Y[4],2)), paste(round(sd_Y,2)), paste(round(tab_Y[1],2)), paste(round(tab_Y[2],2)), paste(round(tab_Y[3],2)), paste(round(tab_Y[5],2)), paste(round(tab_Y[6],2))),
  cbind(paste("X1"), paste(round(tab_x1[4],2)), paste(round(sd_x1,2)), paste(round(tab_x1[1],2)), paste(round(tab_x1[2],2)), paste(round(tab_x1[3],2)), paste(round(tab_x1[5],2)), paste(round(tab_x1[6],2))),
  cbind(paste("X2"), paste(round((tab_x2[2]/sum(tab_Z))*100,2),"% (",tab_x2[2],"/",tab_Z[2],")",sep=""), paste("--"), paste("--"), paste("--"), paste("--"), paste("--"), paste("--")),
  cbind(paste("X3"), paste(round((tab_x3[2]/sum(tab_Z))*100,2),"% (",tab_x3[2],"/",tab_Z[2],")",sep=""), paste("--"), paste("--"), paste("--"), paste("--"), paste("--"), paste("--"))
)

colnames(Table1) <- c("Variable","Mean (proportion)", "SD", "Min", "Q1", "Median", "Q3", "Max")
saveRDS(Table1, file = "tables/Table1.rds")

## TABLE 2
Table2_left<-table(observed.data$RD.obs[observed.data$Z==1],
      observed.data$RY[observed.data$Z==1])
rownames(Table2_left)<-c("I(D <= C) = 0", "I(D <= C) = 1")
colnames(Table2_left)<-c("I(Y <= C) = 0", "I(Y <= C) = 1")
Table2_left<-cbind(Table2_left,c(rowSums(Table2_left)))
Table2_left<-rbind(Table2_left,c(colSums(Table2_left)))

Table2_right<-table(observed.data$RD.obs[observed.data$Z==0],
                   observed.data$RY[observed.data$Z==0])
rownames(Table2_right)<-"I(D <= C) = 0"
colnames(Table2_right)<-c("I(Y <= C) = 0", "I(Y <= C) = 1")
Table2_right<-cbind(Table2_right,c(rowSums(Table2_right)))

saveRDS(Table2_left, file = "tables/Table2_left.rds")
saveRDS(Table2_right, file = "tables/Table2_right.rds")

## Figure 2
status<-ifelse(observed.data$RY==1,1,0)
jpeg("figures/Figure2.jpeg",width = 5000, height = 5000, units="px", res=1200)
print(autoplot(survfit(Surv(observed.data$Y,status)~observed.data$Z), xlim=c(0,35)) +
  theme_light() +
  guides(colour="none") +    
  labs(x = "Months", y = "Survival Probability") +
  labs(colour = "Z"))
dev.off()
