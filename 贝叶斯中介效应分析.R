###贝叶斯中介效应分析
library(MCMCglmm)
a<-read.csv("P5S1.csv")
a <- a[!a$Populations %in% c("3", "4", "45","20"), ]
trees <- read.nexus("44species1000times.nex")
dropTip <- trees[[1]]$tip.label[which(is.na(match(trees[[1]]$tip.label, a$Species)))]
zrTrees <- lapply(trees, drop.tip, dropTip, trim.internal=T)
zrTrees <- lapply(zrTrees, makeNodeLabel, method = "number")
a$Species[which((a$Species %in% zrTrees[[1]]$tip.label) == FALSE)]
zrTrees[[1]]$tip.label[which((zrTrees[[1]]$tip.label %in% a$Species) == FALSE)]
a1 <- 1000
prior1<-list(R=list(V=diag(1),nu=0.002)
             , G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a1),
                      G2=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a1)))
INtree <- inverseA(zrTrees[[1]], nodes="TIPS")
###############降水的年间变异与SMD
#中介效应方程（1）
MEV <- a$SMD_sampling.variances
model_a.start<- MCMCglmm( SMD ~ PreVaramongyears, random = ~ Species+Populations,data=a,mev=MEV,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_a<- model_a.start
for(i in 1:1000){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=model_a$Liab[1,], R=model_a$VCV[1,3], G=list(G1=model_a$VCV[1,1],G2=model_a$VCV[1,2]))
  model_a <- MCMCglmm(SMD ~ PreVaramongyears, random = ~ Species+Populations,data=a, mev=MEV, pl=TRUE, slice=FALSE, nitt=1000, thin=1, burnin=990, prior=prior1,ginverse=list(Species=INtree$Ainv), start=start, verbose=FALSE)
  if (i > 1) {
    # 保存model_a$VCV的十行数据
    row_start <- 10 * i - 9
    row_end <- 10 * i
    model_a.start$VCV[row_start:row_end, ] <- model_a$VCV[1:10, ]
    model_a.start$Sol[row_start:row_end, ] <- model_a$Sol[1:10, ]
    model_a.start$Liab[row_start:row_end, ] <- model_a$Liab[1:10, ]
  }
}
Model_A <- model_a.start
#中介效应方程（2）
MEV <- a$BSNonCooSD^2
model_b.start<- MCMCglmm( BSNonCooMean ~ PreVaramongyears, random = ~ Species+Populations,data=a,mev=MEV,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_b<- model_b.start
for(i in 1:1000){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=model_b$Liab[1,], R=model_b$VCV[1,3], G=list(G1=model_b$VCV[1,1],G2=model_b$VCV[1,2]))
  model_b <- MCMCglmm(BSNonCooMean ~ PreVaramongyears, random = ~ Species+Populations,data=a, mev=MEV, pl=TRUE, slice=FALSE, nitt=1000, thin=1, burnin=990, prior=prior1,ginverse=list(Species=INtree$Ainv), start=start, verbose=FALSE)
  if (i > 1) {
    # 保存model_b$VCV的十行数据
    row_start <- 10 * i - 9
    row_end <- 10 * i
    model_b.start$VCV[row_start:row_end, ] <- model_b$VCV[1:10, ]
    model_b.start$Sol[row_start:row_end, ] <- model_b$Sol[1:10, ]
    model_b.start$Liab[row_start:row_end, ] <- model_b$Liab[1:10, ]
  }
}
Model_B <- model_b.start
#中介效应方程（3）
MEV <- a$SMD_sampling.variances
model_c.start<- MCMCglmm( SMD ~ BSNonCooMean+PreVaramongyears, random = ~ Species+Populations,data=a,mev=MEV,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_c<- model_c.start
for(i in 1:1000){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=model_c$Liab[1,], R=model_c$VCV[1,3], G=list(G1=model_c$VCV[1,1],G2=model_c$VCV[1,2]))
  model_c <- MCMCglmm(SMD ~ BSNonCooMean+PreVaramongyears, random = ~ Species+Populations,data=a, mev=MEV, pl=TRUE, slice=FALSE, nitt=1000, thin=1, burnin=990, prior=prior1,ginverse=list(Species=INtree$Ainv), start=start, verbose=FALSE)
  if (i > 1) {
    # 保存model_c$VCV的十行数据
    row_start <- 10 * i - 9
    row_end <- 10 * i
    model_c.start$VCV[row_start:row_end, ] <- model_c$VCV[1:10, ]
    model_c.start$Sol[row_start:row_end, ] <- model_c$Sol[1:10, ]
    model_c.start$Liab[row_start:row_end, ] <- model_c$Liab[1:10, ]
  }
}
Model_C <- model_c.start

summary(Model_A)
summary(Model_B)
summary(Model_C)
hist(Model_A$Liab)
plot(Model_A$VCV)    	 
plot(Model_A$Sol)

# 提取路径a的后验样本：PreVaramongyears对BSNonCooMean的影响
post_samples_a <-Model_B$Sol[, "PreVaramongyears"]
# 提取路径b的后验样本： BSNonCooMean对SMD的影响
post_samples_b <- Model_C$Sol[, "BSNonCooMean"]
# 计算间接效应的后验样本
indirect_effect_samples <- post_samples_a * post_samples_b
# 计算间接效应的后验总结
indirect_effect_summary <- c(mean = mean(indirect_effect_samples),
                             lower = quantile(indirect_effect_samples, 0.05),
                             upper = quantile(indirect_effect_samples, 0.95))
# 输出间接效应的后验总结
indirect_effect_summary
p_value <- sum(indirect_effect_samples >= 0) / length(indirect_effect_samples)

###############温度的年间变异与SMD
#中介效应方程（1）
MEV <- a$SMD_sampling.variances
model_a.start<- MCMCglmm( SMD ~ TmpVaramongyears, random = ~ Species+Populations,data=a,mev=MEV,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_a<- model_a.start
for(i in 1:1000){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=model_a$Liab[1,], R=model_a$VCV[1,3], G=list(G1=model_a$VCV[1,1],G2=model_a$VCV[1,2]))
  model_a <- MCMCglmm(SMD ~ TmpVaramongyears, random = ~ Species+Populations,data=a, mev=MEV, pl=TRUE, slice=FALSE, nitt=1000, thin=1, burnin=990, prior=prior1,ginverse=list(Species=INtree$Ainv), start=start, verbose=FALSE)
  if (i > 1) {
    # 保存model_a$VCV的十行数据
    row_start <- 10 * i - 9
    row_end <- 10 * i
    model_a.start$VCV[row_start:row_end, ] <- model_a$VCV[1:10, ]
    model_a.start$Sol[row_start:row_end, ] <- model_a$Sol[1:10, ]
    model_a.start$Liab[row_start:row_end, ] <- model_a$Liab[1:10, ]
  }
}
Model_A1 <- model_a.start
#中介效应方程（2）
MEV <- a$BSNonCooSD^2
model_b.start<- MCMCglmm( BSNonCooMean ~ TmpVaramongyears, random = ~ Species+Populations,data=a,mev=MEV,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_b<- model_b.start
for(i in 1:1000){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=model_b$Liab[1,], R=model_b$VCV[1,3], G=list(G1=model_b$VCV[1,1],G2=model_b$VCV[1,2]))
  model_b <- MCMCglmm(BSNonCooMean ~ TmpVaramongyears, random = ~ Species+Populations,data=a, mev=MEV, pl=TRUE, slice=FALSE, nitt=1000, thin=1, burnin=990, prior=prior1,ginverse=list(Species=INtree$Ainv), start=start, verbose=FALSE)
  if (i > 1) {
    # 保存model_b$VCV的十行数据
    row_start <- 10 * i - 9
    row_end <- 10 * i
    model_b.start$VCV[row_start:row_end, ] <- model_b$VCV[1:10, ]
    model_b.start$Sol[row_start:row_end, ] <- model_b$Sol[1:10, ]
    model_b.start$Liab[row_start:row_end, ] <- model_b$Liab[1:10, ]
  }
}
Model_B1 <- model_b.start
#中介效应方程（3）
MEV <- a$SMD_sampling.variances
model_c.start<- MCMCglmm( SMD ~ BSNonCooMean+TmpVaramongyears, random = ~ Species+Populations,data=a,mev=MEV,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_c<- model_c.start
for(i in 1:1000){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=model_c$Liab[1,], R=model_c$VCV[1,3], G=list(G1=model_c$VCV[1,1],G2=model_c$VCV[1,2]))
  model_c <- MCMCglmm(SMD ~ BSNonCooMean+TmpVaramongyears, random = ~ Species+Populations,data=a, mev=MEV, pl=TRUE, slice=FALSE, nitt=1000, thin=1, burnin=990, prior=prior1,ginverse=list(Species=INtree$Ainv), start=start, verbose=FALSE)
  if (i > 1) {
    # 保存model_c$VCV的十行数据
    row_start <- 10 * i - 9
    row_end <- 10 * i
    model_c.start$VCV[row_start:row_end, ] <- model_c$VCV[1:10, ]
    model_c.start$Sol[row_start:row_end, ] <- model_c$Sol[1:10, ]
    model_c.start$Liab[row_start:row_end, ] <- model_c$Liab[1:10, ]
  }
}
Model_C1 <- model_c.start

summary(Model_A1)
summary(Model_B1)
summary(Model_C1)
hist(Model_A1$Liab)
plot(Model_A1$VCV)    	 
plot(Model_A1$Sol)
post_samples_a <-Model_B1$Sol[, "TmpVaramongyears"]
post_samples_b <- Model_C1$Sol[, "BSNonCooMean"]
indirect_effect_samples <- post_samples_a * post_samples_b
indirect_effect_summary <- c(mean = mean(indirect_effect_samples),
                             lower = quantile(indirect_effect_samples, 0.05),
                             upper = quantile(indirect_effect_samples, 0.95))
indirect_effect_summary
p_value <- sum(indirect_effect_samples >= 0) / length(indirect_effect_samples)


##环境的不确定(温度年间变异)与合作繁殖效应（提高成体年存活率）
library(MCMCglmm)
library(dplyr)
Calc.lnRR <- function(P1, P2){
  log(P1) - log(P2)
}
Calc.var.lnRR <- function(P1, N1, P2, N2){
  ((1-P1)/(N1*P1)) + ((1-P2)/(N2*P2))
}

# 数据预处理
MA_dataset <- read.csv("P5S3.csv") %>%
  mutate(
    female.pair....survial = female.pair....survial / 100,
    female.group.....survial = female.group.....survial / 100,
    lnRR = Calc.lnRR(female.group.....survial, female.pair....survial),
    var.lnRR = Calc.var.lnRR(female.group.....survial, N.female.group, female.pair....survial, N.male.pair)
  )
trees <- read.nexus("25species1000times.nex")
dropTip <- trees[[1]]$tip.label[which(is.na(match(trees[[1]]$tip.label, MA_dataset$Species)))]
zrTrees <- lapply(trees, drop.tip, dropTip, trim.internal=TRUE)
zrTrees <- lapply(zrTrees, makeNodeLabel, method = "number")
prior1 <- list(R=list(V=diag(1),nu=0.002), G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))
INtree <- inverseA(zrTrees[[1]], nodes="TIPS")
#中介分析方程（1）
MEV<-MA_dataset$var.lnRR
model_a.start <- MCMCglmm(  lnRR ~ TmpVaramongyears,   random = ~ Species,  data=MA_dataset,mev=MEV,  pl=TRUE,  slice=FALSE,  nitt=103000,  thin=10,  burnin=3000,  prior=prior1,  ginverse=list(Species=INtree$Ainv),  verbose=FALSE)
model_a<- model_a.start
for(i in 1:1000){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=model_a$Liab[1,], R=model_a$VCV[1,2], G=list(G1=model_a$VCV[1,1]))
  model_a <- MCMCglmm(lnRR ~ TmpVaramongyears, random = ~ Species,data=MA_dataset, mev=MEV, pl=TRUE, slice=FALSE, nitt=1000, thin=1, burnin=990, prior=prior1,ginverse=list(Species=INtree$Ainv), start=start, verbose=FALSE)
  if (i > 1) {
    # 保存model_a$VCV的十行数据
    row_start <- 10 * i - 9
    row_end <- 10 * i
    model_a.start$VCV[row_start:row_end, ] <- model_a$VCV[1:10, ]
    model_a.start$Sol[row_start:row_end, ] <- model_a$Sol[1:10, ]
    model_a.start$Liab[row_start:row_end, ] <- model_a$Liab[1:10, ]
  }
}
Model_A2 <- model_a.start
#中介分析方程（2）
model_b.start <- MCMCglmm(   scale(male.pair.....survial) ~ TmpVaramongyears,   random = ~ Species,  data=MA_dataset,   pl=TRUE,  slice=FALSE,  nitt=103000,  thin=10,  burnin=3000,  prior=prior1,  ginverse=list(Species=INtree$Ainv),  verbose=FALSE)
model_b<- model_b.start
for(i in 1:1000){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=model_b$Liab[1,], R=model_b$VCV[1,2], G=list(G1=model_b$VCV[1,1]))
  model_b <- MCMCglmm(scale(male.pair.....survial) ~ TmpVaramongyears, random = ~ Species,data=MA_dataset,  pl=TRUE, slice=FALSE, nitt=1000, thin=1, burnin=990, prior=prior1,ginverse=list(Species=INtree$Ainv), start=start, verbose=FALSE)
  if (i > 1) {
    # 保存model_b$VCV的十行数据
    row_start <- 10 * i - 9
    row_end <- 10 * i
    model_b.start$VCV[row_start:row_end, ] <- model_b$VCV[1:10, ]
    model_b.start$Sol[row_start:row_end, ] <- model_b$Sol[1:10, ]
    model_b.start$Liab[row_start:row_end, ] <- model_b$Liab[1:10, ]
  }
}
Model_B2 <- model_b.start

#中介分析方程（3）
MEV<-MA_dataset$var.lnRR
model_c.start <- MCMCglmm(  lnRR ~ scale(male.pair.....survial)+TmpVaramongyears ,   random = ~ Species,  data=MA_dataset,  mev=MEV,  pl=TRUE,  slice=FALSE,  nitt=103000,  thin=10,  burnin=3000,  prior=prior1,  ginverse=list(Species=INtree$Ainv),  verbose=FALSE)
model_c<- model_c.start
for(i in 1:1000){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=model_c$Liab[1,], R=model_c$VCV[1,2], G=list(G1=model_c$VCV[1,1]))
  model_c <- MCMCglmm(lnRR ~ scale(male.pair.....survial)+TmpVaramongyears ,   random = ~ Species,  data=MA_dataset,  mev=MEV,  pl=TRUE, slice=FALSE, nitt=1000, thin=1, burnin=990, prior=prior1,ginverse=list(Species=INtree$Ainv), start=start, verbose=FALSE)
  if (i > 1) {
    # 保存model_c$VCV的十行数据
    row_start <- 10 * i - 9
    row_end <- 10 * i
    model_c.start$VCV[row_start:row_end, ] <- model_c$VCV[1:10, ]
    model_c.start$Sol[row_start:row_end, ] <- model_c$Sol[1:10, ]
    model_c.start$Liab[row_start:row_end, ] <- model_c$Liab[1:10, ]
  }
}
Model_C2 <- model_c.start

summary(Model_A2)
summary(Model_B2)
summary(Model_C2)
hist(Model_A2$Liab)
plot(Model_A2$VCV)    	 
plot(Model_A2$Sol)
post_samples_a <-Model_B2$Sol[, "TmpVaramongyears"]
post_samples_b <- Model_C2$Sol[, "scale(male.pair.....survial)"]
indirect_effect_samples <- post_samples_a * post_samples_b
indirect_effect_summary <- c(mean = mean(indirect_effect_samples),
                             lower = quantile(indirect_effect_samples, 0.05),
                             upper = quantile(indirect_effect_samples, 0.95))
indirect_effect_summary
p_value <- sum(indirect_effect_samples <= 0) / length(indirect_effect_samples)
##环境的不确定(降水年间变异)与合作繁殖效应（提高成体年存活率）
#中介分析方程（1）
MEV<-MA_dataset$var.lnRR
model_a.start <- MCMCglmm(  lnRR ~ PreVaramongyears,   random = ~ Species,  data=MA_dataset,mev=MEV,  pl=TRUE,  slice=FALSE,  nitt=103000,  thin=10,  burnin=3000,  prior=prior1,  ginverse=list(Species=INtree$Ainv),  verbose=FALSE)
model_a<- model_a.start
for(i in 1:1000){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=model_a$Liab[1,], R=model_a$VCV[1,2], G=list(G1=model_a$VCV[1,1]))
  model_a <- MCMCglmm(lnRR ~ PreVaramongyears, random = ~ Species,data=MA_dataset, mev=MEV, pl=TRUE, slice=FALSE, nitt=1000, thin=1, burnin=990, prior=prior1,ginverse=list(Species=INtree$Ainv), start=start, verbose=FALSE)
  if (i > 1) {
    # 保存model_a$VCV的十行数据
    row_start <- 10 * i - 9
    row_end <- 10 * i
    model_a.start$VCV[row_start:row_end, ] <- model_a$VCV[1:10, ]
    model_a.start$Sol[row_start:row_end, ] <- model_a$Sol[1:10, ]
    model_a.start$Liab[row_start:row_end, ] <- model_a$Liab[1:10, ]
  }
}
Model_A3 <- model_a.start
#中介分析方程（2）
model_b.start <- MCMCglmm(   scale(male.pair.....survial) ~ PreVaramongyears,   random = ~ Species,  data=MA_dataset,   pl=TRUE,  slice=FALSE,  nitt=103000,  thin=10,  burnin=3000,  prior=prior1,  ginverse=list(Species=INtree$Ainv),  verbose=FALSE)
model_b<- model_b.start
for(i in 1:1000){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=model_b$Liab[1,], R=model_b$VCV[1,2], G=list(G1=model_b$VCV[1,1]))
  model_b <- MCMCglmm(scale(male.pair.....survial) ~ PreVaramongyears, random = ~ Species,data=MA_dataset,  pl=TRUE, slice=FALSE, nitt=1000, thin=1, burnin=990, prior=prior1,ginverse=list(Species=INtree$Ainv), start=start, verbose=FALSE)
  if (i > 1) {
    # 保存model_b$VCV的十行数据
    row_start <- 10 * i - 9
    row_end <- 10 * i
    model_b.start$VCV[row_start:row_end, ] <- model_b$VCV[1:10, ]
    model_b.start$Sol[row_start:row_end, ] <- model_b$Sol[1:10, ]
    model_b.start$Liab[row_start:row_end, ] <- model_b$Liab[1:10, ]
  }
}
Model_B3 <- model_b.start

#中介分析方程（3）
MEV<-MA_dataset$var.lnRR
model_c.start <- MCMCglmm(  lnRR ~ scale(male.pair.....survial)+PreVaramongyears ,   random = ~ Species,  data=MA_dataset,  mev=MEV,  pl=TRUE,  slice=FALSE,  nitt=103000,  thin=10,  burnin=3000,  prior=prior1,  ginverse=list(Species=INtree$Ainv),  verbose=FALSE)
model_c<- model_c.start
for(i in 1:1000){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=model_c$Liab[1,], R=model_c$VCV[1,2], G=list(G1=model_c$VCV[1,1]))
  model_c <- MCMCglmm(lnRR ~ scale(male.pair.....survial)+PreVaramongyears ,   random = ~ Species,  data=MA_dataset,  mev=MEV,  pl=TRUE, slice=FALSE, nitt=1000, thin=1, burnin=990, prior=prior1,ginverse=list(Species=INtree$Ainv), start=start, verbose=FALSE)
  if (i > 1) {
    # 保存model_c$VCV的十行数据
    row_start <- 10 * i - 9
    row_end <- 10 * i
    model_c.start$VCV[row_start:row_end, ] <- model_c$VCV[1:10, ]
    model_c.start$Sol[row_start:row_end, ] <- model_c$Sol[1:10, ]
    model_c.start$Liab[row_start:row_end, ] <- model_c$Liab[1:10, ]
  }
}
Model_C3 <- model_c.start

summary(Model_A3)
summary(Model_B3)
summary(Model_C3)
hist(Model_A3$Liab)
plot(Model_A3$VCV)    	 
plot(Model_A3$Sol)
post_samples_a <-Model_B3$Sol[, "PreVaramongyears"]
post_samples_b <- Model_C3$Sol[, "scale(male.pair.....survial)"]
indirect_effect_samples <- post_samples_a * post_samples_b
indirect_effect_summary <- c(mean = mean(indirect_effect_samples),
                             lower = quantile(indirect_effect_samples, 0.05),
                             upper = quantile(indirect_effect_samples, 0.95))
indirect_effect_summary
p_value <- sum(indirect_effect_samples >= 0) / length(indirect_effect_samples)


















