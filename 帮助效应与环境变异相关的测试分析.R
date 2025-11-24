library(MCMCglmm)
library(dplyr)
library(ape)

# 函数定义
Calc.lnRR <- function(P1, P2){
  log(P1) - log(P2)
}

Calc.var.lnRR <- function(P1, N1, P2, N2){
  ((1-P1)/(N1*P1)) + ((1-P2)/(N2*P2))
}

# 数据预处理
MA_dataset <- read.csv("P5DA.csv") %>%
  mutate(
    female.pair....survial = female.pair....survial / 100,
    female.group.....survial = female.group.....survial / 100,
    lnRR = Calc.lnRR(female.group.....survial, female.pair....survial),
    var.lnRR = Calc.var.lnRR(female.group.....survial, N.female.group, female.pair....survial, N.male.pair)
  )

# 操纵树
trees <- read.nexus("44species1000times.nex")
dropTip <- trees[[1]]$tip.label[which(is.na(match(trees[[1]]$tip.label, MA_dataset$Species)))]
zrTrees <- lapply(trees, drop.tip, dropTip, trim.internal=TRUE)
zrTrees <- lapply(zrTrees, makeNodeLabel, method = "number")

# MCMCglmm 模型
prior1 <- list(R=list(V=diag(1),nu=0.002), G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))
INtree <- inverseA(zrTrees[[1]], nodes="TIPS")

fixedModelEnv1 <- MCMCglmm(
  scale(lnRR) ~ scale(PreVaramongyears), 
  random = ~ Species,
  data=MA_dataset,
  mev=MA_dataset$var.lnRR,
  pl=TRUE,
  slice=FALSE,
  nitt=103000,
  thin=10,
  burnin=3000,
  prior=prior1,
  ginverse=list(Species=INtree$Ainv),
  verbose=FALSE
)

summary(fixedModelEnv1)


###################
library(MCMCglmm)
MA_dataset<-read.csv("P5DA.csv")
colnames(MA_dataset)
Calc.lnRR <- function(P1, P2){ES <- log(P1) - log(P2); return(ES)}
Calc.var.lnRR <- function(P1, N1, P2, N2){S2 <- ((1-P1)/(N1*P1)) + ((1-P2)/(N2*P2)); return(S2)}
MA_dataset$FePairS <- MA_dataset$FePairS/100
MA_dataset$FeGroupS <- MA_dataset$FeGroupS/100
# create effect sizes
MA_dataset$FelnRR <- Calc.lnRR(MA_dataset$FeGroupS, MA_dataset$FePairS)
MA_dataset$var.FelnRR <- Calc.var.lnRR(MA_dataset$FeGroupS, MA_dataset$FeGroupSN, MA_dataset$FePairS, MA_dataset$FePairSN)
MA_dataset$MalnRR <- Calc.lnRR(MA_dataset$MaGroupS, MA_dataset$MaPairS)
MA_dataset$var.MalnRR <- Calc.var.lnRR(MA_dataset$MaGroupS, MA_dataset$MaGroupSN, MA_dataset$MaPairS, MA_dataset$MaPairSN)



a_new<-MA_dataset
a_new<-a_new[-c(2,7),] 
trees <- read.nexus("44species1000times.nex")
dropTip <- trees[[1]]$tip.label[which(is.na(match(trees[[1]]$tip.label, a_new$Species)))]
zrTrees <- lapply(trees, drop.tip, dropTip, trim.internal=T)
zrTrees <- lapply(zrTrees, makeNodeLabel, method = "number")
a_new$Species[which((a_new$Species %in% zrTrees[[1]]$tip.label) == FALSE)]
zrTrees[[1]]$tip.label[which((zrTrees[[1]]$tip.label %in% a_new$Species) == FALSE)]
colnames(a_new)
library(dplyr)
# 从 a_new 中选择指定的变量
a_new_numeric <- a_new %>%
  select(BSCooMean, BSCooSD, BSNonCooMean, BSNonCooSD, 
         FePairS, FeGroupS, MaPairS, MaGroupS, SMD, 
         DOV, lnCVR, PreVaramongyears, TmpVaramongyears, 
         AnMeanTmp, TmpSeason, TmpSeasonCV, AnnualPre, 
         PreSeasonCV, FelnRR, MalnRR)

# 查看新数据框的前几行
head(a_new_numeric)


MEV <- a_new$SMDVar
a1 <- 1000
prior1<-list(R=list(V=diag(1),nu=0.002)
             , G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a1)))


INtree <- inverseA(zrTrees[[1]], nodes="TIPS")
model_a<- MCMCglmm( scale(MaPairS) ~ TmpVaramongyears, random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_b<- MCMCglmm( SMD ~ TmpVaramongyears, random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_c<- MCMCglmm( SMD ~ scale(MaPairS)+TmpVaramongyears, random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model<- MCMCglmm( SMD ~ scale(MaPairS), random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
summary(model)

summary(model_a)
hist(model_a$Liab)
plot(model_a$VCV)    	 
plot(model_a$Sol)
summary(model_b)
hist(model_b$Liab)
plot(model_b$VCV)    	 
plot(model_b$Sol)
summary(model_c)
hist(model_c$Liab)
plot(model_c$VCV)    	 
plot(model_c$Sol)

# 提取路径a的后验样本：TmpVaramongyears对scale(MaPairS)的影响
post_samples_a <-model_a$Sol[, "TmpVaramongyears"]
# 提取路径b的后验样本：scale(MaPairS)对SMD的影响
post_samples_b <- model_c$Sol[, "scale(MaPairS)"]
# 计算间接效应的后验样本
indirect_effect_samples <- post_samples_a * post_samples_b
# 计算间接效应的后验总结
indirect_effect_summary <- c(mean = mean(indirect_effect_samples),
                             lower = quantile(indirect_effect_samples, 0.05),
                             upper = quantile(indirect_effect_samples, 0.95))
# 输出间接效应的后验总结
indirect_effect_summary


model_as<- MCMCglmm( scale(MaPairS) ~ TmpVaramongyears, random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_bs<- MCMCglmm( MalnRR ~ TmpVaramongyears, random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_cs<- MCMCglmm( MalnRR ~ scale(MaPairS)+TmpVaramongyears, random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
summary(model_as)
hist(model_as$Liab)
plot(model_as$VCV)    	 
plot(model_as$Sol)
summary(model_bs)
hist(model_bs$Liab)
plot(model_bs$VCV)    	 
plot(model_bs$Sol)
summary(model_cs)
hist(model_cs$Liab)
plot(model_cs$VCV)    	 
plot(model_cs$Sol)
post_samples_as <-model_as$Sol[, "TmpVaramongyears"]
post_samples_bs <- model_cs$Sol[, "scale(MaPairS)"]
# 计算间接效应的后验样本
indirect_effect_samples <- post_samples_as * post_samples_bs
# 计算间接效应的后验总结
indirect_effect_summary <- c(mean = mean(indirect_effect_samples),
                             lower = quantile(indirect_effect_samples, 0.05),
                             upper = quantile(indirect_effect_samples, 0.95))
# 输出间接效应的后验总结
indirect_effect_summary



model_ay<- MCMCglmm( SMD ~ TmpVaramongyears, random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_by<- MCMCglmm( MalnRR ~ TmpVaramongyears, random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_cy<- MCMCglmm( MalnRR ~ SMD+TmpVaramongyears, random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
summary(model_ay)
hist(model_ay$Liab)
plot(model_ay$VCV)    	 
plot(model_ay$Sol)
summary(model_by)
hist(model_by$Liab)
plot(model_by$VCV)    	 
plot(model_by$Sol)
summary(model_cy)
hist(model_cy$Liab)
plot(model_cy$VCV)    	 
plot(model_cy$Sol)
post_samples_ay <-model_ay$Sol[, "TmpVaramongyears"]
post_samples_by <- model_cy$Sol[, "SMD"]
# 计算间接效应的后验样本
indirect_effect_samples <- post_samples_ay * post_samples_by
# 计算间接效应的后验总结
indirect_effect_summary <- c(mean = mean(indirect_effect_samples),
                             lower = quantile(indirect_effect_samples, 0.05),
                             upper = quantile(indirect_effect_samples, 0.95))
# 输出间接效应的后验总结
indirect_effect_summary


prior1<-list(R=list(V=diag(1),nu=0.002)
             , G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a1)))


INtree <- inverseA(zrTrees[[1]], nodes="TIPS")
model_a<- MCMCglmm( BSNonCooMean ~ TmpVaramongyears, random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_b<- MCMCglmm( SMD ~ BSNonCooMean, random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_c<- MCMCglmm( SMD ~ BSNonCooMean+TmpVaramongyears, random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)

summary(model_a)
hist(model_a$Liab)
plot(model_a$VCV)    	 
plot(model_a$Sol)
summary(model_b)
hist(model_b$Liab)
plot(model_b$VCV)    	 
plot(model_b$Sol)
summary(model_c)
hist(model_c$Liab)
plot(model_c$VCV)    	 
plot(model_c$Sol)

# 提取路径a的后验样本：TmpVaramongyears对BSNonCooMean的影响
post_samples_a <-model_a$Sol[, "TmpVaramongyears"]
# 提取路径b的后验样本： BSNonCooMean对SMD的影响
post_samples_b <- model_c$Sol[, "BSNonCooMean"]
# 计算间接效应的后验样本
indirect_effect_samples <- post_samples_a * post_samples_b
# 计算间接效应的后验总结
indirect_effect_summary <- c(mean = mean(indirect_effect_samples),
                             lower = quantile(indirect_effect_samples, 0.05),
                             upper = quantile(indirect_effect_samples, 0.95))
# 输出间接效应的后验总结
indirect_effect_summary


model_a<- MCMCglmm( BSNonCooMean ~ TmpVaramongyears, random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_b<- MCMCglmm( MalnRR ~ BSNonCooMean, random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
model_c<- MCMCglmm( MalnRR ~ BSNonCooMean+TmpVaramongyears, random = ~ Species,data=a_new,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)

summary(model_a)
hist(model_a$Liab)
plot(model_a$VCV)    	 
plot(model_a$Sol)
summary(model_b)
hist(model_b$Liab)
plot(model_b$VCV)    	 
plot(model_b$Sol)
summary(model_c)
hist(model_c$Liab)
plot(model_c$VCV)    	 
plot(model_c$Sol)

# 提取路径a的后验样本：TmpVaramongyears对BSNonCooMean的影响
post_samples_a <-model_a$Sol[, "TmpVaramongyears"]
# 提取路径b的后验样本： BSNonCooMean对SMD的影响
post_samples_b <- model_c$Sol[, "BSNonCooMean"]
# 计算间接效应的后验样本
indirect_effect_samples <- post_samples_a * post_samples_b
# 计算间接效应的后验总结
indirect_effect_summary <- c(mean = mean(indirect_effect_samples),
                             lower = quantile(indirect_effect_samples, 0.05),
                             upper = quantile(indirect_effect_samples, 0.95))
# 输出间接效应的后验总结
indirect_effect_summary









###路径分析试试
model <- '
  MalnRR ~ MaPairS 
  MalnRR ~ TmpVaramongyears
  MalnRR ~  SMD
  MaPairS ~ TmpVaramongyears
  SMD ~  TmpVaramongyears
  BSNonCooMean ~  TmpVaramongyears
  BSNonCooMean ~ MaPairS
  SMD ~ BSNonCooMean
  SMD ~ MaPairS
'
fit <- sem(model, data=a_new_scaled)

semPaths(fit, whatLabels="est", layout="tree", rotation=3)

fit$ParTable$lhs
fit$vnames$ov$chr [1:5]

labels <- c("温度年变异", "BSN的全称", "合作提高出飞数", "非合作雄性存活率", "合作对雄性存活率提升")

semPaths(fit, whatLabels="est", layout="tree", labels=labels)




#install.packages("lavaan")
library(lavaan)
a_new_scaled <- as.data.frame(scale(a_new[,4:36]))
fit <- sem(model, data=a_new_scaled)
summary(fit)

model <- ' 
  # 直接路径
  MaPairS ~ b1*TmpVaramongyears
  SMD ~ b2*TmpVaramongyears
  MalnRR ~ b3*MaPairS + b4*SMD + b5*TmpVaramongyears
  
  # 间接效应
  indirect1 := b1 * b3
  indirect2 := b2 * b4
'

# 拟合模型并进行bootstrap
fit <- sem(model, data=a_new_scaled, se = "bootstrap", bootstrap = 1000)

# 查看结果
summary(fit, fit.measures=TRUE)




install.packages("semPlot")
library(semPlot)
semPaths(fit, whatLabels="est", layout="circle", rotation=3)



library(semPlot)

# 从模型输出中获取估计值和p值
estimates <- parameterEstimates(fit)

# 构建包含估计值和p值的标签
labels <- paste0("Estimate: ", round(estimates$est, 2), 
                 "\np-value: ", round(estimates$pvalue, 3))

# 绘制路径图，并手动添加标签
semPaths(fit, whatLabels="est",
         edge.label.cex = 0.8,
         node.label.cex = 0.8,
         edge.label.color = "blue",
         layout = "tree",
         rotation = 2,
         esize = 2,
         nCharNodes = 0,
         nCharEdges = 0,
         asize = 5,
         edge.labels = labels) # 使用创建的标签







fixedModelEnv1 <- MCMCglmm(
  MalnRR ~ scale(MaPairS)+SMD+TmpVaramongyears+scale(MaPairS):SMD+scale(MaPairS):TmpVaramongyears+SMD:TmpVaramongyears, 
  random = ~ Species,
  data=a_new,
  pl=TRUE,
  slice=FALSE,
  nitt=103000,
  thin=10,
  burnin=3000,
  prior=prior1,
  ginverse=list(Species=INtree$Ainv),
  verbose=FALSE
)

summary(fixedModelEnv1)




summary(fixedModelEnv1)

hist(fixedModelEnv1$Liab)
plot(fixedModelEnv1$VCV)    	 
plot(fixedModelEnv1$Sol)

# 安装和加载所需的包
if (!require(car)) install.packages("car")
library(car)

# 使用普通最小二乘法拟合模型
model_ols <- lm(SMD ~ scale(MaPairS) + scale(TmpVaramongyears), data=a_new)

# 计算方差膨胀因子
vif(model_ols)











which(a_new$Species == "Erythropygia_coryphaeus")

plot( a_new$PreVaramongyears,a_new$lnRR, main="Scatterplot", xlab="PreVaramongyears", ylab="lnRR")

install.packages("Hmisc")
library(Hmisc)
# 移除非数值型列
a_new_numeric <- a_new[sapply(a_new, is.numeric)]

result <- rcorr(as.matrix(a_new_numeric))


# 显示结果
print(result$r)  # 相关系数矩阵
print(result$P)  # 相关性显著性的 p 值矩阵


# 输出相关系数矩阵到 CSV 文件
write.csv(result$r, file = "correlation_matrix.csv")

# 输出 p 值矩阵到 CSV 文件
write.csv(result$P, file = "pvalue_matrix.csv")

#install.packages("corrplot")
library(corrplot)
# 可视化相关系数矩阵
corrplot(result$r, method = "circle")
corrplot(result$P, method = "circle")

# 设置图形参数，创建一个 1 行 2 列的图形布局
par(mfrow = c(1, 2))

# 可视化相关系数矩阵
corrplot(result$r, method = "circle", title = "Correlation Matrix")

# 可视化 p 值矩阵
corrplot(result$P, method = "circle", title = "P-value Matrix", is.corr = FALSE)

# 如果您想要一个更详细的图表，可以添加颜色和标签
corrplot(result$r, method = "color", type = "upper", order = "hclust", 
         addCoef.col = "black", tl.col = "black", tl.srt = 45)




#########################################
library(MCMCglmm)
MA_dataset<-as.data.frame(read.csv("P5S1.csv"))
a<-MA_dataset
a[1:5,]
which(names(a) == "NonCooBR")
a[,10:11] <- lapply(a[,10:11], function(x) ifelse(is.na(as.numeric(x)), NA, x))
a[,10:11] <- lapply(a[,10:11], as.numeric)
na_rows <- which(is.na(a$NonCooBR))
a_new <- a[-na_rows, ]
a_new <- a
dim(a_new)
#a_new$lnBR<- log(a_new$NonCooBR/a_new$CooBR)
#a_new$CooBR_NonCooBR_diff = a_new$CooBR - a_new$NonCooBR
#write.csv(a_new,"BRlevel.csv")
trees <- read.nexus("44species1000times.nex")
dropTip <- trees[[1]]$tip.label[which(is.na(match(trees[[1]]$tip.label, a_new$Species)))]
zrTrees <- lapply(trees, drop.tip, dropTip, trim.internal=T)
zrTrees <- lapply(zrTrees, makeNodeLabel, method = "number")
a$Species[which((a_new$Species %in% zrTrees[[1]]$tip.label) == FALSE)]
zrTrees[[1]]$tip.label[which((zrTrees[[1]]$tip.label %in% a_new$Species) == FALSE)]

colnames(a)
a1 <- 1000
prior1<-list(R=list(V=diag(1),nu=0.002)
             , G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a1)))


INtree <- inverseA(zrTrees[[1]], nodes="TIPS")
fixedModelEnv1<- MCMCglmm( DOV ~ PreVaramongyears, random = ~ Species,data=a_new ,pl=TRUE, slice=FALSE, nitt=103000, thin=10, burnin=3000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
summary(fixedModelEnv1)

hist(fixedModelEnv1$Liab)
plot(fixedModelEnv1$VCV)    	 
plot(fixedModelEnv1$Sol)


library(ggplot2)
ggplot(a_new, aes(x = NonCooBR, y = PreVaramongyears)) +
  geom_point(color = "blue") +  # 绘制蓝色的点
  labs(title = "Scatter Plot of SMD vs NonCooBR",
       x = "NonCooBR",
       y = "PreVaramongyears")+
  theme_minimal()  # 使用简洁的主题












