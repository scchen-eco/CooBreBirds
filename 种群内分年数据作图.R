library(meta)
library(metafor)
library(MCMCglmm)
setwd("~")
a<-read.csv("P5S2.csv")
trees <- read.nexus("44species1000times.nex")
dropTip <- trees[[1]]$tip.label[which(is.na(match(trees[[1]]$tip.label, a$Species)))]
zrTrees <- lapply(trees, drop.tip, dropTip, trim.internal=T)
zrTrees <- lapply(zrTrees, makeNodeLabel, method = "number")
a$Species[which((a$Species %in% zrTrees[[1]]$tip.label) == FALSE)]
zrTrees[[1]]$tip.label[which((zrTrees[[1]]$tip.label %in% a$Species) == FALSE)]
head(a)
rdat <- escalc(measure = "SMDH", 
               m1i = a$BSCooMean, n1i = a$BSCoon, sd1i = a$BSCooSD, 
               m2i = a$BSNonCooMean, n2i = a$BSNonCoon, sd2i = a$BSNonCooSD)
srdat <- summary(rdat, digits = 2)
a<- cbind(a, srdat)
write.csv(a,"a_SMD.csv")
a1 <- 1000
prior1 <- list(R = list(V = diag(1), nu = 0.001),
               G = list(
                 G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1) * a1),  # For Species
                 G2 = list(V = diag(2), nu = 1, alpha.mu = c(0, 0), alpha.V = diag(2) * a1)  # For us(1+BSNonCooMean):Studies
               ))
MEV <- a$vi
INtree <- inverseA(zrTrees[[1]], nodes="TIPS")


#增强迭代
interceptModel.start <- MCMCglmm( yi ~ BSNonCooMean, random = ~ Species+us(1+ BSNonCooMean):Studies,data=a, mev=MEV,pl=TRUE, slice=FALSE, nitt=10030000, thin=100, burnin=30000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
summary(interceptModel.start)
hist(interceptModel.start$Liab)
plot(interceptModel.start$VCV)     		# species and units close to 0
plot(interceptModel.start$Sol)		   	# intercept estimate well mixed


##不控制MEV
interceptModel.start <- MCMCglmm( yi ~ BSNonCooMean, random = ~ Species+us(1+ BSNonCooMean):Studies,data=a,pl=TRUE, slice=FALSE, nitt=1030000, thin=100, burnin=30000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)
summary(interceptModel.start)



#模型重构：
prior1 <- list(
  R = list(V = diag(1), nu = 0.002),
  G = list(
    G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1) * a1),
    G2 = list(V = diag(2), nu = 1, alpha.mu = c(0, 0), alpha.V = diag(2) * a1),  # For us(1+BSNonCooMean):Studies
    G3 = list(V = 1, fix=1)
  )
)
a$MEV<-a$vi
fixedModel.start <- MCMCglmm(yi ~ BSNonCooMean, random = ~ Species+us(1+ BSNonCooMean):Studies+idh(sqrt(MEV)):units,data=a, pl=TRUE, slice=FALSE, nitt=1030000, thin=100, burnin=30000, prior=prior1, ginverse=list(Species=INtree$Ainv), verbose=FALSE)  ##使用mev进行权重
Modela.predict <- fixedModel.start
summary(Modela.predict)

a_new<-a
p1<-predict(Modela.predict, newdata=a_new, marginal= ~ Species+us(1+ BSNonCooMean):Studies+idh(sqrt(MEV)):units,
            type="response", interval="prediction", level=0.95)
a <- cbind(a, p1)

#pr1<-
ggplot(a, aes(x = BSNonCooMean, y = yi, color = Species)) +
  geom_point(shape = 19, size = 1) +  # 添加散点
  geom_line(aes(y = fit),linewidth = 1.1) +  # 添加拟合线
  geom_ribbon(aes(ymin = lwr, ymax = upr,fill = Species),color=NA,alpha = 0.3) +  # 添加线的区间
 
  
   #scale_color_manual(values = c("#376092", "red"), guide = "none") +
  #scale_fill_manual(values = c("blue", "red"), guide = "none") +
  #labs(x = "BSNonCooMean", y = "yi", title = "Scatter Plot with Linear Fit") +  # 添加标签和标题
  theme_classic()
  
scale_x_continuous(breaks = seq(0, 10, 2),limits = c(0, 10),expand = expansion(mult = c(0.05, 0.06))) +
  scale_y_continuous(breaks = seq(0, 10, 2), limits = c(0, 10),expand = expansion(mult = c(0.05, 0.06)))+
  theme(
    axis.line = element_line(linewidth =0.6),  # 调整坐标轴线条粗细
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.6, color = "black"),  # 调整x轴刻度标线的长度和粗细
    axis.ticks.y = element_line(linewidth = 0.6, color = "black") ,  # 调整y轴刻度标线的长度和粗细
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x=element_text(vjust=1,size=15,face = "bold"),
    axis.text.y=element_text(vjust=1,size=15,face = "bold")
  )
pr1
pr1 + 
  #geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#5C255C", size = 0.8)+
  labs(x = NULL, y = NULL, title = NULL)  # 删除坐标轴和标题



scale_x_continuous(breaks = seq(0, 7.5, 0.5), limits = c(0, 7.5)) +
  scale_y_continuous(breaks = seq(0, 9, 0.5), limits = c(0, 9))
head(a)
dim(a)
a<-read.csv("a_SMD.csv")
###############
ggplot(a, aes(x = BSNonCooMean, y = yi, color = Species)) + 
  geom_point() +  # 添加点
  geom_smooth(method = "lm", se = FALSE) +  # 添加线性回归线，不显示置信区间
  #geom_abline(slope = -0.43519, intercept = 1.56469, color = "blue", linetype = "dashed") +  # 添加固定斜率和截距的直线
  geom_segment(aes(x = 0, y = 1.56469, xend =6, yend = -1.04645), color = "black", size = 1) +  # 添加黑色实线，指定起点和终点
  theme_classic()  # 应用经典主题

x <- 6
y <- -0.43519 * x + 1.56469
print(y)
###########
# 计算斜率范围对应的y值范围
a$y_min <- -0.90979 * a$BSNonCooMean + 4.12266
a$y_max <- 0.04137 * a$BSNonCooMean -0.56571

ggplot(a, aes(x = BSNonCooMean, y = yi, color = Species)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  #geom_abline(slope = -0.90979, intercept = 4.12266, color = "black") + 
  geom_abline(slope = -0.43519, intercept = 1.56469, color = "black",size=1) + 
  #geom_abline(slope = 0.04137, intercept = -0.56571, color = "black") + 
  theme_classic()+
  theme(legend.text = element_text(size = 12, face = "italic"))





