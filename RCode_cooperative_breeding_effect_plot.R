library(dplyr)
library(tidyr)
library(metafor)
a <- read.csv("NESTlEVEL.csv")
head(a)
a$Gethp <- as.factor(a$Gethp)  # Coo 表示是否合作，1 表示合作，0 表示非合作
a$brood <- as.numeric(as.character(a$brood))  # brood 表示繁殖力
# 分组计算每年合作和非合作巢的统计值
annual_stats <- a %>%
  group_by(year, Gethp) %>%  # 按年份和是否合作分组
  summarise(
    mean_RS = mean(brood, na.rm = TRUE),  
    sd_RS = sd(brood, na.rm = TRUE),      
    n = n()                            
  )
#print(annual_stats)
#将长数据整理为宽数据
# 首先确保 Gethp 列是数值型
annual_stats <- annual_stats %>%
  mutate(Gethp = as.numeric(as.character(Gethp)))  # 将 Gethp 转换为数值

# 将长数据转换为宽数据
wide_stats <- annual_stats %>%
  pivot_wider(
    names_from = Gethp,  # 按照 Gethp（合作状态）展开列
    values_from = c(mean_RS, sd_RS, n),  # 使用 mean_RS, sd_RS 和 n 创建新列
    names_prefix = "Coo_",  # 给新列加上后缀
    names_sep = "_"
  )  %>%
  rename(
    Coo_mean_RS = `mean_RS_Coo_1`,       # 重命名合作巢的列
    Non_Coo_mean_RS = `mean_RS_Coo_0`,   # 重命名非合作巢的列
    Coo_sd_RS = `sd_RS_Coo_1`,           # 重命名合作巢的标准差列
    Non_Coo_sd_RS = `sd_RS_Coo_0`,       # 重命名非合作巢的标准差列
    Coo_n = `n_Coo_1`,                   # 重命名合作巢的样本量列
    Non_Coo_n = `n_Coo_0`                # 重命名非合作巢的样本量列
  )

# 使用 escalc 函数计算标准化均差 (SMDH)
rdat <- escalc(measure = "SMDH", 
               m1i = wide_stats$Coo_mean_RS,    # 合作巢的均值
               n1i = wide_stats$Coo_n,          # 合作巢的样本量
               sd1i = wide_stats$Coo_sd_RS,     # 合作巢的标准差
               m2i = wide_stats$Non_Coo_mean_RS, # 非合作巢的均值
               n2i = wide_stats$Non_Coo_n,      # 非合作巢的样本量
               sd2i = wide_stats$Non_Coo_sd_RS) # 非合作巢的标准差

# 输出结果的摘要
srdat <- summary(rdat, digits = 2)


a_broodsize<- cbind(wide_stats, srdat)
# a_broodsize<- as.data.frame(a_broodsize)

Model_DX <- rma.mv(yi ~ 1, V = vi, random = list(~ 1 | year), data = a_broodsize)
summary(Model_DX)
# 提取效应均值（合作繁殖效应的 SMD）
mean_effect <- coef(Model_DX)[1]  

# 提取随机截距的方差（即年间方差）
var_component <- Model_DX$sigma2[1]

# 计算标准差
sd_effect <- sqrt(var_component)

# 计算变异系数（CV），乘以 100 变为百分比形式
cv <- (sd_effect / mean_effect) * 100
cv



#########################################################################
### colors to be used in the plot
colp <- "#6b58a6"
coll <- "#a7a9ac"
k <- nrow(a_broodsize)
### generate point sizes
psize <- weights(Model_DX)
psize <- 1.2 + (psize - min(psize)) / (max(psize) - min(psize))

### get the weights and format them as will be used in the forest plot
weightsDX <- fmtx(weights(Model_DX), digits=1)

### adjust the margins
#par(mar=c(2.7,3.2,2.5,1.3), mgp=c(3,0,0), tcl=0.15) 
#下左上右；坐标轴标题距离坐标轴的距离，坐标轴标签距离坐标轴的距离，坐标轴刻度线的距离；设置坐标轴刻度线的长度
### forest plot with extra annotations
forest_DX_broodsize<- forest(a_broodsize$yi, a_broodsize$vi,slab = a_broodsize$year,   xlim=c(-4,4), ylim=c(-1,k+3), alim=c(-2,2), cex=0.88,
                             header=FALSE, pch=18, psize=psize, efac=0, refline=NA, lty=c(1,0), xlab="",
                             rowadj=-.07,
                             annotate = FALSE # 去除右侧注释
)
### add the vertical reference line at 0
segments(0, -1, 0, k+1.6, col=coll)
### add the vertical reference line at the pooled estimate
segments(coef(Model_DX), 0, coef(Model_DX), k+1, col=colp, lty="33", lwd=0.8)
### redraw the CI lines and points in the chosen color
segments(a_broodsize$ci.lb, k:1, a_broodsize$ci.ub, k:1, col=colp, lwd=1.5)
points(a_broodsize$yi, k:1, pch=18, cex=psize*1.15, col="white")
points(a_broodsize$yi, k:1, pch=18, cex=psize, col=colp)
### add the summary polygon
addpoly(Model_DX, row=0, mlab="Total (95% CI)", efac=2, col=colp, border=colp)
### add the horizontal line at the top
abline(h=k+1.6, col=coll)
### redraw the x-axis in the chosen color
#axis(side=1, at=seq(-3,3,by=0.5), col=coll, labels=FALSE)
plot1 <- recordPlot()
# 随时重新调用
replayPlot(plot1)

#######################################################################幼体存活
a<- read.csv("DX_offss.csv")
# 计算每年合作与非合作存活率均值和样本量
survival_stats <- a %>%
  group_by(Year, Gethp) %>%
  summarise(
    mean_survival = mean(SUR, na.rm = TRUE),  # 存活率均值
    sample_size = n()                          # 样本量
  ) %>%
  ungroup()
print(survival_stats)
# 提取合作和非合作的存活率和样本量
cooperative <- survival_stats %>% filter(Gethp == 1)
non_cooperative <- survival_stats %>% filter(Gethp == 0)

lnRR_data <- cooperative %>%
  select(Year, mean_survival, sample_size) %>%
  rename(Coo_mean_survival = mean_survival, Coo_n = sample_size) %>%
  left_join(non_cooperative %>% select(Year, mean_survival, sample_size) %>%
              rename(Non_Coo_mean_survival = mean_survival, Non_Coo_n = sample_size),
            by = "Year")
# 函数定义
Calc.lnRR <- function(P1, P2){
  log(P1) - log(P2)
}

Calc.var.lnRR <- function(P1, N1, P2, N2){
  ((1-P1)/(N1*P1)) + ((1-P2)/(N2*P2))
}
# 计算 lnRR 和其方差
lnRR_data <- lnRR_data %>%
  mutate(
    # 使用你提供的公式计算 lnRR
    lnRR = Calc.lnRR(Coo_mean_survival, Non_Coo_mean_survival),  
    # 计算方差
    var_lnnRR = Calc.var.lnRR(Coo_mean_survival, Coo_n, Non_Coo_mean_survival, Non_Coo_n)
  )

# 查看 lnRR 结果
print(lnRR_data)

# 剔除lnRR为-Inf的行
lnRR_data_clean <- lnRR_data %>% filter(lnRR != -Inf)
lnRR_data_clean

# 在这里构建模型
Model_DX <- rma.mv(lnRR ~ 1, V = var_lnnRR, random = list(~ 1 | Year), data = lnRR_data_clean)

# 打印模型摘要
summary(Model_DX)

lnRR_CV_summary <- lnRR_data_clean %>%
  summarise(
    mean_lnRR = mean(lnRR, na.rm = TRUE),
    sd_lnRR = sd(lnRR, na.rm = TRUE),
    CV = sd_lnRR / mean_lnRR * 100
  )

print(lnRR_CV_summary)

### colors to be used in the plot
colp <- "#6b58a6"
coll <- "#a7a9ac"
k <- nrow(lnRR_data_clean)
### generate point sizes
psize <- weights(Model_DX)
psize <- 1.2 + (psize - min(psize)) / (max(psize) - min(psize))

### get the weights and format them as will be used in the forest plot
weightsDX <- fmtx(weights(Model_DX), digits=1)

### adjust the margins
#par(mar=c(2.7,3.2,2.5,1.3), mgp=c(3,0,0), tcl=0.15) 
#下左上右；坐标轴标题距离坐标轴的距离，坐标轴标签距离坐标轴的距离，坐标轴刻度线的距离；设置坐标轴刻度线的长度


### forest plot with extra annotations
forest_DX_offs<- forest(lnRR_data_clean$lnRR, lnRR_data_clean$var_lnnRR,slab = lnRR_data_clean$Year,   xlim=c(-4,4), ylim=c(-1,k+3), alim=c(-3,3), cex=0.88,
                        header=FALSE, pch=18, psize=psize, efac=0, refline=NA, lty=c(1,0), xlab="",
                        rowadj=-.07,
                        annotate = FALSE # 去除右侧注释
)

### add the vertical reference line at 0
segments(0, -1, 0, k+1.6, col=coll)
### add the vertical reference line at the pooled estimate
segments(coef(Model_DX), 0, coef(Model_DX), k+1, col=colp, lty="33", lwd=0.8)
### redraw the CI lines and points in the chosen color
points(lnRR_data_clean$lnRR, k:1, pch=18, cex=psize*1.15, col="white")
points(lnRR_data_clean$lnRR, k:1, pch=18, cex=psize, col=colp)
### add the summary polygon
addpoly(Model_DX, row=0, mlab="Total (95% CI)", efac=2, col=colp, border=colp)
### add the horizontal line at the top
abline(h=k+1.6, col=coll)

### redraw the x-axis in the chosen color
#axis(side=1, at=seq(-3,3,by=0.5), col=coll, labels=FALSE)
# 记录绘图
plot2 <- recordPlot()
# 随时重新调用
replayPlot(plot2)

#######################################################################成体存活
# 读取数据
a <- read.csv("P_SUR_Age_TQ_Ex.csv")
# 计算每年合作与非合作存活率均值和样本量
survival_stats <- a %>%
  group_by(Year, Recieved) %>%
  summarise(
    mean_survival = mean(SUR, na.rm = TRUE),  # 存活率均值
    sample_size = n()                          # 样本量
  ) %>%
  ungroup()
# 查看计算结果
print(survival_stats)
# 提取合作和非合作的存活率和样本量
cooperative <- survival_stats %>% filter(Recieved == 1)
non_cooperative <- survival_stats %>% filter(Recieved == 0)
# 确保两个数据框有相同的年份
lnRR_data <- cooperative %>%
  select(Year, mean_survival, sample_size) %>%
  rename(Coo_mean_survival = mean_survival, Coo_n = sample_size) %>%
  left_join(non_cooperative %>% select(Year, mean_survival, sample_size) %>%
              rename(Non_Coo_mean_survival = mean_survival, Non_Coo_n = sample_size),
            by = "Year")
# 函数定义
Calc.lnRR <- function(P1, P2){
  log(P1) - log(P2)
}

Calc.var.lnRR <- function(P1, N1, P2, N2){
  ((1-P1)/(N1*P1)) + ((1-P2)/(N2*P2))
}
# 计算 lnRR 和其方差
lnRR_data <- lnRR_data %>%
  mutate(
    # 使用你提供的公式计算 lnRR
    lnRR = Calc.lnRR(Coo_mean_survival, Non_Coo_mean_survival),  
    # 计算方差
    var_lnnRR = Calc.var.lnRR(Coo_mean_survival, Coo_n, Non_Coo_mean_survival, Non_Coo_n)
  )

# 在这里构建模型
Model_DX <- rma.mv(lnRR ~ 1, V = var_lnnRR, random = list(~ 1 | Year), data = lnRR_data)
# 打印模型摘要
summary(Model_DX)
mean_effect <- coef(Model_DX)[1]  
var_component <- Model_DX$sigma2[1]
sd_effect <- sqrt(var_component)
cv <- (sd_effect / mean_effect) * 100
cv

### colors to be used in the plot
colp <- "#6b58a6"
coll <- "#a7a9ac"
k <- nrow(lnRR_data)
### generate point sizes
psize <- weights(Model_DX)
psize <- 1.2 + (psize - min(psize)) / (max(psize) - min(psize))
### get the weights and format them as will be used in the forest plot
weightsDX <- fmtx(weights(Model_DX), digits=1)
### adjust the margins
##par(mar=c(2.7,3.2,2.5,1.3), mgp=c(3,0,0), tcl=0.15) 
#下左上右；坐标轴标题距离坐标轴的距离，坐标轴标签距离坐标轴的距离，坐标轴刻度线的距离；设置坐标轴刻度线的长度
### forest plot with extra annotations
forest_DX_adult<- forest(lnRR_data$lnRR, lnRR_data$var_lnnRR,slab = lnRR_data$Year,   xlim=c(-4,4), ylim=c(-1,k+3), alim=c(-1.5,1.5), cex=0.88,
                         header=FALSE, pch=18, psize=psize, efac=0, refline=NA, lty=c(1,0), xlab="",
                         rowadj=-.07,
                         annotate = FALSE # 去除右侧注释
)
### add the vertical reference line at 0
segments(0, -1, 0, k+1.6, col=coll)

### add the vertical reference line at the pooled estimate
segments(coef(Model_DX), 0, coef(Model_DX), k+1, col=colp, lty="33", lwd=0.8)

### redraw the CI lines and points in the chosen color
points(lnRR_data$lnRR, k:1, pch=18, cex=psize*1.15, col="white")
points(lnRR_data$lnRR, k:1, pch=18, cex=psize, col=colp)

### add the summary polygon
addpoly(Model_DX, row=0, mlab="Total (95% CI)", efac=2, col=colp, border=colp)

### add the horizontal line at the top
abline(h=k+1.6, col=coll)

### redraw the x-axis in the chosen color
# axis(side=1, at=seq(-1.5,1.5,by=0.5), col=coll, labels=FALSE)
# 记录绘图
plot3 <- recordPlot()
# 随时重新调用
replayPlot(plot3)

##################################################################### 对Ksai的合作繁殖效应
a<- read.csv("dx_regression_Pt.csv")
head(a)
# 分组计算每年合作和非合作巢的统计值
annual_stats <- a %>%
  filter(state != "H")%>%
  group_by(year, Gethp) %>%  # 按年份和是否合作分组
  summarise(
    mean_ksai = mean(ksai, na.rm = TRUE),  
    sd_ksai = sd(ksai, na.rm = TRUE),      
    n = n()                            
  )
print(annual_stats)
# 将长数据转换为宽数据
wide_stats <- annual_stats %>%
  pivot_wider(
    names_from = Gethp,  # 按照 Gethp（合作状态）展开列
    values_from = c(mean_ksai, sd_ksai, n),  # 使用 mean_ksai, sd_ksai 和 n 创建新列
    names_prefix = "Coo_",  # 给新列加上后缀
    names_sep = "_"
  )  %>%
  rename(
    Coo_mean_ksai = `mean_ksai_Coo_1`,       # 重命名合作巢的列
    Non_Coo_mean_ksai = `mean_ksai_Coo_0`,   # 重命名非合作巢的列
    Coo_sd_ksai = `sd_ksai_Coo_1`,           # 重命名合作巢的标准差列
    Non_Coo_sd_ksai = `sd_ksai_Coo_0`,       # 重命名非合作巢的标准差列
    Coo_n = `n_Coo_1`,                   # 重命名合作巢的样本量列
    Non_Coo_n = `n_Coo_0`                # 重命名非合作巢的样本量列
  )

# 使用 escalc 函数计算标准化均差 (SMDH)
rdat <- escalc(measure = "SMDH", 
               m1i = wide_stats$Coo_mean_ksai,    # 合作巢的均值
               n1i = wide_stats$Coo_n,          # 合作巢的样本量
               sd1i = wide_stats$Coo_sd_ksai,     # 合作巢的标准差
               m2i = wide_stats$Non_Coo_mean_ksai, # 非合作巢的均值
               n2i = wide_stats$Non_Coo_n,      # 非合作巢的样本量
               sd2i = wide_stats$Non_Coo_sd_ksai) # 非合作巢的标准差

# 输出结果的摘要
srdat <- summary(rdat, digits = 2)
a_ksai<- cbind(wide_stats, srdat)
Model_DX <- rma.mv(yi ~ 1, V = vi, random = list(~ 1 | year), data = a_ksai)
summary(Model_DX)
mean_effect <- coef(Model_DX)[1]  
var_component <- Model_DX$sigma2[1]
sd_effect <- sqrt(var_component)
cv <- (sd_effect / mean_effect) * 100
cv


#########################################################################
### colors to be used in the plot
colp <- "#6b58a6"
coll <- "#a7a9ac"
k <- nrow(a_ksai)
### generate point sizes
psize <- weights(Model_DX)
psize <- 1.2 + (psize - min(psize)) / (max(psize) - min(psize))

### get the weights and format them as will be used in the forest plot
weightsDX <- fmtx(weights(Model_DX), digits=1)

### adjust the margins
#par(mar=c(2.7,3.2,2.5,1.3), mgp=c(3,0,0), tcl=0.15) 
#下左上右；坐标轴标题距离坐标轴的距离，坐标轴标签距离坐标轴的距离，坐标轴刻度线的距离；设置坐标轴刻度线的长度
### forest plot with extra annotations
forest_DX_ksaisize<- forest(a_ksai$yi, a_ksai$vi,slab = a_ksai$year,   xlim=c(-4,4), ylim=c(-1,k+3), alim=c(-2,2), cex=0.88,
                            header=FALSE, pch=18, psize=psize, efac=0, refline=NA, lty=c(1,0), xlab="",
                            rowadj=-.07,
                            annotate = FALSE # 去除右侧注释
)
### add the vertical reference line at 0
segments(0, -1, 0, k+1.6, col=coll)
### add the vertical reference line at the pooled estimate
segments(coef(Model_DX), 0, coef(Model_DX), k+1, col=colp, lty="33", lwd=0.8)
### redraw the CI lines and points in the chosen color
segments(a_ksai$ci.lb, k:1, a_ksai$ci.ub, k:1, col=colp, lwd=1.5)
points(a_ksai$yi, k:1, pch=18, cex=psize*1.15, col="white")
points(a_ksai$yi, k:1, pch=18, cex=psize, col=colp)
### add the summary polygon
addpoly(Model_DX, row=0, mlab="Total (95% CI)", efac=2, col=colp, border=colp)
### add the horizontal line at the top
abline(h=k+1.6, col=coll)
### redraw the x-axis in the chosen color
#axis(side=1, at=seq(-3,3,by=0.5), col=coll, labels=FALSE)
plot4 <- recordPlot()
# 随时重新调用
replayPlot(plot4)



#####################################################################排版
# 设置绘图区域为1行3列
par(mfrow = c(1, 4), mar = c(5, 4, 4, 2) + 0.1)
#分批运行上面代码
plot1
plot2
plot3
plot4