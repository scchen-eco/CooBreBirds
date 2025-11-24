library(MCMCglmm)
library(dplyr)
#招募数和成体存活
a<- read.csv("dx_regression_Pt.csv")
a<- a[a$state != "H", ]
head(a)
a$Recruits1 <- a$Recruits * 2
a$Experience <- ifelse(is.na(a$Age), NA, ifelse(a$Age > 1, 1, 0))
a$Experience <- replace(a$Experience, !is.na(a$Experience), scale(a$Experience[!is.na(a$Experience)]))
a$Experience[is.na(a$Experience)] <- median(a$Experience, na.rm = TRUE)
a$Ter.z <- scale(a$Ter.)
#绘制合作与非合作各自的ksai随降水变化的图
precip_data <- read.csv("cleaned_dataday_pre.csv")
temp_data <- read.csv("cleaned_climlongday.csv")
# 将日期列转换为正确的日期格式 (dd/mm/yyyy)
precip_data$Date <- as.Date(precip_data$Date, format = "%d/%m/%Y")
temp_data$Date <- as.Date(temp_data$Date, format = "%d/%m/%Y")
#合并环境数据
temp_nov_dec <- temp_data %>%
  filter(format(Date, "%m") %in% c("11")) %>%  # 筛选11月
  mutate(year = format(Date, "%Y")) %>%  # 提取年份
  group_by(year) %>%
  summarise(avg_temp_nov_dec = mean(tmpvalue, na.rm = TRUE))  # 计算平均温度
precip_apr_may <- precip_data %>%
  filter(format(Date, "%m") %in% c("04", "05")) %>%  # 筛选4月和5月
  mutate(year = format(Date, "%Y")) %>%  # 提取年份
  group_by(year) %>%
  summarise(total_precip_apr_may = sum(Pre, na.rm = TRUE))  # 计算总降水量
a$year <- as.character(a$year)  # 确保年份格式一致
# 合并温度数据和降水数据到a中
a_merged <- a %>%
  left_join(temp_nov_dec, by = "year") %>%  # 合并平均温度
  left_join(precip_apr_may, by = "year")  # 合并总降水量
a_merged$TempNov.z <- scale(a_merged$avg_temp_nov_dec)
a_merged$PreAM.z <- scale(a_merged$total_precip_apr_may)
#设置 prior（允许 residual covariance）
prior <- list(
  R = list(V = diag(2), nu = 0.002),  # residual covariance
  G = list(
    G1 = list(V = diag(2), nu = 0.002),  # for year
    G2 = list(V = diag(2), nu = 0.002)   # for id
  )
)

model_SurRecr_trade<- MCMCglmm(
  cbind(Survival, Recruits1) ~ trait * (sex + Experience + Ter.z + Gethp ),
  rcov = ~ us(trait):units,
  random = ~ us(trait):PreAM.z + us(trait):id,
  family = c("categorical", "poisson"),
  data = a_merged,
  prior = prior,
  nitt = 720000, burnin = 144000, thin = 30,
  verbose = FALSE
)
summary(model_SurRecr_trade)

model_SurRecr_trade1<- MCMCglmm(
  cbind(Survival, Recruits1) ~ trait * (sex + Experience + Ter.z + Gethp + TempNov.z +PreAM.z),
  rcov = ~ us(trait):units,
  random = ~ us(trait):PreAM.z:Gethp+ us(trait):id,
  family = c("categorical", "poisson"),
  data = a_merged,
  prior = prior,
  nitt = 7200, burnin = 1440, thin = 30,
  verbose = FALSE
)
summary(model_SurRecr_trade1)







