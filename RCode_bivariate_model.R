########################################双变量模型构建##################################################
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
length(a_merged$id)
length(unique(a_merged$fn))
head(a_merged)
table(a_merged$sex)
unique_individuals <- a_merged[!duplicated(a_merged$id), c("id", "sex")]

# 2. 使用 table() 函数统计 sex 列的分布
table(unique_individuals$sex)

#设置 prior（允许 residual covariance）
prior <- list(
  R = list(V = diag(2), nu = 0.002),  # residual covariance
  G = list(
    G1 = list(V = diag(2), nu = 0.002),  # for year
    G2 = list(V = diag(2), nu = 0.002)   # for id
  )
)

model_SurRecr<- MCMCglmm(
  cbind(Survival, Recruits1) ~ trait * (sex + Experience + Ter.z + Gethp + TempNov.z + PreAM.z +
                                          Gethp:TempNov.z + Gethp:PreAM.z +
                                          TempNov.z:PreAM.z + Gethp:TempNov.z:PreAM.z) - 1,
  rcov = ~ us(trait):units,
  random = ~ us(trait):year + us(trait):id,
  family = c("categorical", "poisson"),
  data = a_merged,
  prior = prior,
  nitt = 720000, burnin = 144000, thin = 30,
  verbose = FALSE
)
summary(model_SurRecr)
saveRDS(model_SurRecr, file = "model_SurRecr.rds")

data_helped <- a_merged %>% filter(Gethp == 1)
data_nonhelped <- a_merged %>% filter(Gethp == 0)
nrow(data_helped)
nrow(data_nonhelped)
model_SurRecr_helped <- MCMCglmm(
  cbind(Survival, Recruits1) ~ trait * (sex + Experience + Ter.z + TempNov.z + PreAM.z +
                                          TempNov.z:PreAM.z) - 1,
  rcov = ~ us(trait):units,
  random = ~ us(trait):year + us(trait):id,
  family = c("categorical", "poisson"),
  data = data_helped,
  prior = prior,
  nitt = 720000, burnin = 144000, thin = 30,
  verbose = FALSE
)
summary(model_SurRecr_helped )
saveRDS(model_SurRecr_helped, file = "model_SurRecr_helped.rds")

model_SurRecr_nonhelped <- MCMCglmm(
  cbind(Survival, Recruits1) ~ trait * (sex + Experience + Ter.z + TempNov.z + PreAM.z +
                                          TempNov.z:PreAM.z) - 1,
  rcov = ~ us(trait):units,
  random = ~ us(trait):year + us(trait):id,
  family = c("categorical", "poisson"),
  data = data_nonhelped,
  prior = prior,
  nitt = 720000, burnin = 144000, thin = 30,
  verbose = FALSE
)
summary(model_SurRecr_nonhelped)
saveRDS(model_SurRecr_nonhelped, file = "model_SurRecr_nonhelped.rds")


#出飞窝雏数和幼体存活的双变量模型——个体水平的分析
#数据读取
Nest_Fle<- read.csv("NESTlEVEL.csv")
offspring_Sur<-read.csv("DX_offss.csv")  
precip_data <- read.csv("cleaned_dataday_pre.csv")
temp_data <- read.csv("cleaned_climlongday.csv")

#数据检查
head(Nest_Fle)
head(offspring_Sur)
fn_vals <- unique(offspring_Sur$fn)
missing_fn <- fn_vals[!fn_vals %in% Nest_Fle$ID]
missing_fn
id_vals <- unique(Nest_Fle$ID)
missing_id <- id_vals[!id_vals %in% offspring_Sur$fn]
missing_id

#数据处理
a_OffSurFle_clean <- offspring_Sur %>%
  left_join(Nest_Fle %>% select(ID, brood), by = c("fn" = "ID")) %>%
  # 合并11月温度
  left_join(
    temp_data %>%
      mutate(Date = as.Date(Date, format = "%d/%m/%Y")) %>%
      filter(format(Date, "%m") == "11") %>%
      mutate(year = as.integer(format(Date, "%Y"))) %>%
      group_by(year) %>%
      summarise(avg_temp_nov = mean(tmpvalue, na.rm = TRUE)),
    by = c("Year" = "year")
  ) %>%
  # 合并4-5月降水
  left_join(
    precip_data %>%
      mutate(Date = as.Date(Date, format = "%d/%m/%Y")) %>%
      filter(format(Date, "%m") %in% c("04", "05")) %>%
      mutate(year = as.integer(format(Date, "%Y"))) %>%
      group_by(year) %>%
      summarise(total_precip_apr_may = sum(Pre, na.rm = TRUE)),
    by = c("Year" = "year")
  ) %>%
  mutate(
    Exp. = case_when(
      !is.na(Exp.) & Exp. >= 1 ~ 1,
      !is.na(Exp.) & Exp. == 0 ~ 0,
      TRUE ~ Exp.
    ),
    Exp. = ifelse(is.na(Exp.), mean(Exp., na.rm = TRUE), Exp.),
    Exp. = scale(Exp.)[, 1]
  ) %>%
  filter(!is.na(SUR))%>%
  filter(Sex != "x") %>%
  mutate(
    Terr. = scale(Terr.)[, 1],
    avg_temp_nov = scale(avg_temp_nov)[, 1],
    total_precip_apr_may = scale(total_precip_apr_may)[, 1]
  )
head(a_OffSurFle_clean)
# 设定先验，允许残差协方差矩阵
prior <- list(
  R = list(V = diag(2), nu = 0.002),  # 残差协方差矩阵
  G = list(
    G1 = list(V = diag(2), nu = 0.002)  # 巢ID随机效应
  )
)

model_OffSurFle <- MCMCglmm(
  cbind(SUR, brood) ~ trait * (Gethp + Exp. + Terr. + Sex + avg_temp_nov + total_precip_apr_may +
                                 avg_temp_nov:total_precip_apr_may +
                                 Gethp:avg_temp_nov + Gethp:total_precip_apr_may +
                                 Gethp:avg_temp_nov:total_precip_apr_may),
  random = ~ us(trait):fn ,  # 添加年份随机效应
  rcov = ~ us(trait):units,
  family = c("categorical", "poisson"),
  data = a_OffSurFle_clean,
  prior = prior,
  nitt = 720000, burnin = 144000, thin = 30,
  verbose = FALSE
)
saveRDS(model_OffSurFle, file = "model_OffSurFle.rds")
summary(model_OffSurFle)


# Gethp = 1 的个体
data_helped <- a_OffSurFle_clean %>% filter(Gethp == 1)

# Gethp = 0 的个体
data_nonhelped <- a_OffSurFle_clean %>% filter(Gethp == 0)

model_OffSurFle_helped <- MCMCglmm(
  cbind(SUR, brood) ~ trait * (Exp. + Terr. + Sex + avg_temp_nov + total_precip_apr_may +
                                 avg_temp_nov:total_precip_apr_may),
  random = ~ us(trait):fn,
  rcov = ~ us(trait):units,
  family = c("categorical", "poisson"),
  data = data_helped,
  prior = prior,
  nitt = 720000, burnin = 144000, thin = 30,
  verbose = FALSE
)
summary(model_OffSurFle_helped )
saveRDS(model_OffSurFle_helped, file = "model_OffSurFle_helped.rds")


model_OffSurFle_nonhelped <- MCMCglmm(
  cbind(SUR, brood) ~ trait * (Exp. + Terr. + Sex + avg_temp_nov + total_precip_apr_may +
                                 avg_temp_nov:total_precip_apr_may),
  random = ~ us(trait):fn,
  rcov = ~ us(trait):units,
  family = c("categorical", "poisson"),
  data = data_nonhelped,
  prior = prior,
  nitt = 720000, burnin = 144000, thin = 30,
  verbose = FALSE
)
summary(model_OffSurFle_nonhelped)
saveRDS(model_OffSurFle_nonhelped, file = "model_OffSurFle_nonhelped.rds")

##############双变量模型中第二特征值效应的提取########
# 设置文件夹路径（替换为你的实际路径）
model_folder <- "F:/MCMC_bivariate_model"
# 读取所有 .rds 文件的文件名（也可以根据实际扩展名改成 .RData 或其他）
model_files <- list.files(path = model_folder, pattern = "\\.rds$", full.names = TRUE)
# 显示文件名
print(model_files)
# 创建一个空列表用于存储模型
models <- list()
# 逐个文件读取模型并存入列表，列表名为去掉路径和扩展名的文件名
for (file in model_files) {
  model_name <- tools::file_path_sans_ext(basename(file))
  models[[model_name]] <- readRDS(file)
}
# 查看模型列表中有哪些元素（模型名）
names(models)

##对model_SurRecr中第二特征值的提取
# 提取后验样本
post <- model_SurRecr$Sol
colnames(post)
var_list <- colnames(post)[!grepl("^traitRecruits1:", colnames(post))]
var_list <- var_list[-c(1,2)]  
print(var_list)
# 创建一个函数来提取每个变量对第二性状的总效应
# 总效应 = 主效应 + traitRecruits1:主效应
get_total_effect <- function(varname) {
  main <- post[, varname]
  interaction <- post[, paste0("traitRecruits1:", varname)]
  total <- main + interaction
  return(total)
}
# 存储结果
results <- data.frame(
  Variable = character(),
  Mean = numeric(),
  Lower = numeric(),
  Upper = numeric(),
  pMCMC = numeric(),
  stringsAsFactors = FALSE
)

# 逐个变量处理
for (var in var_list) {
  total_effect <- get_total_effect(var)
  hpd <- HPDinterval(as.mcmc(total_effect))
  pval <- 2 * min(mean(total_effect > 0), mean(total_effect < 0))
  
  results <- rbind(results, data.frame(
    Variable = var,
    Mean = mean(total_effect),
    Lower = hpd[1],
    Upper = hpd[2],
    pMCMC = pval
  ))
}
# 查看结果
print(results)

#model_SurRecr_helped
model_SurRecr_helped<-readRDS( "model_SurRecr_helped.rds") 
# 提取后验样本
post <- model_SurRecr_helped$Sol
colnames(post)
var_list <- colnames(post)[!grepl("^traitRecruits1:", colnames(post))]
var_list <- var_list[-c(1,2)]  
print(var_list)
# 创建一个函数来提取每个变量对第二性状的总效应
# 总效应 = 主效应 + traitRecruits1:主效应
get_total_effect <- function(varname) {
  main <- post[, varname]
  interaction <- post[, paste0("traitRecruits1:", varname)]
  total <- main + interaction
  return(total)
}
# 存储结果
results <- data.frame(
  Variable = character(),
  Mean = numeric(),
  Lower = numeric(),
  Upper = numeric(),
  pMCMC = numeric(),
  stringsAsFactors = FALSE
)

# 逐个变量处理
for (var in var_list) {
  total_effect <- get_total_effect(var)
  hpd <- HPDinterval(as.mcmc(total_effect))
  pval <- 2 * min(mean(total_effect > 0), mean(total_effect < 0))
  
  results <- rbind(results, data.frame(
    Variable = var,
    Mean = mean(total_effect),
    Lower = hpd[1],
    Upper = hpd[2],
    pMCMC = pval
  ))
}
# 查看结果
print(results)

#model_SurRecr_nonhelped
model_SurRecr_nonhelped<-readRDS( "model_SurRecr_nonhelped.rds") 
# 提取后验样本
post <- model_SurRecr_nonhelped$Sol
colnames(post)
var_list <- colnames(post)[!grepl("^traitRecruits1:", colnames(post))]
var_list <- var_list[-c(1,2)]  
print(var_list)
# 创建一个函数来提取每个变量对第二性状的总效应
# 总效应 = 主效应 + traitRecruits1:主效应
get_total_effect <- function(varname) {
  main <- post[, varname]
  interaction <- post[, paste0("traitRecruits1:", varname)]
  total <- main + interaction
  return(total)
}
# 存储结果
results <- data.frame(
  Variable = character(),
  Mean = numeric(),
  Lower = numeric(),
  Upper = numeric(),
  pMCMC = numeric(),
  stringsAsFactors = FALSE
)

# 逐个变量处理
for (var in var_list) {
  total_effect <- get_total_effect(var)
  hpd <- HPDinterval(as.mcmc(total_effect))
  pval <- 2 * min(mean(total_effect > 0), mean(total_effect < 0))
  
  results <- rbind(results, data.frame(
    Variable = var,
    Mean = mean(total_effect),
    Lower = hpd[1],
    Upper = hpd[2],
    pMCMC = pval
  ))
}
# 查看结果
print(results)


###model_OffSurFle
# 提取后验样本
post <- model_OffSurFle$Sol
colnames(post)
var_list <- colnames(post)[!grepl("^traitbrood:", colnames(post))]
var_list <- var_list[-c(1,2)]  
print(var_list)
# 创建一个函数来提取每个变量对第二性状的总效应
# 总效应 = 主效应 + traitbrood:主效应
get_total_effect <- function(varname) {
  main <- post[, varname]
  interaction <- post[, paste0("traitbrood:", varname)]
  total <- main + interaction
  return(total)
}
# 存储结果
results <- data.frame(
  Variable = character(),
  Mean = numeric(),
  Lower = numeric(),
  Upper = numeric(),
  pMCMC = numeric(),
  stringsAsFactors = FALSE
)
# 逐个变量处理
for (var in var_list) {
  total_effect <- get_total_effect(var)
  hpd <- HPDinterval(as.mcmc(total_effect))
  pval <- 2 * min(mean(total_effect > 0), mean(total_effect < 0))
  
  results <- rbind(results, data.frame(
    Variable = var,
    Mean = mean(total_effect),
    Lower = hpd[1],
    Upper = hpd[2],
    pMCMC = pval
  ))
}
# 查看结果
print(results)

###model_OffSurFle_helped

post <- model_OffSurFle_helped$Sol
colnames(post)
var_list <- colnames(post)[!grepl("^traitbrood:", colnames(post))]
var_list <- var_list[-c(1,2)]  
print(var_list)
# 创建一个函数来提取每个变量对第二性状的总效应
# 总效应 = 主效应 + traitbrood:主效应
get_total_effect <- function(varname) {
  main <- post[, varname]
  interaction <- post[, paste0("traitbrood:", varname)]
  total <- main + interaction
  return(total)
}
# 存储结果
results <- data.frame(
  Variable = character(),
  Mean = numeric(),
  Lower = numeric(),
  Upper = numeric(),
  pMCMC = numeric(),
  stringsAsFactors = FALSE
)
# 逐个变量处理
for (var in var_list) {
  total_effect <- get_total_effect(var)
  hpd <- HPDinterval(as.mcmc(total_effect))
  pval <- 2 * min(mean(total_effect > 0), mean(total_effect < 0))
  
  results <- rbind(results, data.frame(
    Variable = var,
    Mean = mean(total_effect),
    Lower = hpd[1],
    Upper = hpd[2],
    pMCMC = pval
  ))
}
# 查看结果
print(results)

###model_OffSurFle_nonhelped
post <- model_OffSurFle_nonhelped$Sol
colnames(post)
var_list <- colnames(post)[!grepl("^traitbrood:", colnames(post))]
var_list <- var_list[-c(1,2)]  
print(var_list)
# 创建一个函数来提取每个变量对第二性状的总效应
# 总效应 = 主效应 + traitbrood:主效应
get_total_effect <- function(varname) {
  main <- post[, varname]
  interaction <- post[, paste0("traitbrood:", varname)]
  total <- main + interaction
  return(total)
}
# 存储结果
results <- data.frame(
  Variable = character(),
  Mean = numeric(),
  Lower = numeric(),
  Upper = numeric(),
  pMCMC = numeric(),
  stringsAsFactors = FALSE
)
# 逐个变量处理
for (var in var_list) {
  total_effect <- get_total_effect(var)
  hpd <- HPDinterval(as.mcmc(total_effect))
  pval <- 2 * min(mean(total_effect > 0), mean(total_effect < 0))
  
  results <- rbind(results, data.frame(
    Variable = var,
    Mean = mean(total_effect),
    Lower = hpd[1],
    Upper = hpd[2],
    pMCMC = pval
  ))
}
# 查看结果
print(results)

#比较合作与非合作个体观测层面上双变量协变的差异
vcv_helped <- model_SurRecr_helped$VCV[, "traitRecruits1:traitSurvival.1.units"]
vcv_nonhelped <- model_SurRecr_nonhelped$VCV[, "traitRecruits1:traitSurvival.1.units"]
# 后验差值
vcv_diff <- vcv_helped - vcv_nonhelped
# 后验概率：差值是否显著偏离 0
p_value <- 2 * min(mean(vcv_diff > 0), mean(vcv_diff < 0))  # 双尾检验p = 0.810625
HPDinterval(vcv_diff) #[-48.81585, 65.03693]



