library(MCMCglmm)
library(dplyr)

Nest_Fle<- read.csv("NESTlEVEL.csv")
offspring_Sur<-read.csv("DX_offss.csv")  
precip_data <- read.csv("cleaned_dataday_pre.csv")
temp_data <- read.csv("cleaned_climlongday.csv")
length(unique(Nest_Fle$ID)) #窝雏数数据的巢数644
length(unique(offspring_Sur$ID))#幼体存活数据2426个
length(unique(offspring_Sur$fn))#幼体存活数据分属541巢

#数据检查
head(Nest_Fle)
head(offspring_Sur)
fn_vals <- unique(offspring_Sur$fn)
missing_fn <- fn_vals[!fn_vals %in% Nest_Fle$ID]
missing_fn
id_vals <- unique(Nest_Fle$ID)
missing_id <- id_vals[!id_vals %in% offspring_Sur$fn]
missing_id
#繁殖失败的有104巢，成功出飞的有540巢
Nest_Fle %>%
    mutate(brood_status = ifelse(brood == 0, "brood_zero", "brood_non_zero")) %>%
    count(brood_status)


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
str(a_OffSurFle_clean)

#巢水平数据制作
nest_level_data <- a_OffSurFle_clean %>%
  group_by(fn, Year) %>%
  summarise(
    SUR_total = sum(SUR, na.rm = TRUE),
    brood = first(brood),
    Gethp = first(Gethp),
    Exp = first(Exp.),
    Terr = first(Terr.),
    avg_temp_nov = first(avg_temp_nov),
    total_precip_apr_may = first(total_precip_apr_may),
    .groups = "drop"
  )

head(nest_level_data)
# 统计 brood 列中 0 和 非 0 的数量
cat("brood 列的统计结果 (TRUE = 0, FALSE = 非 0):\n")
print(table(nest_level_data$brood == 0))

# 统计 SUR_total 列中 0 和 非 0 的数量
cat("\nSUR_total 列的统计结果 (TRUE = 0, FALSE = 非 0):\n")
print(table(nest_level_data$SUR_total == 0))

# Gethp = 1 的个体
data_helped <- as.data.frame(nest_level_data %>% filter(Gethp == 1))
nrow(data_helped)
# Gethp = 0 的个体
data_nonhelped <- as.data.frame(nest_level_data %>% filter(Gethp == 0))
nrow(data_nonhelped)

# 设定先验，允许残差协方差矩阵
prior <- list(
  R = list(V = diag(2), nu = 0.002),  # 残差协方差矩阵
  G = list(
    G1 = list(V = diag(2), nu = 0.002)  # 随机效应
  )
)

model_OffSurFle_Nest <- MCMCglmm(
  cbind(SUR_total, brood) ~ trait * (Gethp + Exp + Terr + avg_temp_nov + total_precip_apr_may +
                                 avg_temp_nov:total_precip_apr_may +
                                 Gethp:avg_temp_nov + Gethp:total_precip_apr_may +
                                 Gethp:avg_temp_nov:total_precip_apr_may),
  random = ~ us(trait):Year ,  
  rcov = ~ us(trait):units,
  family = c("poisson", "poisson"),
  data = as.data.frame(nest_level_data),
  prior = prior,
  nitt = 720000, burnin = 144000, thin = 30,
  verbose = FALSE
)
summary(model_OffSurFle_Nest)


model_OffSurFle_Nest_helped <- MCMCglmm(
  cbind(SUR_total, brood) ~ trait * (Exp + Terr + avg_temp_nov + total_precip_apr_may +
                                 avg_temp_nov:total_precip_apr_may),
  random = ~ us(trait):Year,
  rcov = ~ us(trait):units,
  family = c("poisson", "poisson"),
  data = data_helped,
  prior = prior,
  nitt = 720000, burnin = 144000, thin = 30,
  verbose = FALSE
)
summary(model_OffSurFle_Nest_helped)


model_OffSurFle_Nest_Non <- MCMCglmm(
  cbind(SUR_total, brood) ~ trait * (Exp + Terr + avg_temp_nov + total_precip_apr_may +
                                 avg_temp_nov:total_precip_apr_may),
  random = ~ us(trait):Year,
  rcov = ~ us(trait):units,
  family = c("poisson", "poisson"),
  data = data_nonhelped,
  prior = prior,
  nitt = 720000, burnin = 144000, thin = 30,
  verbose = FALSE
)
summary(model_OffSurFle_Nest_Non)

###model_OffSurFle_Nest
# 提取后验样本
post <- model_OffSurFle_Nest$Sol
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

###model_OffSurFle_Nest_helped

post <- model_OffSurFle_Nest_helped$Sol
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

###model_OffSurFle_Nest_Non
post <- model_OffSurFle_Nest_Non$Sol
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

saveRDS(model_OffSurFle_Nest, file = "model_OffSurFle_Nest.rds")
saveRDS(model_OffSurFle_Nest_helped, file = "model_OffSurFle_Nest_helped.rds")
saveRDS(model_OffSurFle_Nest_Non, file = "model_OffSurFle_Nest_Non.rds")


