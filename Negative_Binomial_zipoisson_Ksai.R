############################################################################################
#绘制合作与非合作各自的ksai随降水变化的图
precip_data <- read.csv("cleaned_dataday_pre.csv")
temp_data <- read.csv("cleaned_climlongday.csv")

# 将日期列转换为正确的日期格式 (dd/mm/yyyy)
precip_data$Date <- as.Date(precip_data$Date, format = "%d/%m/%Y")
temp_data$Date <- as.Date(temp_data$Date, format = "%d/%m/%Y")

#合并环境数据
library(dplyr)
# 计算每年11月的平均温度
temp_nov_dec <- temp_data %>%
  filter(format(Date, "%m") %in% c("11")) %>%  # 筛选11月
  mutate(year = format(Date, "%Y")) %>%  # 提取年份
  group_by(year) %>%
  summarise(avg_temp_nov_dec = mean(tmpvalue, na.rm = TRUE))  # 计算平均温度
# 计算每年4月和5月的总降水量
precip_apr_may <- precip_data %>%
  filter(format(Date, "%m") %in% c("04", "05")) %>%  # 筛选4月和5月
  mutate(year = format(Date, "%Y")) %>%  # 提取年份
  group_by(year) %>%
  summarise(total_precip_apr_may = sum(Pre, na.rm = TRUE))  # 计算总降水量
# 合并到数据集a中
a<- read.csv("dx_regression_Pt.csv")
a<- a[a$state != "H", ]
# 确保a数据框中有年份列
a$year <- as.character(a$year)  # 确保年份格式一致
# 合并温度数据和降水数据到a中
a_merged <- a %>%
  left_join(temp_nov_dec, by = "year") %>%  # 合并平均温度
  left_join(precip_apr_may, by = "year")  # 合并总降水量
# 查看合并后的数据框
head(a_merged)
median_temp <- median(a_merged$avg_temp_nov_dec, na.rm = TRUE)
median_precip <- median(a_merged$total_precip_apr_may, na.rm = TRUE)
# 创建高温/低温和高降水/低降水分类
a_merged <- a_merged %>%
  mutate(temp_category = ifelse(avg_temp_nov_dec >= median_temp, "High Temp", "Low Temp"),
         precip_category = ifelse(total_precip_apr_may >= median_precip, "High Precip", "Low Precip"))
# 组合四组环境条件
a_merged <- a_merged %>%
  mutate(environment = paste(temp_category, precip_category, sep = " - "))
head(a_merged)
a_merged <- a_merged %>%
  mutate(
    Age= scale(Age),
    Ter = scale(Ter.),
    avg_temp_nov_dec= scale(avg_temp_nov_dec),
    total_precip_apr_may = scale(total_precip_apr_may)
  )
a_merged$ksai_1 <- a_merged$Survival + 2 * a_merged$Recruits
print(unique(a_merged$ksai_1))
#############################################ksai#############################################
a_merged_ExH_numeric<- a_merged %>% 
  filter(state != "H")%>%
  mutate(ksai_1 = as.integer(ksai_1),
         Recruits1 = Recruits * 2,
         Gethp = as.numeric(Gethp))
a_merged_ExH_numeric$Age[is.na(a_merged_ExH_numeric$Age)] <- median(a_merged_ExH_numeric$Age, na.rm = TRUE)
colSums(is.na(a_merged_ExH_numeric))
a_merged_ExH_numeric_Non <- a_merged_ExH_numeric %>%
  filter(Gethp == 0)
a_merged_ExH_numeric_Coo <- a_merged_ExH_numeric %>%
  filter(Gethp == 1)

library(MCMCglmm)
##全模型
prior_zipoisson <- list(
  R = list(V = diag(2), nu = 0.002, fix = 2),  # 残差协方差结构固定
  G = list(
    G1 = list(V = diag(2), nu = 2),  # us(trait):fn
    G2 = list(V = diag(2), nu = 2)   # us(trait):id
  )
)

# 拟合零膨胀泊松模型
modelKsai_zipoisson_full <- MCMCglmm(
  fixed = ksai_1 ~ trait - 1 + 
    at.level(trait, 1):(Gethp + Age + Ter. + sex + 
                          scale(total_precip_apr_may) + scale(avg_temp_nov_dec) + 
                          Gethp:scale(total_precip_apr_may) + 
                          Gethp:scale(avg_temp_nov_dec) + 
                          scale(total_precip_apr_may):scale(avg_temp_nov_dec) + 
                          Gethp:scale(total_precip_apr_may):scale(avg_temp_nov_dec)) +
    at.level(trait, 2):(Gethp + Age + Ter. + sex + 
                          scale(total_precip_apr_may) + scale(avg_temp_nov_dec) + 
                          Gethp:scale(total_precip_apr_may) + 
                          Gethp:scale(avg_temp_nov_dec) + 
                          scale(total_precip_apr_may):scale(avg_temp_nov_dec) + 
                          Gethp:scale(total_precip_apr_may):scale(avg_temp_nov_dec)),
  
  random = ~ us(trait):fn + us(trait):id,  # 加入个体ID作为随机效应
  rcov   = ~ idh(trait):units,             # 残差结构
  family = "zipoisson",
  data   = a_merged_ExH_numeric,
  prior  = prior_zipoisson,
  pl     = TRUE,
  slice  = FALSE,
  nitt   = 2000000,
  burnin = 500000,
  thin   = 1000,
  verbose = FALSE
)

summary(modelKsai_zipoisson_full)

# 拟合负二项分布模型
# 设置先验：允许残差估计，并为三个随机项（ID、fn）分别设置G结构
prior_nb <- list(
  R = list(V = 1, nu = 0.002),  # 残差结构，必须可估
  G = list(
    G1 = list(V = 1, nu = 0.002),   
    G2 = list(V = 1, nu = 0.002)   
  )
)
# 拟合近似负二项分布的模型
modelKsai_nb <- MCMCglmm(
  fixed = ksai_1 ~ Gethp + Age + Ter. + sex + 
    scale(total_precip_apr_may) + scale(avg_temp_nov_dec) + 
    Gethp:scale(total_precip_apr_may) +
    Gethp:scale(avg_temp_nov_dec) +
    scale(total_precip_apr_may):scale(avg_temp_nov_dec) +
    Gethp:scale(total_precip_apr_may):scale(avg_temp_nov_dec),
  
  random = ~  id + fn,                 # 多个随机效应项
  family = "poisson",                        # 使用可估残差模拟负二项
  data   = a_merged_ExH_numeric,
  prior  = prior_nb,
  nitt   = 2000000,
  burnin = 500000,
  thin   = 1000,
  verbose = FALSE
)

##在只包含非合作的数据集中，看看非合作的跟随效应


modelKsai_NonCoo<-  MCMCglmm(
  fixed = ksai_1 ~ trait - 1 + 
    at.level(trait, 1):(Age + Ter. + sex + 
                          scale(total_precip_apr_may) + scale(avg_temp_nov_dec) + 
                          Gethp:scale(total_precip_apr_may) + 
                          Gethp:scale(avg_temp_nov_dec) + 
                          scale(total_precip_apr_may):scale(avg_temp_nov_dec) + 
                          Gethp:scale(total_precip_apr_may):scale(avg_temp_nov_dec)) +
    at.level(trait, 2):(Age + Ter. + sex + 
                          scale(total_precip_apr_may) + scale(avg_temp_nov_dec) + 
                          Gethp:scale(total_precip_apr_may) + 
                          Gethp:scale(avg_temp_nov_dec) + 
                          scale(total_precip_apr_may):scale(avg_temp_nov_dec) + 
                          Gethp:scale(total_precip_apr_may):scale(avg_temp_nov_dec)),
  
  random = ~ us(trait):fn + us(trait):id,  # 加入个体ID作为随机效应
  rcov   = ~ idh(trait):units,             # 残差结构
  family = "zipoisson",
  data   = a_merged_ExH_numeric_Non,
  prior  = prior_zipoisson,
  pl     = TRUE,
  slice  = FALSE,
  nitt   = 2000000,
  burnin = 500000,
  thin   = 1000,
  verbose = FALSE
)

# 拟合负二项分布模型
# 拟合近似负二项分布的模型
modelKsai_nb_NonCoo<- MCMCglmm(
  fixed = ksai_1 ~ Age + Ter. + sex + 
    scale(total_precip_apr_may) + scale(avg_temp_nov_dec) + 
    Gethp:scale(total_precip_apr_may) +
    Gethp:scale(avg_temp_nov_dec) +
    scale(total_precip_apr_may):scale(avg_temp_nov_dec) +
    Gethp:scale(total_precip_apr_may):scale(avg_temp_nov_dec),
  
  random = ~  id + fn,                 # 多个随机效应项
  family = "poisson",                        # 使用可估残差模拟负二项
  data   = a_merged_ExH_numeric_Non,
  prior  = prior_nb,
  nitt   = 2000000,
  burnin = 500000,
  thin   = 1000,
  verbose = FALSE
)

##在只包含合作的数据集中，看看合作的跟随效应


modelKsai_Coo<-  MCMCglmm(
  fixed = ksai_1 ~ trait - 1 + 
    at.level(trait, 1):(Age + Ter. + sex + 
                          scale(total_precip_apr_may) + scale(avg_temp_nov_dec) + 
                          Gethp:scale(total_precip_apr_may) + 
                          Gethp:scale(avg_temp_nov_dec) + 
                          scale(total_precip_apr_may):scale(avg_temp_nov_dec) + 
                          Gethp:scale(total_precip_apr_may):scale(avg_temp_nov_dec)) +
    at.level(trait, 2):(Age + Ter. + sex + 
                          scale(total_precip_apr_may) + scale(avg_temp_nov_dec) + 
                          Gethp:scale(total_precip_apr_may) + 
                          Gethp:scale(avg_temp_nov_dec) + 
                          scale(total_precip_apr_may):scale(avg_temp_nov_dec) + 
                          Gethp:scale(total_precip_apr_may):scale(avg_temp_nov_dec)),
  
  random = ~ us(trait):fn + us(trait):id,  # 加入个体ID作为随机效应
  rcov   = ~ idh(trait):units,             # 残差结构
  family = "zipoisson",
  data   = a_merged_ExH_numeric_Coo ,
  prior  = prior_zipoisson,
  pl     = TRUE,
  slice  = FALSE,
  nitt   = 2000000,
  burnin = 500000,
  thin   = 1000,
  verbose = FALSE
)

# 拟合负二项分布模型
# 拟合近似负二项分布的模型
modelKsai_nb_Coo<- MCMCglmm(
  fixed = ksai_1 ~ Age + Ter. + sex + 
    scale(total_precip_apr_may) + scale(avg_temp_nov_dec) + 
    Gethp:scale(total_precip_apr_may) +
    Gethp:scale(avg_temp_nov_dec) +
    scale(total_precip_apr_may):scale(avg_temp_nov_dec) +
    Gethp:scale(total_precip_apr_may):scale(avg_temp_nov_dec),
  
  random = ~  id + fn,                 # 多个随机效应项
  family = "poisson",                        # 使用可估残差模拟负二项
  data   = a_merged_ExH_numeric_Coo ,
  prior  = prior_nb,
  nitt   = 2000000,
  burnin = 500000,
  thin   = 1000,
  verbose = FALSE
)

#种群整体降水对ksai的效应趋势

modelKsai_population<-  MCMCglmm(
  fixed = ksai_1 ~ trait - 1 + 
    at.level(trait, 1):(Gethp + Age + Ter. + sex + 
                          scale(total_precip_apr_may) + scale(avg_temp_nov_dec) + 
                          scale(total_precip_apr_may):scale(avg_temp_nov_dec)) +
    at.level(trait, 2):(Gethp + Age + Ter. + sex + 
                          scale(total_precip_apr_may) + scale(avg_temp_nov_dec) + 
                          scale(total_precip_apr_may):scale(avg_temp_nov_dec)),
  
  random = ~ us(trait):fn + us(trait):id,  # 加入个体ID作为随机效应
  rcov   = ~ idh(trait):units,             # 残差结构
  family = "zipoisson",
  data   = a_merged_ExH_numeric,
  prior  = prior_zipoisson,
  pl     = TRUE,
  slice  = FALSE,
  nitt   = 2000000,
  burnin = 500000,
  thin   = 1000,
  verbose = FALSE
)

# 拟合负二项分布模型
# 拟合近似负二项分布的模型
modelKsai_nb_population<- MCMCglmm(
  fixed = ksai_1 ~ Gethp + Age + Ter. + sex + 
    scale(total_precip_apr_may) + scale(avg_temp_nov_dec) + 
    
    scale(total_precip_apr_may):scale(avg_temp_nov_dec),
  
  random = ~  id + fn,                 # 多个随机效应项
  family = "poisson",                        # 使用可估残差模拟负二项
  data   = a_merged_ExH_numeric,
  prior  = prior_nb,
  nitt   = 2000000,
  burnin = 500000,
  thin   = 1000,
  verbose = FALSE
)

# 确保 Gethp 是一个因子，以便 MCMCglmm 正确处理其系数
a_merged_ExH_numeric$Gethp <- factor(a_merged_ExH_numeric$Gethp)

# 构建整合模型：Gethp 作为固定效应
modelKsai_nb_Combined_simple <- MCMCglmm(
  fixed = ksai_1 ~ Gethp + Age + Ter. + sex,  # <--- Gethp 是关键的比较变量
  random = ~ id + year,  
  family = "poisson",    
  data  = a_merged_ExH_numeric,             # <--- 使用完整数据
  prior = prior_nb,
  nitt  = 2000000,
  burnin = 500000,
  thin  = 1000,
  verbose = FALSE
)

# 查看结果
summary(modelKsai_nb_Combined_simple)

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble) # 用于更可靠地创建数据框

# --- I. MCMC 样本计算与摘要 ---

# 1. 提取固定效应的后验样本
mcmc_samples <- modelKsai_nb_Combined_simple$Sol

# 2. 计算非合作组 (Gethp=0) 的后验分布 (log 尺度)
log_ksai_noncoop <- mcmc_samples[, "(Intercept)"]

# 3. 计算合作组 (Gethp=1) 的后验分布 (log 尺度)
log_ksai_coop <- mcmc_samples[, "(Intercept)"] + mcmc_samples[, "Gethp1"]

# 4. 对 log 尺度上的后验样本进行反向链接 (指数化)
ksai_noncoop_posterior <- exp(log_ksai_noncoop)
ksai_coop_posterior <- exp(log_ksai_coop)

# 5. 计算最终的后验均值 (Estimated Mean ξ) 和 CI
estimated_mean_ksai <- data.frame(
  Group = c("Non-cooperative", "Cooperative"),
  Post_Mean_ksai = c(mean(ksai_noncoop_posterior), mean(ksai_coop_posterior)),
  Lower_95_CI = c(quantile(ksai_noncoop_posterior, 0.025), quantile(ksai_coop_posterior, 0.025)),
  Upper_95_CI = c(quantile(ksai_noncoop_posterior, 0.975), quantile(ksai_coop_posterior, 0.975))
)

print(estimated_mean_ksai)

# --- II. 数据整理 (绘制密度图) ---

posterior_df <- bind_rows(
  tibble(
    Group = "Non_cooperative", 
    ksai_posterior = ksai_noncoop_posterior
  ),
  tibble(
    Group = "Cooperative",
    ksai_posterior = ksai_coop_posterior
  )
) %>%
  mutate(
      Group = factor(Group, levels = c("Non_cooperative", "Cooperative"))
  )

posterior_means <- posterior_df %>%
  group_by(Group) %>%
  summarise(
    mean_ksai = mean(ksai_posterior),
    .groups = 'drop'
  )

# --- III. 绘图 (Figure 2a) ---
library(ggplot2)
library(dplyr)
library(tidyr)

color_mapping <- c("Non_cooperative" = "steelblue", "Cooperative" = "tomato")
ggplot(posterior_df, aes(x = ksai_posterior, color = Group, fill = Group)) +
  geom_density(alpha = 0.3, linewidth = 1, bw = 0.1) + 
  geom_vline(data = posterior_means, aes(xintercept = mean_ksai, color = Group), 
             linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  
  labs(
    x = expression(xi), # x轴标签简化为 xi，保持风格
    y = "Density",
        title = expression("Posterior distribution of mean " * xi)
  ) +
      theme_classic(base_size = 14) +
   theme(
    legend.title = element_blank(),
    legend.position = c(0.8, 0.8), # 图例位置
     )

###对合作组和非合作组的模型进行重构#
#将气候因子去掉，以年为随机效应来估算环境随机性,个体ID被作为随机效应加以控制
modelKsai_nb_Coo_simple <- MCMCglmm(
  fixed = ksai_1 ~ Age + Ter. + sex,
  random = ~ id + year,  
  family = "poisson",  
  data   = a_merged_ExH_numeric_Coo,
  prior  = prior_nb,
  nitt   = 2000000,
  burnin = 500000,
  thin   = 1000,
  verbose = FALSE
)

summary(modelKsai_nb_Coo_simple)


modelKsai_nb_NonCoo_simple <- MCMCglmm(
  fixed = ksai_1 ~ Age + Ter. + sex,
  random = ~ id + year,  
  family = "poisson",  
  data   = a_merged_ExH_numeric_Non,
  prior  = prior_nb,
  nitt   = 2000000,
  burnin = 500000,
  thin   = 1000,
  verbose = FALSE
)

summary(modelKsai_nb_NonCoo_simple)

head(modelKsai_nb_Coo_simple$VCV)

# 提取 posterior samples
year_coo  <- modelKsai_nb_Coo_simple$VCV[,"year"]
year_non  <- modelKsai_nb_NonCoo_simple$VCV[,"year"]
units_coo <- modelKsai_nb_Coo_simple$VCV[,"units"]
units_non <- modelKsai_nb_NonCoo_simple$VCV[,"units"]
r_coo  <- modelKsai_nb_Coo_simple$Sol[, "(Intercept)"]
r_non  <- modelKsai_nb_NonCoo_simple$Sol[, "(Intercept)"]

sigma_env_coo  <- (exp(year_coo)  - 1) * exp(2 * r_coo + year_coo)
sigma_env_non  <- (exp(year_non)  - 1) * exp(2 * r_non + year_non)
sigma_stat_coo <- (exp(units_coo) - 1) * exp(2 * r_coo + units_coo)
sigma_stat_non <- (exp(units_non) - 1) * exp(2 * r_non + units_non)

# 差值分布
year_diff  <- sigma_env_coo - sigma_env_non     # 合作 - 非合作（环境随机性差值）
units_diff <- sigma_stat_coo - sigma_stat_non   # 合作 - 非合作（统计随机性差值）

# 后验差值大于0的比例
mean(year_diff)  
mean(units_diff)
mean(year_diff > 0)   # > 0 表示合作年间异质性更大
mean(units_diff > 0)  # > 0 表示合作统计方差更大
HPDinterval(as.mcmc(year_diff)) #合作的环境随机性更大
HPDinterval(as.mcmc(units_diff)) #两者的统计随机性不分上下
# 单尾概率（合作组更大）
pMCMC_year_one_tail  <- 1 - mean(year_diff > 0)
pMCMC_year_one_tail 
pMCMC_units_one_tail <- 1 - mean(units_diff > 0)
pMCMC_units_one_tail

# 种群规模N
a <- read.csv("dx_regression_Pt.csv")
# 每年样本数量
table_per_year <- table(a$year)
print(table_per_year)
# 平均每年样本数
N <- mean(as.numeric(table_per_year))
cat("平均每年样本数为：", round(N, 2), "\n")
# 计算每组的随机增长率 a
a_coo <- r_coo - year_coo / 2 - units_coo / (2 * N)
a_non <- r_non - year_non / 2 - units_non / (2 * N)

# 合作组的随机增长率
mean_coo <- mean(a_coo)
ci_coo <- HPDinterval(as.mcmc(a_coo))

# 非合作组的随机增长率
mean_non <- mean(a_non)
ci_non <- HPDinterval(as.mcmc(a_non))

# 差值（合作 - 非合作）
a_diff <- a_coo - a_non

# 比较两者差异
mean(a_diff > 0)  # >0 表示合作组的随机增长率更高
HPDinterval(as.mcmc(a_diff))  # 95%置信区间

# 单尾 p 值（合作组更高）
pMCMC_a_one_tail <- 1 - mean(a_diff > 0)

# 输出
cat("合作组随机增长率显著高于非合作组的比例：", mean(a_diff > 0), "\n")
cat("95% HPD区间：\n")
print(HPDinterval(as.mcmc(a_diff)))
cat("单尾 pMCMC 值：", round(pMCMC_a_one_tail, 4), "\n")

library(dplyr)

# 假设母数据集 a_merged_ExH_numeric 已经加载到 R 环境中。

# 1. 计算冬季温度（avg_temp_nov_dec）的全局均值
# 使用 na.rm = TRUE 确保在有缺失值的情况下也能正确计算均值
temp_cutoff <- mean(a_merged_ExH_numeric$avg_temp_nov_dec, na.rm = TRUE)

# 2. 创建高温（暖冬）数据集
# 包含所有 avg_temp_nov_dec 大于或等于均值的观测值
a_merged_warm_winter <- a_merged_ExH_numeric %>%
  filter(avg_temp_nov_dec >= temp_cutoff)

# 3. 创建低温（冷冬）数据集
# 包含所有 avg_temp_nov_dec 小于均值的观测值
a_merged_cold_winter <- a_merged_ExH_numeric %>%
  filter(avg_temp_nov_dec < temp_cutoff)

# 4. (可选) 检查子集大小
cat("母数据集观测数量:", nrow(a_merged_ExH_numeric), "\n")
cat("高温（暖冬）数据集观测数量:", nrow(a_merged_warm_winter), "\n")
cat("低温（冷冬）数据集观测数量:", nrow(a_merged_cold_winter), "\n")

# 2. 拟合暖冬（高温）子模型 
cat("Fitting Warm Winter Model...\n")
modelKsai_warm_winter <- MCMCglmm(
  fixed = ksai_1 ~ Gethp + Age + Ter. + sex +
    scale(total_precip_apr_may) + 
    Gethp:scale(total_precip_apr_may),
  
  random = ~ id + fn,
  family = "poisson",
  data   = a_merged_warm_winter, # 使用暖冬数据集
  prior  = prior_nb,             # 使用你定义的先验
  nitt   = 2000000,
  burnin = 500000,
  thin   = 1000,
  verbose = FALSE
)
cat("Warm Winter Model Fitting Complete.\n")


# 3. 拟合冷冬（低温）子模型
cat("Fitting Cold Winter Model...\n")
modelKsai_cold_winter <- MCMCglmm(
  fixed = ksai_1 ~ Gethp + Age + Ter. + sex +
    scale(total_precip_apr_may) + 
    Gethp:scale(total_precip_apr_may),
  
  random = ~ id + fn,
  family = "poisson",
  data   = a_merged_cold_winter, # 使用冷冬数据集
  prior  = prior_nb,             # 使用你定义的先验
  nitt   = 2000000,
  burnin = 500000,
  thin   = 1000,
  verbose = FALSE
)
cat("Cold Winter Model Fitting Complete.\n")


# 4. 检查结果
# 使用 summary() 函数检查每个模型的收敛性和系数
summary(modelKsai_warm_winter)
summary(modelKsai_cold_winter)

library(MCMCglmm)
library(dplyr)
library(ggplot2)

# ----------------------------------------------------------------------
# 预测 a_merged_warm_winter
a_merged_warm_winter$ksai_pred_mean<- predict(modelKsai_warm_winter, newdata = a_merged_warm_winter, type = "response")
head(a_merged_warm_winter)
# 绘图
# 1. 准备绘图数据框
plot_df_obs <- a_merged_warm_winter %>%
  mutate(
    Cooperation = factor(Gethp, levels = c(1, 0), labels = c("Cooperative", "Non-cooperative"))
  )

# 自定义颜色
color_values <- c(
  "Cooperative" = "#FF5757",      # 红色
  "Non-cooperative" = "#0070C0"    # 蓝色
)

# 2. 绘图：使用 geom_smooth 自动拟合回归线
library(ggplot2)
library(dplyr)

# 假设 plot_df_obs 包含所有数据点（合作和非合作）
# 假设 color_values 已经定义

ggplot(plot_df_obs, aes(x = total_precip_apr_may, y = ksai_pred_mean)) +
  
  # 1. 实际数据点
  geom_point(aes(color = Cooperation, shape = Cooperation), # 添加 shape 映射
             alpha = 0.3, size = 1) +
  
  # 2. 回归线和置信区间
  geom_smooth(aes(color = Cooperation), # 移除 fill 映射
              fill = "gray70",          # 统一使用灰色阴影
              method = "lm", 
              se = TRUE,              # 显示置信区间
              linewidth = 1.2, 
              linetype = "solid") + 
  
  # 手动设置样式
  scale_color_manual(values = color_values, name = "Phenotype") +
  # 移除 scale_fill_manual
  
  # 新增：手动设置点型
  scale_shape_manual(values = c("Cooperative" = 19, "Non-cooperative" = 17), name = "Phenotype") +
  
  # 标签与主题
  labs(
    x = "Spring Rainfall (mm, April–May)",
    y = expression(paste("Predicted Cooperation Effect (", xi, ")")), 
    title = "Cooperation Effect in Warmer Winters"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold")
  )




saveRDS(modelKsai_zipoisson_full,  file = "modelKsai_zipoisson_full.rds")
saveRDS(modelKsai_nb,              file = "modelKsai_nb.rds")
saveRDS(modelKsai_NonCoo,          file = "modelKsai_NonCoo.rds")
saveRDS(modelKsai_nb_NonCoo,       file = "modelKsai_nb_NonCoo.rds")
saveRDS(modelKsai_Coo,             file = "modelKsai_Coo.rds")
saveRDS(modelKsai_nb_Coo,          file = "modelKsai_nb_Coo.rds")
saveRDS(modelKsai_population,      file = "modelKsai_population.rds")
saveRDS(modelKsai_nb_population,   file = "modelKsai_nb_population.rds")
saveRDS(modelKsai_nb_Coo_simple,   file = "modelKsai_nb_Coo_simple.rds")
saveRDS(modelKsai_nb_NonCoo_simple,file = "modelKsai_nb_NonCoo_simple.rds")