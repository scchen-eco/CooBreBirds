library(MCMCglmm)
library(ggplot2)
library(dplyr)


#环境随机性和统计随机性的频数图绘制

modelKsai_nb_Coo_simple         <- readRDS("modelKsai_nb_Coo_simple.rds")
modelKsai_nb_NonCoo_simple      <- readRDS("modelKsai_nb_NonCoo_simple.rds")

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
mean(a_diff)
HPDinterval(as.mcmc(a_diff))  # 95%置信区间

# 单尾 p 值（合作组更高）
pMCMC_a_one_tail <- 1 - mean(a_diff > 0)

# 输出
cat("合作组随机增长率显著高于非合作组的比例：", mean(a_diff > 0), "\n")
cat("95% HPD区间：\n")
print(HPDinterval(as.mcmc(a_diff)))
cat("单尾 pMCMC 值：", round(pMCMC_a_one_tail, 4), "\n")



library(ggplot2)
library(dplyr)

#统计随机性密度图
df_stat <- data.frame(
  value = c(sigma_stat_coo, sigma_stat_non),
  group = rep(c("Cooperative", "Non-cooperative"), each = length(sigma_env_coo))
)

# 计算均值用于画虚线
#means <- df_stat %>%
#  group_by(group) %>%
#  summarise(mean_value = mean(value))

# 绘图
ggplot(df_stat, aes(x = value, color = group, fill = group)) +
  geom_density(alpha = 0.3, linewidth = 1,bw = 0.01) +
 # geom_vline(data = means, aes(xintercept = mean_value, color = group),
             #linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = c("tomato", "steelblue")) +
  scale_fill_manual(values = c("tomato", "steelblue")) +
  labs(x = expression(sigma[demographic]),
       y = "Density",
       title = expression("Posterior distribution of " * sigma[demographic])) +
  #coord_cartesian(xlim = c(0, 1)) +  # 限制 x 轴范围
  theme_classic(base_size = 14)

#环境随机性密度图
df_env <- data.frame(
  value = c(sigma_env_coo, sigma_env_non),
  group = rep(c("Cooperative", "Non-cooperative"), each = length(sigma_env_coo))
)

# 计算均值用于画虚线
#means <- df_env %>%
#  group_by(group) %>%
#  summarise(mean_value = mean(value))

# 绘图
ggplot(df_env, aes(x = value, color = group, fill = group)) +
  geom_density(alpha = 0.3, linewidth = 1,bw = 0.1) +
  # geom_vline(data = means, aes(xintercept = mean_value, color = group),
  #linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = c("tomato", "steelblue")) +
  scale_fill_manual(values = c("tomato", "steelblue")) +
  labs(x = expression(sigma[env]),
       y = "Density",
       title = expression("Posterior distribution of " * sigma[env])) +
  coord_cartesian(xlim = c(0, 1)) +  # 限制 x 轴范围
  theme_classic(base_size = 14)

#种群随机增长率的估计绘图
df_a <- data.frame(
  value = c(a_coo, a_non),
  group = rep(c("Cooperative", "Non-cooperative"), each = length(a_coo))
)
means <- df_a %>%
  group_by(group) %>%
  summarise(mean_value = mean(value))

# 绘图
ggplot(df_a, aes(x = value, color = group, fill = group)) +
  geom_density(alpha = 0.3, linewidth = 1, bw = 0.15) +  # 可微调 bw 控制平滑度
  geom_vline(data = means, aes(xintercept = mean_value, color = group),
             linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = c("tomato", "steelblue")) +
  scale_fill_manual(values = c("tomato", "steelblue")) +
  labs(x = expression(ln(a)[sto]),
       y = "Density",
       title = expression("Posterior distribution of environmental stochastic growth rate " * ln(a)[sto])) +
  coord_cartesian(xlim = range(df_a$value)) +  # 自动适应范围，也可手动指定
  theme_classic(base_size = 14)

# 整理数据框
df_r <- data.frame(
  value = c(r_coo, r_non),
  group = rep(c("Cooperative", "Non-cooperative"), each = length(r_coo))
)

# 计算均值用于画虚线（可选）
means_r <- df_r %>%
  group_by(group) %>%
  summarise(mean_value = mean(value))

# 绘图
ggplot(df_r, aes(x = value, color = group, fill = group)) +
  geom_density(alpha = 0.3, linewidth = 1, bw = 0.1) +
  geom_vline(data = means_r, aes(xintercept = mean_value, color = group),
             linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = c("tomato", "steelblue")) +
  scale_fill_manual(values = c("tomato", "steelblue")) +
  labs(x = expression(r),
       y = "Density",
       title = expression("Posterior distribution of ln(mean) " * r)) +
  coord_cartesian(xlim = range(df_r$value)) +
  theme_classic(base_size = 14)
