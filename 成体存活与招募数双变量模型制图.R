#成体存活与招募数的双变量模型绘图
#读入模型（模型制作见“RCode_bivariate_model.R”）
#model
library(MCMCglmm)
model_SurRecr<-readRDS( "model_SurRecr.rds") # data from：RCode_bivariate_model.R
summary(model_SurRecr)
head(a_merged)

#进行模型预测值计算（首先整理数据结构，然后提取模型后验系数，然后进行计算）
# Step 1：扩展 a_merged 为长格式，加上 trait 列
n <- nrow(a_merged)
data_long <- rbind(a_merged, a_merged)
data_long$trait <- rep(c("traitRecruits1", "traitSurvival"), each = n)

# Step 2：构建 design matrix（使用模型的固定效应结构）
X <- model.matrix(model_SurRecr$Fixed$formula, data = data_long)

# Step 3：提取后验固定效应系数（n_iter × n_coef）
beta_samples <- model_SurRecr$Sol

# Step 4：进行预测（线性预测值）
linpred_mat <- X %*% t(beta_samples)  # 2n × n_iter 的矩阵

# Step 5：计算后验均值和置信区间（logit/log尺度）
linpred_mean <- rowMeans(linpred_mat)
linpred_CI   <- apply(linpred_mat, 1, quantile, probs = c(0.025, 0.975))

# Step 6：转换为原始尺度
inv_logit <- function(x) exp(x) / (1 + exp(x))
mu_mean  <- c(
  sapply(1:n, function(i) mean(inv_logit(linpred_mat[i, ]))),            # Survival
  sapply((n+1):(2*n), function(i) mean(exp(linpred_mat[i, ])))           # Recruits1
)

mu_lower <- c(
  sapply(1:n, function(i) quantile(inv_logit(linpred_mat[i, ]), 0.025)),
  sapply((n+1):(2*n), function(i) quantile(exp(linpred_mat[i, ]), 0.025))
)

mu_upper <- c(
  sapply(1:n, function(i) quantile(inv_logit(linpred_mat[i, ]), 0.975)),
  sapply((n+1):(2*n), function(i) quantile(exp(linpred_mat[i, ]), 0.975))
)

# Step 7：整合结果（风格一致）
prediction_df <- data.frame(
  id        = rep(1:n, 2),
  trait     = rep(c("Survival", "Recruits1"), each = n),
  mu_mean   = mu_mean,
  mu_lower  = mu_lower,
  mu_upper  = mu_upper
)

# Step 8：拆分预测值，合并回宽格式
pred_surv   <- prediction_df[prediction_df$trait == "Survival", ]
pred_recr   <- prediction_df[prediction_df$trait == "Recruits1", ]

a_merged_pred <- a_merged
a_merged_pred$Surv_pred  <- pred_surv$mu_mean
a_merged_pred$Recr_pred  <- pred_recr$mu_mean

# （可选）也可以加上 CI
a_merged_pred$Surv_lower <- pred_surv$mu_lower
a_merged_pred$Surv_upper <- pred_surv$mu_upper
a_merged_pred$Recr_lower <- pred_recr$mu_lower
a_merged_pred$Recr_upper <- pred_recr$mu_upper

# 查看结果
head(a_merged_pred) 
write.csv(a_merged_pred, file = "a_merged_pred.csv", row.names = FALSE) 


library(dplyr)
library(ggplot2)

# Step 1: 按 TempNov.z 均值分类（高温 / 低温）
temp_mean <- mean(a_merged_pred$TempNov.z, na.rm = TRUE)

plot_df <- a_merged_pred %>%
  mutate(
    cooperation = factor(Gethp, levels = c(1, 0), labels = c("Cooperative", "Non-cooperative")),
    temp_group = ifelse(TempNov.z > temp_mean, "High Temp", "Low Temp"),
    temp_group = factor(temp_group, levels = c("Low Temp", "High Temp"))
  )

# Step 2: 自定义颜色（合作 vs 非合作）
color_values <- c(
  "Cooperative" = "#F05D60",       # 红色
  "Non-cooperative" = "#1C71B6"    # 蓝色
)

# Step 3: 自定义线型（冷冬 vs 暖冬）
linetype_values <- c(
  "Low Temp" = "solid",
  "High Temp" = "dashed"
)

# Step 4: 绘图
ggplot(plot_df, aes(x = PreAM.z, y = Surv_pred)) +
  
  # 点图
  geom_point(aes(shape = cooperation,
                 size = cooperation,
                 color = cooperation),
             alpha = 0.5) +
  
  # 回归线
  geom_smooth(aes(color = cooperation,
                  linetype = temp_group),
              method = "lm", se = TRUE, linewidth = 1) +
  
  # 手动设置样式
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  scale_shape_manual(values = c("Cooperative" = 16, "Non-cooperative" = 17)) +
  scale_size_manual(values = c("Cooperative" = 1.8, "Non-cooperative" = 1.5)) +
  
  # 标签与主题
  labs(
    x = "Standardized Precipitation (April–May)",
    y = "Predicted Adult Survival",
    color = "Cooperation",
    linetype = "Temperature",
    shape = "Cooperation",
    size = "Cooperation"
  ) +
  theme_classic(base_size = 13) +  # 干净坐标系
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold")
  )



library(dplyr)
library(ggplot2)

# Step 1: 按 TempNov.z 均值分类（高温 / 低温）
temp_mean <- mean(a_merged_pred$TempNov.z, na.rm = TRUE)

plot_df <- a_merged_pred %>%
  mutate(
    cooperation = factor(Gethp, levels = c(1, 0),
                         labels = c("Cooperative", "Non-cooperative")),
    temp_group = ifelse(TempNov.z > temp_mean, "High Temp", "Low Temp"),
    temp_group = factor(temp_group, levels = c("Low Temp", "High Temp"))
  )

# Step 2: 自定义颜色（合作 vs 非合作）
color_values <- c(
  "Cooperative" = "#F05D60",       # 红色
  "Non-cooperative" = "#1C71B6"    # 蓝色
)

# Step 3: 自定义线型（冷冬 vs 暖冬）
linetype_values <- c(
  "Low Temp" = "solid",
  "High Temp" = "dashed"
)

# Step 4: 绘图
ggplot(plot_df, aes(x = PreAM.z, y = Recr_pred)) +
  
  # 点图（只用颜色，不要 linetype）
  geom_point(aes(shape = cooperation,
                 size = cooperation,
                 color = cooperation),
             alpha = 0.6) +
  
  # 回归线（颜色 = 合作，线型 = 冷/暖冬）
  geom_smooth(aes(color = cooperation,
                  linetype = temp_group),
              method = "lm", se = TRUE, linewidth = 1) +
  
  # 手动设置样式
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  scale_shape_manual(values = c("Cooperative" = 16,
                                "Non-cooperative" = 17)) +
  scale_size_manual(values = c("Cooperative" = 1.8,
                               "Non-cooperative" = 1.5)) +
  
  # 标签与主题
  labs(
    x = "Standardized Precipitation (April–May)",
    y = "Predicted Number of Recruits",
    color = "Cooperation",
    linetype = "Temperature",
    shape = "Cooperation",
    size = "Cooperation"
  ) +
  theme_classic(base_size = 13) +  # 干净坐标系
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold")
  )


#数据检查

# 查看 Recr_pred 的最大值、最小值和分布
summary(a_merged_pred$Recr_pred)
hist(a_merged_pred$Recr_pred, breaks = 50, main = "Histogram of Recr_pred")


library(ggplot2)
library(tidyr)
library(dplyr)
# 转成长格式便于绘图
beta_long <- as.data.frame(beta_samples) %>%
  mutate(iteration = 1:n()) %>%
  pivot_longer(cols = -iteration, names_to = "parameter", values_to = "value")

# 画所有参数的密度图
ggplot(beta_long, aes(x = value)) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  facet_wrap(~ parameter, scales = "free") +
  theme_minimal() +
  labs(title = "Posterior Distributions of Model Parameters")

library(dplyr)
library(ggplot2)

# ----------------------------------------------------
# Figure 5b 
# Gethp (合作状态 1=合作, 0=非合作), TempNov.z (温度Z分), Surv_pred (预测存活率)
# ----------------------------------------------------

plot_df <- a_merged_pred %>%
  mutate(
    cooperation = factor(Gethp, levels = c(1, 0), labels = c("Cooperative", "Non-cooperative"))
  )

# 自定义颜色（合作 vs 非合作）
color_values <- c(
  "Cooperative" = "#F05D60",      # 红色
  "Non-cooperative" = "#1C71B6"    # 蓝色
)

# 绘图：将 TempNov.z 作为 X 轴梯度
ggplot(plot_df, aes(x = TempNov.z, y = Surv_pred)) +
  
  # 点图 (保持不变，但 X 轴已切换为 TempNov.z)
  geom_point(aes(shape = cooperation,
                 size = cooperation,
                 color = cooperation),
             alpha = 0.5) +
  
  # 回归线：只按 cooperation 分组，使用实线（默认）
  geom_smooth(aes(color = cooperation),
              method = "lm", 
              se = TRUE, 
              linewidth = 1,
              linetype = "solid") + # 明确指定实线
  
  # 手动设置样式
  scale_color_manual(values = color_values) +
  scale_shape_manual(values = c("Cooperative" = 16, "Non-cooperative" = 17)) +
  scale_size_manual(values = c("Cooperative" = 1.8, "Non-cooperative" = 1.5)) +
  
  # 标签与主题：X 轴改为温度
  labs(
    x = "Standardized Temperature (November)", # 根据你的变量 TempNov.z 调整标签
    y = "Predicted Adult Survival",
    color = "Cooperation",
    shape = "Cooperation",
    size = "Cooperation"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold")
  )
