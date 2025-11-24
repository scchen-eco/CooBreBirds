#################################################幼体存活与出飞窝雏数的双变量模型绘图####################################################
library(ggplot2)
#模型和数据见代码：“双变量模型出飞数和存活巢水平.R”
#model
model_OffSurFle_Nest<-readRDS("model_OffSurFle_Nest.rds") #data：nest_level_data
summary(model_OffSurFle_Nest)
head(nest_level_data)
# 构建数据结构，并进行数据预测
# Step 1：扩展 nest_level_data 为长格式，加上 trait 列
n <- nrow(nest_level_data)
data_long <- rbind(nest_level_data, nest_level_data)
data_long$trait <- rep(c("traitbrood", "traitSUR_total"), each = n)

# Step 2：构建 design matrix（使用模型的固定效应结构）
X <- model.matrix(model_OffSurFle_Nest$Fixed$formula, data = data_long)

# Step 3：提取后验固定效应系数（n_iter × n_coef）19200个后验系数（sample size）
beta_samples <- model_OffSurFle_Nest$Sol

# Step 4：进行预测（线性预测值）
linpred_mat <- X %*% t(beta_samples)  # 得到 2n × n_iter 的预测矩阵

# Step 5：计算后验均值和置信区间（logμ尺度）
linpred_mean <- rowMeans(linpred_mat)
linpred_CI   <- apply(linpred_mat, 1, quantile, probs = c(0.025, 0.975))

# Step 6：转换为原始尺度（exp，因为 Poisson 模型）
mu_mean  <- exp(linpred_mean)
mu_lower <- exp(linpred_CI[1, ])
mu_upper <- exp(linpred_CI[2, ])

# Step 7：整合结果
prediction_df <- data.frame(
  id        = rep(1:n, 2),
  trait     = rep(c("SUR_total", "brood"), each = n),
  mu_mean   = mu_mean,
  mu_lower  = mu_lower,
  mu_upper  = mu_upper
)
# Step 8：拆分预测值，合并回宽格式
pred_SUR   <- prediction_df[prediction_df$trait == "SUR_total", ]
pred_brood <- prediction_df[prediction_df$trait == "brood", ]

nest_level_data_pred <- nest_level_data
nest_level_data_pred$SUR_pred   <- pred_SUR$mu_mean
nest_level_data_pred$brood_pred <- pred_brood$mu_mean
head(nest_level_data_pred)


###########绘图
library(dplyr)
library(ggplot2)

# 添加高/低降水标签
precip_mean <- mean(nest_level_data_pred$total_precip_apr_may, na.rm = TRUE)

plot_df <- nest_level_data_pred %>%
  mutate(
    cooperation = factor(Gethp, levels = c(1, 0), labels = c("Cooperative", "Non-cooperative")),
    precip_group = ifelse(total_precip_apr_may > precip_mean, "High Precip", "Low Precip"),
    precip_group = factor(precip_group, levels = c("Low Precip", "High Precip")),
    color_group = paste(precip_group, cooperation, sep = "_")
  )

# 自定义颜色（深色合作，浅色非合作）
color_values <- c(
  "Low Precip_Cooperative" = "#1f77b4",        # 深蓝
  "High Precip_Cooperative" = "#ff7f0e",       # 深橙
  "Low Precip_Non-cooperative" = "#aec7e8",    # 浅蓝
  "High Precip_Non-cooperative" = "#ffbb78"    # 浅橙
)

# 绘图
ggplot(plot_df, aes(x = avg_temp_nov, y = SUR_pred)) +
  
  # 点图
  geom_point(aes(shape = cooperation,
                 size = cooperation,
                 color = color_group),
             alpha = 0.6) +
  
  # 回归线：分别画合作与非合作两组
  geom_smooth(data = filter(plot_df, cooperation == "Cooperative"),
              aes(color = color_group, linetype = cooperation),
              method = "lm", se = TRUE, linewidth = 1.0) +
  
  geom_smooth(data = filter(plot_df, cooperation == "Non-cooperative"),
              aes(color = color_group, linetype = cooperation),
              method = "lm", se = TRUE, linewidth = 0.8) +
  
  # 手动设置样式
  scale_color_manual(
    values = color_values,
    labels = c(
      "Low Precip_Cooperative" = "Low Precip (Coop.)",
      "High Precip_Cooperative" = "High Precip (Coop.)",
      "Low Precip_Non-cooperative" = "Low Precip (Non-coop.)",
      "High Precip_Non-cooperative" = "High Precip (Non-coop.)"
    )
  ) +
  scale_shape_manual(values = c("Cooperative" = 16, "Non-cooperative" = 17)) +
  scale_linetype_manual(values = c("Cooperative" = "solid", "Non-cooperative" = "dashed")) +
  scale_size_manual(values = c("Cooperative" = 1.8, "Non-cooperative" = 1.5)) +
  
  # 标签与主题
  labs(
    x = "Average November Temperature",
    y = "Predicted Juvenile Survival",
    color = "Group",
    shape = "Cooperation",
    linetype = "Cooperation",
    size = "Cooperation"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

library(dplyr)
library(ggplot2)

# Step 1: 按 avg_temp_nov 均值分类（高温 / 低温）
temp_mean <- mean(nest_level_data_pred$avg_temp_nov, na.rm = TRUE)

plot_df <- nest_level_data_pred %>%
  mutate(
    cooperation = factor(Gethp, levels = c(1, 0),
                         labels = c("Cooperative", "Non-cooperative")),
    temp_group = ifelse(avg_temp_nov > temp_mean, "High Temp", "Low Temp"),
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
ggplot(plot_df, aes(x = total_precip_apr_may, y = SUR_pred)) +
  
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
    x = "Total Precipitation (April–May)",
    y = "Predicted Juvenile Survival",
    color = "Cooperation",
    linetype = "Temperature",
    shape = "Cooperation",
    size = "Cooperation"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold")
  )


library(dplyr)
library(ggplot2)

# Step 1: 创建高/低温分组 + 合作组标签 + 颜色分组
# Step 1: 按 avg_temp_nov 均值分类（高温 / 低温）
temp_mean <- mean(nest_level_data_pred$avg_temp_nov, na.rm = TRUE)

plot_df <- nest_level_data_pred %>%
  mutate(
    cooperation = factor(Gethp, levels = c(1, 0),
                         labels = c("Cooperative", "Non-cooperative")),
    temp_group = ifelse(avg_temp_nov > temp_mean, "High Temp", "Low Temp"),
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

# Step 4: 绘图（y轴 = brood_pred）
ggplot(plot_df, aes(x = total_precip_apr_may, y = brood_pred)) +
  
  # 点图（颜色 = 合作/非合作）
  geom_point(aes(shape = cooperation,
                 size = cooperation,
                 color = cooperation),
             alpha = 0.6) +
  
  # 回归线（颜色 = 合作/非合作，线型 = 冷冬/暖冬）
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
    x = "Total Precipitation (April–May)",
    y = "Predicted Fledglings",
    color = "Cooperation",
    linetype = "Temperature",
    shape = "Cooperation",
    size = "Cooperation"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold")
  )



#figure4a 4b

library(dplyr)
library(coda)

# 假设模型和数据对象已加载
model_OffSurFle_Nest_helped <- readRDS("model_OffSurFle_Nest_helped.rds")

# --- 合作组数据准备 ---
model_Coo <- model_OffSurFle_Nest_helped
data_Coo <- data_helped
n_Coo <- nrow(data_Coo)

# Step 1：扩展数据为长格式 (假设 trait 1 是 brood, trait 2 是 SUR)
data_long_Coo <- rbind(data_Coo, data_Coo)
data_long_Coo$trait <- rep(c("traitbrood", "traitSUR_total"), each = n_Coo)

# Step 2：构建 design matrix
X_Coo <- model.matrix(model_Coo$Fixed$formula, data = data_long_Coo)

# Step 3：提取后验固定效应系数
beta_samples_Coo <- model_Coo$Sol

# Step 4：进行预测（线性预测值）
linpred_mat_Coo <- X_Coo %*% t(beta_samples_Coo)

# Step 5：计算后验均值和置信区间（logμ尺度）
linpred_mean_Coo <- rowMeans(linpred_mat_Coo)
linpred_CI_Coo   <- apply(linpred_mat_Coo, 1, quantile, probs = c(0.025, 0.975))

# Step 6：转换为原始尺度（exp）
mu_mean_Coo  <- exp(linpred_mean_Coo)
mu_lower_Coo <- exp(linpred_CI_Coo[1, ])
mu_upper_Coo <- exp(linpred_CI_Coo[2, ])

# Step 7：整合结果
prediction_df_Coo <- data.frame(
  id      = rep(1:n_Coo, 2),
  trait   = rep(c("SUR_total", "brood"), each = n_Coo),
  mu_mean  = mu_mean_Coo,
  mu_lower = mu_lower_Coo,
  mu_upper = mu_upper_Coo
)

# Step 8：提取 brood size 的预测值，并合并回原始数据框
pred_brood_Coo <- prediction_df_Coo %>% filter(trait == "brood")

data_helped_pred <- data_Coo
data_helped_pred$brood_pred_mean  <- pred_brood_Coo$mu_mean
data_helped_pred$brood_pred_lower <- pred_brood_Coo$mu_lower
data_helped_pred$brood_pred_upper <- pred_brood_Coo$mu_upper
pred_SUR_Coo <- prediction_df_Coo %>% 
  filter(trait == "SUR_total")
data_helped_pred$SUR_pred_mean  <- pred_SUR_Coo$mu_mean
data_helped_pred$SUR_pred_lower <- pred_SUR_Coo$mu_lower
data_helped_pred$SUR_pred_upper <- pred_SUR_Coo$mu_upper

head(data_helped_pred) # Figure 4a 的绘图数据准备完毕

# 假设模型和数据对象已加载
model_OffSurFle_Nest_Non <- readRDS("model_OffSurFle_Nest_Non.rds")
# data_nonhelped

# --- 非合作组数据准备 ---
model_NonCoo <- model_OffSurFle_Nest_Non
data_NonCoo <- data_nonhelped
n_NonCoo <- nrow(data_NonCoo)

# Step 1：扩展数据为长格式
data_long_NonCoo <- rbind(data_NonCoo, data_NonCoo)
data_long_NonCoo$trait <- rep(c("brood", "SUR_total"), each = n_NonCoo)

# Step 2：构建 design matrix
X_NonCoo <- model.matrix(model_NonCoo$Fixed$formula, data = data_long_NonCoo)

# Step 3：提取后验固定效应系数
beta_samples_NonCoo <- model_NonCoo$Sol

# Step 4：进行预测（线性预测值）
linpred_mat_NonCoo <- X_NonCoo %*% t(beta_samples_NonCoo)

# Step 5：计算后验均值和置信区间（logμ尺度）
linpred_mean_NonCoo <- rowMeans(linpred_mat_NonCoo)
linpred_CI_NonCoo   <- apply(linpred_mat_NonCoo, 1, quantile, probs = c(0.025, 0.975))

# Step 6：转换为原始尺度（exp）
mu_mean_NonCoo  <- exp(linpred_mean_NonCoo)
mu_lower_NonCoo <- exp(linpred_CI_NonCoo[1, ])
mu_upper_NonCoo <- exp(linpred_CI_NonCoo[2, ])

# Step 7：整合结果
prediction_df_NonCoo <- data.frame(
  id      = rep(1:n_NonCoo, 2),
  trait   = rep(c("SUR_total", "brood"), each = n_NonCoo),
  mu_mean  = mu_mean_NonCoo,
  mu_lower = mu_lower_NonCoo,
  mu_upper = mu_upper_NonCoo
)

# Step 8：提取 brood size 的预测值，并合并回原始数据框
pred_brood_NonCoo <- prediction_df_NonCoo %>% filter(trait == "brood")

data_nonhelped_pred <- data_NonCoo
data_nonhelped_pred$brood_pred_mean  <- pred_brood_NonCoo$mu_mean
data_nonhelped_pred$brood_pred_lower <- pred_brood_NonCoo$mu_lower
data_nonhelped_pred$brood_pred_upper <- pred_brood_NonCoo$mu_upper
pred_SUR_NonCoo <- prediction_df_NonCoo %>% 
  filter(trait == "SUR_total")
data_nonhelped_pred$SUR_pred_mean  <- pred_SUR_NonCoo$mu_mean
data_nonhelped_pred$SUR_pred_lower <- pred_SUR_NonCoo$mu_lower
data_nonhelped_pred$SUR_pred_upper <- pred_SUR_NonCoo$mu_upper

head(data_nonhelped_pred) # Figure 4b 的绘图数据准备完毕




library(dplyr)
library(ggplot2)

# 定义颜色 (仅使用合作组的红色)
color_cooperative <- "#F05D60"

# Step 1: 准备合作组数据（将 X 和 Y 轴变量名标准化）
# 移除 Gethp 和 temp_group 的不必要操作
plot_df_4a <- data_helped_pred %>%
  rename(
    brood_pred = brood_pred_mean, # 使用预测均值作为 Y 轴
    total_precip_apr_may = total_precip_apr_may
  )

# Step 2: 绘图（y轴 = brood_pred）
ggplot(plot_df_4a, aes(x = total_precip_apr_may, y = brood_pred)) +
  
  # 点图 (固定颜色和形状)
  geom_point(color = color_cooperative, 
             shape = 16, 
             size = 2, 
             alpha = 0.6) +
  
  geom_smooth(method = "lm", 
              se = TRUE, 
              linewidth = 1,
              color = color_cooperative,
              fill = color_cooperative, # 使用颜色填充 CI 区域
              alpha = 0.2) + # 降低 CI 填充的透明度
  
  # 标签与主题
  labs(
    x = "Spring Rainfall (mm)",
    y = "Predicted Brood Size (Estimated Mean)",
    title = "Figure 4a: Cooperative Phenotypes"
  ) +
  theme_classic(base_size = 13) +
  theme(
    # 移除所有关于图例的设置，因为图中只有一个趋势
    legend.position = "none",
    axis.title = element_text(face = "bold")
  )

library(dplyr)
library(ggplot2)
library(patchwork)

# 定义颜色 
color_cooperative <- "#F05D60"     # 红色
color_noncooperative <- "#1C71B6"  # 蓝色

# --- Figure 4a: 合作组 (Cooperative) ---

plot_df_4a <- data_helped_pred %>%
  rename(
    brood_pred = brood_pred_mean,
    total_precip_apr_may = total_precip_apr_may
  )

plot_4a <- ggplot(plot_df_4a, aes(x = total_precip_apr_may, y = brood_pred)) +
  
  # 1. 回归线和 CI
  geom_smooth(method = "lm", 
              se = TRUE, 
              linewidth = 1,
              color = color_cooperative,
              fill = color_cooperative,
              alpha = 0.2) +
  
  # 2. 散点图
  geom_point(color = color_cooperative, 
             shape = 16, 
             size = 2, 
             alpha = 0.6) +
  
  # 3. 统一 Y 轴范围
  coord_cartesian(ylim = c(3, 7)) + 
  
  # 4. 标签与主题
  labs(
    x = "Spring Rainfall (mm)",
    y = "Predicted Brood Size (Estimated Mean)",
    title = "Figure 4a: Cooperative Phenotypes"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", axis.title = element_text(face = "bold"))


# --- Figure 4b: 非合作组 (Non-cooperative) ---

plot_df_4b <- data_nonhelped_pred %>%
  rename(
    brood_pred = brood_pred_mean,
    total_precip_apr_may = total_precip_apr_may
  )

plot_4b <- ggplot(plot_df_4b, aes(x = total_precip_apr_may, y = brood_pred)) +
  
  # 1. 回归线和 CI
  geom_smooth(method = "lm", 
              se = TRUE, 
              linewidth = 1,
              color = color_noncooperative,
              fill = color_noncooperative,
              alpha = 0.2) +
  
  # 2. 散点图
  geom_point(color = color_noncooperative, 
             shape = 17, 
             size = 2, 
             alpha = 0.6) +
  
  # 3. 统一 Y 轴范围
  coord_cartesian(ylim = c(3, 7)) + 
  
  # 4. 标签与主题
  labs(
    x = "Spring Rainfall (mm)",
    y = "Predicted Brood Size (Estimated Mean)",
    title = "Figure 4b: Non-cooperative Phenotypes"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", axis.title = element_text(face = "bold"))


# 组合图表并显示
Combined_Figure_4 <- plot_4a + plot_4b
Combined_Figure_4

library(dplyr)
library(ggplot2)

# 假设 nest_level_data_pred 是你提供的包含 SUR_pred 的数据框

# --- Step 1: 准备绘图数据 (只创建合作分组) ---
plot_df <- nest_level_data_pred %>%
  # 确保 SUR_pred 是 Y 轴的预测值
  mutate(
    cooperation = factor(Gethp, levels = c(1, 0),
                         labels = c("Cooperative", "Non-cooperative"))
    # avg_temp_nov 是 X 轴
  )

# --- Step 2: 自定义颜色（合作 vs 非合作）---
color_values <- c(
  "Cooperative" = "#F05D60",    # 红色
  "Non-cooperative" = "#1C71B6" # 蓝色
)

# --- Step 3: 绘图 (X轴 = 温度，Y轴 = 存活率) ---
ggplot(plot_df, aes(x = avg_temp_nov, y = SUR_pred, color = cooperation)) +
  
  # 点图 (颜色和形状根据合作/非合作分组)
  geom_point(aes(shape = cooperation, size = cooperation),
             alpha = 0.6) +
  
  # 回归线 (颜色根据 cooperation 分组，使用 lm 方法和 CI)
  geom_smooth(method = "lm", 
              se = TRUE, 
              linewidth = 1,
              # fill 继承 color 的值，确保 CI 填充色与线色一致
              aes(fill = cooperation)) + 
  
  # 手动设置样式
  scale_color_manual(values = color_values) +
  scale_fill_manual(values = color_values) + # 确保 CI 填充也使用自定义颜色
  scale_shape_manual(values = c("Cooperative" = 16, "Non-cooperative" = 17)) +
  scale_size_manual(values = c("Cooperative" = 1.8, "Non-cooperative" = 1.5)) +
  
  # 标签与主题 (修正 X/Y 轴标签)
  labs(
    x = "Average November Temperature (z-score)",
    y = "Predicted Juvenile Survival",
    title = "Predicted Survival vs. November Temperature by Cooperation Status",
    color = "Cooperation",
    shape = "Cooperation",
    size = "Cooperation"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold")
  )