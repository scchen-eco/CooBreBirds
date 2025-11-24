#####################################################################################################
#data_from Negative_Binomial_zipoisson_Ksai.R
modelKsai_nb               <- readRDS("modelKsai_nb.rds") #data:a_merged_ExH_numeric 
modelKsai_nb_NonCoo        <- readRDS("modelKsai_nb_NonCoo.rds") #data:a_merged_ExH_numeric_Non 
modelKsai_nb_Coo           <- readRDS("modelKsai_nb_Coo.rds") #data:a_merged_ExH_numeric_Coo 
#整合数据
a_merged_ExH_numeric_pred <- a_merged_ExH_numeric
a_merged_ExH_numeric_pred$ksai_1_pred<- predict(modelKsai_nb, newdata = a_merged_ExH_numeric, type = "response")
head(a_merged_ExH_numeric_pred)
#非合作数据
a_merged_ExH_numeric_Non_pred <- a_merged_ExH_numeric_Non 
a_merged_ExH_numeric_Non_pred$ksai_1_pred<- predict(modelKsai_nb_NonCoo, newdata = a_merged_ExH_numeric_Non, type = "response")
head(a_merged_ExH_numeric_Non_pred)
#合作数据
a_merged_ExH_numeric_Coo_pred <- a_merged_ExH_numeric_Coo  
a_merged_ExH_numeric_Coo_pred$ksai_1_pred<- predict(modelKsai_nb_Coo, newdata = a_merged_ExH_numeric_Coo , type = "response")
head(a_merged_ExH_numeric_Coo_pred)
library(ggplot2)
library(dplyr) 

# 将 Gethp 转为标签
a_merged_ExH_numeric_pred <- a_merged_ExH_numeric_pred %>%
  mutate(Cooperation = ifelse(Gethp == 1, "Cooperative", "Non-cooperative"))

# 绘图
ggplot(a_merged_ExH_numeric_pred, aes(x = environment, y = ksai_1_pred, fill = Cooperation)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, 
               width = 0.7, alpha = 0.3,
               position = position_dodge(0.8)) +
  labs(
    x = "Environment category",
    y = "Predicted ksai",
    fill = "Group"
  ) +
  scale_fill_manual(values = c("Cooperative" = "#FF5757", "Non-cooperative" = "#0070C0")) +
  theme_classic(base_size = 13) +   # 经典主题：只保留XY轴
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top"
  )


library(ggplot2)
library(dplyr)
# 非合作
# temp_category 确保为 factor，控制顺序（可选）
a_merged_ExH_numeric_Non_pred <- a_merged_ExH_numeric_Non_pred %>%
  mutate(temp_category = factor(temp_category, levels = c("Low Temp", "High Temp")))

# 绘图
ggplot(a_merged_ExH_numeric_Non_pred, aes(x = total_precip_apr_may, y = ksai_1_pred, color = temp_category)) +
  geom_point(alpha = 0.3, size = 1.2) +  # 可选散点层
  geom_smooth(method = "lm", se = TRUE, linewidth =  1.2) +
  labs(
    x = "Total precipitation (April–May)",
    y = "Predicted ksai (non-cooperative)",
    color = "Temperature category"
  ) +
  scale_color_manual(values = c("Low Temp" = "#1f77b4", "High Temp" = "#ff7f0e")) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# 合作
# temp_category 确保为 factor，控制顺序（可选）
a_merged_ExH_numeric_Coo_pred <- a_merged_ExH_numeric_Coo_pred %>%
  mutate(temp_category = factor(temp_category, levels = c("Low Temp", "High Temp")))

# 绘图
ggplot(a_merged_ExH_numeric_Coo_pred, aes(x = total_precip_apr_may, y = ksai_1_pred, color = temp_category)) +
  geom_point(alpha = 0.3, size = 1.2) +  # 可选散点层
  geom_smooth(method = "lm", se = TRUE, linewidth =  1.2) +
  labs(
    x = "Total precipitation (April–May)",
    y = "Predicted ksai (cooperative)",
    color = "Temperature category"
  ) +
  scale_color_manual(values = c("Low Temp" = "#1f77b4", "High Temp" = "#ff7f0e")) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# 合作与非合作
# 添加合作标签
coo_df <- a_merged_ExH_numeric_Coo_pred %>%
  mutate(cooperation = "Cooperative")

non_df <- a_merged_ExH_numeric_Non_pred %>%
  mutate(cooperation = "Non-cooperative")

# 合并并整理变量类型
combined_df <- bind_rows(coo_df, non_df) %>%
  mutate(
    cooperation = factor(cooperation, levels = c("Cooperative", "Non-cooperative")),
    temp_category = factor(temp_category, levels = c("Low Temp", "High Temp"))
  )

# 自定义颜色：合作为深色，非合作为浅色
color_values <- c(
  "Low Temp_Cooperative" = "#1f77b4",      # 深蓝
  "High Temp_Cooperative" = "#ff7f0e",     # 深橙/红
  "Low Temp_Non-cooperative" = "#aec7e8",  # 浅蓝
  "High Temp_Non-cooperative" = "#ffbb78"  # 浅橙
)

# 添加一个变量用于颜色映射
combined_df$color_group <- paste(combined_df$temp_category, combined_df$cooperation, sep = "_")

# 绘图
ggplot(combined_df, aes(x = total_precip_apr_may, y = ksai_1_pred)) +
  
  # 点图（颜色由 color_group 控制）
  geom_point(aes(shape = cooperation,
                 size = cooperation,
                 color = color_group)) +
  
  # 回归线（按合作组分别添加，颜色由 color_group 控制）
  geom_smooth(data = filter(combined_df, cooperation == "Cooperative"),
              aes(color = color_group, linetype = cooperation),
              method = "lm", se = TRUE, linewidth = 1.0) +
  
  geom_smooth(data = filter(combined_df, cooperation == "Non-cooperative"),
              aes(color = color_group, linetype = cooperation),
              method = "lm", se = TRUE, linewidth = 0.8) +
  
  # 手动缩放设置
  scale_color_manual(values = color_values,
                     labels = c(
                       "Low Temp_Cooperative" = "Low Temp (Coop.)",
                       "High Temp_Cooperative" = "High Temp (Coop.)",
                       "Low Temp_Non-cooperative" = "Low Temp (Non-coop.)",
                       "High Temp_Non-cooperative" = "High Temp (Non-coop.)"
                     )) +
  scale_shape_manual(values = c("Cooperative" = 16, "Non-cooperative" = 17)) +
  scale_linetype_manual(values = c("Cooperative" = "solid", "Non-cooperative" = "dashed")) +
  scale_size_manual(values = c("Cooperative" = 1.8, "Non-cooperative" = 1.5)) +
  
  # 图例与主题
  labs(
    x = "Total precipitation (April–May)",
    y = "Predicted ksai",
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



#Alternatively
# 自定义颜色：红/蓝代表合作/非合作
color_values <- c(
  "Cooperative" = "#F05D60",       # 红色
  "Non-cooperative" = "#1C71B6"    # 蓝色
)

# 自定义线型：冷冬/暖冬
linetype_values <- c(
  "Low Temp" = "solid",
  "High Temp" = "dashed"
)

ggplot(combined_df, aes(x = total_precip_apr_may, y = ksai_1_pred)) +
  
  # 点图
  geom_point(aes(shape = cooperation,
                 size = cooperation,
                 color = cooperation,
                 linetype = temp_category)) +
  
  # 回归线
  geom_smooth(aes(color = cooperation,
                  linetype = temp_category),
              method = "lm", se = TRUE, linewidth = 1) +
  
  # 手动缩放设置
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  scale_shape_manual(values = c("Cooperative" = 16, "Non-cooperative" = 17)) +
  scale_size_manual(values = c("Cooperative" = 1.8, "Non-cooperative" = 1.5)) +
  
  # 图例与主题
  labs(
    x = "Total precipitation (April–May)",
    y = "Predicted ksai",
    color = "Cooperation",
    linetype = "Temperature",
    shape = "Cooperation",
    size = "Cooperation"
  ) +
  theme_classic(base_size = 13) +   # 经典主题：只保留XY轴
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top"
  )







































