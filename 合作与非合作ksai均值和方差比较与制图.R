library(ggplot2)
library(tidyr)
library(dplyr)
library(car) 
#制图 合作与非合作的ksai的对比，箱线图
individual_contribution <- read.csv("dx_regression_Pt.csv") 
head(individual_contribution)
# 去掉 state 为 "H" 的行
individual_contribution <- individual_contribution %>%
  filter(state != "H")
#不区分年进行比较
t.test(ksai ~ Gethp, data = individual_contribution, var.equal = FALSE) # Welch's t 检验 方差不等

#区分年进行配对的t检验
# 每年每个Gethp组（0/1）计算 ksai 平均值
yearly_mean <- individual_contribution %>%
  group_by(year, Gethp) %>%
  summarise(annual_ksai = mean(ksai, na.rm = TRUE), .groups = "drop")
#统计合作与非合作年ksai的均值和SE
# 计算多年内每组（合作/非合作）年均 ksai 的均值与标准误
# 统计合作与非合作年ksai的均值、SE 和方差
# 计算多年内每组（合作/非合作）年均 ksai 的均值、标准差、方差和标准误
summary_stats <- yearly_mean %>%
  group_by(Gethp) %>%
  summarise(
    mean_annual_ksai = mean(annual_ksai, na.rm = TRUE),
    sd = sd(annual_ksai, na.rm = TRUE),
    n = n(),
    SE = sd / sqrt(n),
    CV = sd / mean_annual_ksai,
    .groups = "drop"
  ) %>%
  # 添加年间方差 (Variance = sd^2)
  mutate(
    Variance = sd^2,
    Group = ifelse(Gethp == 1, "Cooperative", "Non-cooperative")
  ) %>%
  # 重新排列和选择列，将 Variance 包含在内
  select(Group, mean_annual_ksai, sd, Variance, SE, CV, n)

print(summary_stats)

paired_data <- yearly_mean %>%
  select(year, Gethp, annual_ksai) %>%
  tidyr::pivot_wider(names_from = Gethp, values_from = annual_ksai,
                     names_prefix = "Gethp_") %>%
  filter(!is.na(Gethp_0) & !is.na(Gethp_1))  # 确保两组都有数据

t.test(paired_data$Gethp_1, paired_data$Gethp_0, paired = TRUE) # Paired t-test
# 提取合作组 (Gethp = 1) 和非合作组 (Gethp = 0) 的年均 ksai
ksai_cooperative <- yearly_mean %>%
  filter(Gethp == 1) %>%
  pull(annual_ksai)

ksai_noncooperative <- yearly_mean %>%
  filter(Gethp == 0) %>%
  pull(annual_ksai)

# 1. 使用 F 检验比较方差
# var.test(x, y) 用于检验两个样本的方差是否相等 (原假设 H0: var(x) = var(y))
variance_test_result <- var.test(ksai_cooperative, ksai_noncooperative)

# 打印 F 检验结果
print("--- F 检验 (比较合作与非合作年间 ksai 年均值的方差) ---")
print(variance_test_result)
print("---------------------------------------------------")

# ----------------------------------------------------
# 比较合作与非合作巢 ksai 个体间方差的同质性
# 原假设 (H0): 合作与非合作巢的 ksai 方差相等
# ----------------------------------------------------
levene_test_result <- car::leveneTest(ksai ~ factor(Gethp), 
                                      data = individual_contribution, 
                                      center = median) # 使用中位数 (median) 更稳健

print("--- Levene 检验 (比较合作与非合作巢 ksai 年内个体间方差) ---")
print(levene_test_result)
print("-----------------------------------------------------------------")


# 补充: 计算两个组 ksai 个体间的实际方差值 (Variance)
individual_variance <- individual_contribution %>%
  group_by(Gethp) %>%
  summarise(
    n_individuals = n(),
    Individual_Variance = var(ksai, na.rm = TRUE),
    Individual_SD = sd(ksai, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Group = ifelse(Gethp == 1, "Cooperative", "Non-cooperative"))

print("--- 合作与非合作巢 ksai 个体间方差的具体值 ---")
print(individual_variance)
print("--------------------------------------------------")



#可视化
# 转成长格式
plot_data <- paired_data %>%
  pivot_longer(cols = c(Gethp_0, Gethp_1),
               names_to = "condition", values_to = "ksai") %>%
  mutate(
    condition = ifelse(condition == "Gethp_0", "Noncooperative", "Cooperative"),
    condition = factor(condition, levels = c("Noncooperative", "Cooperative"))
  )

# 画图
ggplot(plot_data, aes(x = condition, y = ksai)) +
  geom_boxplot(aes(fill = condition), outlier.shape = NA, width = 0.5, alpha = 0.3) +
  geom_point(aes(color = condition, group = year), size = 2) +
  geom_line(aes(group = year), color = "gray", linewidth = 0.4) +
  scale_fill_manual(values = c("Noncooperative" = "#0070C0", "Cooperative" = "#FF5757")) +
  scale_color_manual(values = c("Noncooperative" = "#0070C0", "Cooperative" = "#FF5757")) +
  labs(x = "", y = "Yearly mean ksai") +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.x = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(0.4, "cm")
  )

#点图
ggplot(summary_stats, aes(x = Group, y = mean_annual_ksai, fill = Group)) +
  # 绘制点图，将其颜色 (color) 也映射到 Group
  geom_point(aes(color = Group), size = 4) + 
  
  # 绘制误差棒，将其颜色 (color) 也映射到 Group
  geom_errorbar(aes(ymin = mean_annual_ksai - sd, ymax = mean_annual_ksai + sd, color = Group), 
                width = 0.2, linewidth = 1) +
  
  # 保持 scale_fill_manual（虽然在这里可能没有填充，但习惯性保留，或者如果你改成柱状图会用到）
  scale_fill_manual(values = c("Non-cooperative" = "#0070C0", "Cooperative" = "#FF5757")) +
  
  # 现在 scale_color_manual 会起作用了，因为它有对应的美学映射
  scale_color_manual(values = c("Non-cooperative" = "#0070C0", "Cooperative" = "#FF5757")) +
  
  labs(x = "", y = expression("Mean Annual "*xi), 
       title = "Figure 1b: Between-Year Variance Comparison (SD)") +
  theme_classic() +
  theme(legend.position = "none")

# Figure 1b: 年间均值和方差比较 (强调方差的柱状图)

ggplot(summary_stats, aes(x = Group, y = mean_annual_ksai)) + # <--- 移除 fill=Group
  
  # 1. 绘制柱状图：无填充，银灰色边框
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           width = 0.5,
           fill = "white",       # <--- 柱子内部填充白色
           color = "grey70") +  # <--- 柱子边框为银灰色
  
  # 2. 绘制彩色点
  geom_point(aes(color = Group), # <--- color 映射到 Group
             size = 4) +
  
  # 3. 绘制彩色误差棒 (使用 SD 代表年间变异性)
  geom_errorbar(aes(ymin = mean_annual_ksai - sd, 
                    ymax = mean_annual_ksai + sd, 
                    color = Group), # <--- color 映射到 Group
                position = position_dodge(width = 0.5), 
                width = 0.15, linewidth = 1) +
  
  # 4. 设置彩色点和误差棒的颜色
  scale_color_manual(values = c("Non-cooperative" = "#0070C0", "Cooperative" = "#FF5757")) +
  
  # 5. （可选）如果你不想有图例，保留 scale_fill_manual 也没关系，因为没有 fill 美学映射
  # scale_fill_manual(values = c("Non-cooperative" = "#0070C0", "Cooperative" = "#FF5757")) +
  
  # 6. 设置标签和主题
  labs(x = "", y = expression("Mean Annual "*xi), 
       title = "Figure 1b: Between-Year Variance Comparison (SD)") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.x = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(0.4, "cm")
  )

# Figure 1c: 年内个体间方差比较 (小提琴图)

# 确保 individual_contribution 包含原始 ksai 和 Gethp
individual_plot_data <- individual_contribution %>%
  mutate(Group = ifelse(Gethp == 1, "Cooperative", "Non-cooperative"),
         Group = factor(Group, levels = c("Non-cooperative", "Cooperative")))

ggplot(individual_plot_data, aes(x = Group, y = ksai, fill = Group)) +
  # 1. 绘制小提琴图：用宽度展示数据密度，即方差
  geom_violin(trim = TRUE, alpha = 0.5) + 
  
  # 2. 叠加箱线图：用线条展示中位数和四分位数范围
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.8) + 
  
  # 3. 设置填充颜色
  scale_fill_manual(values = c("Non-cooperative" = "#0070C0", "Cooperative" = "#FF5757")) +
  
  # 4. 设置标签和主题
  labs(x = "", y = expression("Individual "*xi), 
       title = "Figure 1c: Within-Year (Individual) Variance Comparison") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.x = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(0.4, "cm")
  )

