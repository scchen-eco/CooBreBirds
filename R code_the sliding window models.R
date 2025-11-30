##################滑动窗口模型####################
library(climwin)
library(tidyverse)
library(lme4)
library(lmerTest)
##################繁殖力##################
climbio<-read.csv("NESTlEVEL.csv") 
#climbio <- climbio %>%mutate(Experience1 = ifelse(Experience1 >= 1, 1, Experience1))%>% 
#  filter(brood !=0,!is.na(Experience1) & !is.na(Terr.)) #过滤掉完全繁殖失败的巢,未知繁殖经验和栖息地质量的巢
climbio <- climbio %>%mutate(Experience1 = ifelse(Experience1 >= 1, 1, Experience1))%>% 
  filter(brood !=0) #过滤掉完全繁殖失败的巢,未知繁殖经验和栖息地质量的巢
head(climbio)
dim(climbio)
table(climbio$Gethp)
#降水对繁殖力
climlong<-read.csv("cleaned_dataday_pre.csv")
Broodsize_Pre<- slidingwin(xvar = list(Pre= climlong$Pre),
                           cdate = climlong$Date,
                           bdate = climbio$Date,
                           baseline = glm(brood ~  Gethp + scale(Terr.) +scale(Experience1), data = climbio, family = poisson),
                           cinterval = "month",
                           range = c(6, 0),
                           type = "absolute", refday = c(31, 07),
                           stat = "mean",
                           func = "lin")
#随机化
Broodsize_RandPre <- randwin(repeats = 1000,
                             xvar = list(Pre= climlong$Pre),
                             cdate = climlong$Date,
                             bdate = climbio$Date,
                             baseline = glm(brood ~  Gethp+scale(Terr.) +scale(Experience1), data = climbio, family = poisson),
                             cinterval = "month",
                             range = c(6, 0),
                             type = "absolute", refday = c(31, 07),
                             stat = "mean",
                             func = "lin")

#温度对繁殖力
climlong<-read.csv("cleaned_climlongday.csv")
Broodsize_Tmp<- slidingwin(xvar = list(Tmp= climlong$tmpvalue),
                           cdate = climlong$Date,
                           bdate = climbio$Date,
                           baseline = glm(brood ~  Gethp + scale(Terr.) +scale(Experience1), data = climbio, family = poisson),
                           cinterval = "month",
                           range = c(6, 0),
                           type = "absolute", refday = c(31, 07),
                           stat = "mean",
                           func = "lin")
#随机化
Broodsize_RandTmp <- randwin(repeats = 1000,
                             xvar = list(Tmp= climlong$tmpvalue),
                             cdate = climlong$Date,
                             bdate = climbio$Date,
                             baseline = glm(brood ~  Gethp+scale(Terr.) +scale(Experience1), data = climbio, family = poisson),
                             cinterval = "month",
                             range = c(6, 0),
                             type = "absolute", refday = c(31, 07),
                             stat = "mean",
                             func = "lin")

#################合作繁殖对繁殖力
#降水
climlong<-read.csv("cleaned_dataday_pre.csv")
climbio$climate <- 1
Broodsize_Pre_Coo <- slidingwin(xvar = list(Pre= climlong$Pre),
                                cdate = climlong$Date,
                                bdate = climbio$Date,
                                baseline = glm(brood ~  Gethp + scale(Terr.) +scale(Experience1)+Gethp:climate, data = climbio, family = poisson),
                                cinterval = "month",
                                range = c(6, 0),
                                type = "absolute", refday = c(31, 07),
                                stat = "mean",
                                func = "lin")
#随机化
Broodsize_Randpre_Coo <- randwin_modified(repeats = 1000,
                                          xvar = list(Pre= climlong$Pre),
                                          cdate = climlong$Date,
                                          bdate = climbio$Date,
                                          baseline = glm(brood ~  Gethp + scale(Terr.) +scale(Experience1)+Gethp:climate, data = climbio, family = poisson),
                                          cinterval = "month",
                                          range = c(6, 0),
                                          type = "absolute", refday = c(31, 07),
                                          stat = "mean",
                                          func = "lin") #该随机化过程报错，需要使用循环方法记录

#寻找错误原因
results_list <- list()  # 存储每次重复的结果
for (i in 1:10) {
  # 使用 tryCatch 捕获每次重复的错误
  result <- tryCatch({
    randwin(repeats = 1, 
            xvar = list(Pre = climlong$Pre),
            cdate = climlong$Date,
            bdate = climbio$Date,
            baseline = glm(brood ~ Gethp + scale(Terr.) +scale(Experience1)+ Gethp:climate, 
                           data = climbio, family = poisson),
            cinterval = "month",
            range = c(6, 0),
            type = "absolute",
            refday = c(31, 07),
            stat = "mean",
            func = "lin")
  }, error = function(e) {
    # 打印错误并返回 NA
    message("Error in randomization number ", i, ": ", e)
    return(NULL)  # 可以用 NULL 标记出错
  })
  
  # 将每次成功运行的结果添加到结果列表中
  if (!is.null(result)) {
    results_list[[i]] <- result
  }
}

# 检查所有结果的列名是否一致
for (i in seq_along(results_list)) {
  print(paste("Randomization number", i, "structure:"))
  print(str(results_list[[i]]))  # 检查每次随机化的结果
}

final_result <- do.call(rbind, results_list)# 合并结果
final_result #成功捕获错误原因，在于随机化数据后交互项不收敛。

#考虑到采样成功的概率，进行大量采样，保存每一次采样结果
results_list <- list()  
# 列表中的所有可能列名
all_columns <- c("deltaAICc", "WindowOpen", "WindowClose", "Function", "Furthest", 
                 "Closest", "Statistics", "Type", "K", "ModWeight", "sample.size",
                 "Reference.day", "Reference.month", "Custom.mod", "Gethp", "Terr.", 
                 "Experience1", "Gethp:climate", "GethpSE", "Terr.SE", "Experience1SE", 
                 "Gethp:climateSE", "Randomised", "Repeat", "WeightDist")

for (i in 1:12000) {
  # 使用 tryCatch 捕获每次重复的错误
  result <- tryCatch({
    randwin(repeats = 1, 
            xvar = list(Pre = climlong$Pre),
            cdate = climlong$Date,
            bdate = climbio$Date,
            baseline = glm(brood ~ Gethp + scale(Terr.) +scale(Experience1) + Gethp:climate, 
                           data = climbio, family = poisson),
            cinterval = "month",
            range = c(6, 0),
            type = "absolute",
            refday = c(31, 07),
            stat = "mean",
            func = "lin")
  }, error = function(e) {
    message("Error in randomization number ", i, ": ", e)
    return(NULL)  # 可以用 NULL 标记出错
  })
  
  # 确保结果存在
  if (!is.null(result)) {
    # 手动添加缺失的列并填充NA
    missing_cols <- setdiff(all_columns, names(result[[1]]))
    for (col in missing_cols) {
      result[[1]][[col]] <- NA
    }
    results_list[[i]] <- result[[1]]  # 只保存第一部分的数据
  }
}
# 合并结果
# 保留所有结果，包括失败的，统一为同样结构
results_list_fixed <- lapply(results_list, function(x) {
  if (is.null(x)) {
    # 构造全 NA 的一行数据框，保留列名
    return(as.data.frame(matrix(NA, nrow = 1, ncol = length(all_columns),
                                dimnames = list(NULL, all_columns))))
  }
  # 填补缺失列
  missing_cols <- setdiff(all_columns, names(x))
  for (col in missing_cols) {
    x[[col]] <- NA
  }
  # 按照统一顺序排列
  x <- x[all_columns]
  return(x)
})
# 合并为一个数据框
final_result <- bind_rows(results_list_fixed)
#去掉那些采样失败的滑窗模型
cleaned_result <- final_result %>%
  filter(!is.na(`Gethp:climate`)) 
Broodsize_RandPre_Coo <- cleaned_result
dim(Broodsize_RandPre_Coo) #检查有效数据量

################合作繁殖对繁殖力
#温度
climlong<-read.csv("cleaned_climlongday.csv")
Broodsize_Tmp_Coo<- slidingwin(xvar = list(Tmp= climlong$tmpvalue),
                               cdate = climlong$Date,
                               bdate = climbio$Date,
                               baseline = glm(brood ~  Gethp + scale(Terr.) +scale(Experience1)+Gethp:climate, data = climbio, family = poisson),
                               cinterval = "month",
                               range = c(6, 0),
                               type = "absolute", refday = c(31, 07),
                               stat = "mean",
                               func = "lin")
#随机化
results_list <- list()  
# 列表中的所有可能列名
all_columns <- c("deltaAICc", "WindowOpen", "WindowClose", "Function", "Furthest", 
                 "Closest", "Statistics", "Type", "K", "ModWeight", "sample.size",
                 "Reference.day", "Reference.month", "Custom.mod", "Gethp", "Terr.", 
                 "Experience1", "Gethp:climate", "GethpSE", "Terr.SE", "Experience1SE", 
                 "Gethp:climateSE", "Randomised", "Repeat", "WeightDist")

for (i in 1:12000) {
  # 使用 tryCatch 捕获每次重复的错误 # 此处3000次循环即可
  result <- tryCatch({
    randwin(repeats = 1, 
            xvar = list(Tmp= climlong$tmpvalue),
            cdate = climlong$Date,
            bdate = climbio$Date,
            baseline = glm(brood ~ Gethp + scale(Terr.) +scale(Experience1) + Gethp:climate, 
                           data = climbio, family = poisson),
            cinterval = "month",
            range = c(6, 0),
            type = "absolute",
            refday = c(31, 07),
            stat = "mean",
            func = "lin")
  }, error = function(e) {
    message("Error in randomization number ", i, ": ", e)
    return(NULL)  # 可以用 NULL 标记出错
  })
  
  # 确保结果存在
  if (!is.null(result)) {
    # 手动添加缺失的列并填充NA
    missing_cols <- setdiff(all_columns, names(result[[1]]))
    for (col in missing_cols) {
      result[[1]][[col]] <- NA
    }
    results_list[[i]] <- result[[1]]  # 只保存第一部分的数据
  }
}
# 合并结果
# 保留所有结果，包括失败的，统一为同样结构
results_list_fixed <- lapply(results_list, function(x) {
  if (is.null(x)) {
    # 构造全 NA 的一行数据框，保留列名
    return(as.data.frame(matrix(NA, nrow = 1, ncol = length(all_columns),
                                dimnames = list(NULL, all_columns))))
  }
  # 填补缺失列
  missing_cols <- setdiff(all_columns, names(x))
  for (col in missing_cols) {
    x[[col]] <- NA
  }
  # 按照统一顺序排列
  x <- x[all_columns]
  return(x)
})
# 合并为一个数据框
final_result <- bind_rows(results_list_fixed)
dim(final_result)
#去掉那些采样失败的滑窗模型
cleaned_result <- final_result %>%
  filter(!is.na(`Gethp:climate`)) 
Broodsize_RandTmp_Coo <- cleaned_result
dim(Broodsize_RandTmp_Coo) #检查有效数据量


#################幼体存活窗口#################
climbio<-read.csv("DX_offss.csv")  #这已经去掉了完全繁殖失败的巢，但是没有约束区域
climbio <- climbio %>%mutate(Exp. = ifelse(Exp. >= 1, 1, Exp.),Gethp = as.numeric(Gethp))%>% 
  filter(Sex !="x")#将性别不明的个体去掉
climbio $Exp.[is.na(climbio$Exp.)] <- median(climbio $Exp., na.rm = TRUE)
climlong<-read.csv("cleaned_climlongday.csv")
table(climbio$Gethp)

#后代存活与气候滑窗_温度
Offs_SurWintmp<- slidingwin(xvar = list(Tmp= climlong$tmpvalue),
                            cdate = climlong$Date,
                            bdate = climbio$Date,
                            baseline = glm(SUR  ~  Gethp + scale(Terr.) +Sex + scale(Exp.) , data = climbio,family = binomial),
                            cinterval = "month",
                            range = c(10, 0),
                            type = "absolute", refday = c(30, 05),
                            stat = "mean",
                            func = "lin")
#随机化模型
Offs_SurRandtmp <- randwin(repeats = 1000, 
                           xvar = list(Tmp= climlong$tmpvalue),
                           cdate = climlong$Date,
                           bdate = climbio$Date,
                           baseline = glm(SUR  ~  Gethp + scale(Terr.) +Sex + scale(Exp.) , data = climbio,family = binomial),
                           cinterval = "month",
                           range = c(10, 0),
                           type = "absolute", refday = c(30, 05),
                           stat = "mean",
                           func = "lin")

#合作繁殖后代存活效应滑窗法_温度
climbio$climate <- 1
Offs_SurWintmp_Coo<- slidingwin(xvar = list(Tmp= climlong$tmpvalue),
                                cdate = climlong$Date,
                                bdate = climbio$Date,
                                baseline = glm(SUR  ~  Gethp +scale(Terr.) +Sex + scale(Exp.)+Gethp:climate , data = climbio,family = binomial),
                                cinterval = "month",
                                range = c(10, 0),
                                type = "absolute", refday = c(30, 05),
                                stat = "mean",
                                func = "lin")
##随机化模型
library(plyr)
basewin <- climwin:::basewin
Offs_SurRandtmp_Coo <- randwin_modified(repeats = 12000, 
                                        xvar = list(Tmp= climlong$tmpvalue),
                                        cdate = climlong$Date,
                                        bdate = climbio$Date,
                                        baseline = glm(SUR  ~  Gethp +scale(Terr.) +Sex + scale(Exp.)+Gethp:climate , data = climbio,family = binomial),
                                        cinterval = "month",
                                        range = c(10, 0),
                                        type = "absolute", refday = c(30, 05),
                                        stat = "mean",
                                        func = "lin")#1000次随机化中只有1个有效的，因此加大重复次数。

#后代存活与气候滑窗_降水
climlong<-read.csv("cleaned_dataday_pre.csv")
Offs_SurWinpre<- slidingwin(xvar = list(Pre= climlong$Pre),
                            cdate = climlong$Date,
                            bdate = climbio$Date,
                            baseline = glm(SUR  ~  Gethp + scale(Terr.) +Sex + scale(Exp.) , data = climbio,family = binomial),
                            cinterval = "month",
                            range = c(15, 0),
                            type = "absolute", refday = c(30, 05),
                            stat = "mean",
                            func = "lin")
##随机化模型
Offs_SurRandpre <- randwin_modified(repeats = 1000, 
                                    xvar = list(Pre= climlong$Pre),
                                    cdate = climlong$Date,
                                    bdate = climbio$Date,
                                    baseline = glm(SUR  ~  Gethp + scale(Terr.) +Sex + scale(Exp.) , data = climbio,family = binomial),
                                    cinterval = "month",
                                    range = c(15, 0),
                                    type = "absolute", refday = c(30, 05),
                                    stat = "mean",
                                    func = "lin") 

#合作繁殖后代存活效应滑窗_降水
climbio$climate <- 1
Offs_SurWinpre_Coo<- slidingwin(xvar = list(Pre= climlong$Pre),
                                cdate = climlong$Date,
                                bdate = climbio$Date,
                                baseline = glm(SUR  ~  Gethp +scale(Terr.) +Sex + scale(Exp.)+Gethp:climate , data = climbio,family = binomial),
                                cinterval = "month",
                                range = c(15, 0),
                                type = "absolute", refday = c(30, 05),
                                stat = "mean",
                                func = "lin")
#随机化
Offs_SurRandpre_Coo <- randwin(repeats = 5, 
                               xvar = list(Pre= climlong$Pre),
                               cdate = climlong$Date,
                               bdate = climbio$Date,
                               baseline = glm(SUR  ~  Gethp +scale(Terr.) +Sex + scale(Exp.)+Gethp:climate , data = climbio,family = binomial),
                               cinterval = "month",
                               range = c(15, 0),
                               type = "absolute", refday = c(30, 05),
                               stat = "mean",
                               func = "lin") #随机化过程中交互项被弃，不成功。
# 获取两组列名
cols1 <- colnames(Offs_SurWinpre_Coo[[1]]$Dataset)
cols2 <- colnames(Offs_SurRandpre_Coo[[1]])
#考虑到采样成功的概率，进行大量采样，保存每一次采样结果
results_list <- list()  
# 列表中的所有可能列名
all_columns <- union(cols1, cols2)

for (i in 1:12000) {
  # 使用 tryCatch 捕获每次重复的错误 
  result <- tryCatch({
    randwin(repeats = 1, 
            xvar = list(Pre= climlong$Pre),
            cdate = climlong$Date,
            bdate = climbio$Date,
            baseline = glm(SUR  ~  Gethp +scale(Terr.) +Sex + scale(Exp.)+Gethp:climate , data = climbio,family = binomial),
            cinterval = "month",
            range = c(15, 0),
            type = "absolute", refday = c(30, 05),
            stat = "mean",
            func = "lin")
  }, error = function(e) {
    message("Error in randomization number ", i, ": ", e)
    return(NULL)  # 可以用 NULL 标记出错
  })
  
  # 确保结果存在
  if (!is.null(result)) {
    # 手动添加缺失的列并填充NA
    missing_cols <- setdiff(all_columns, names(result[[1]]))
    for (col in missing_cols) {
      result[[1]][[col]] <- NA
    }
    results_list[[i]] <- result[[1]]  # 只保存第一部分的数据
  }
}
# 合并结果

# 保留所有结果，包括失败的，统一为同样结构
results_list_fixed <- lapply(results_list, function(x) {
  if (is.null(x)) {
    # 构造全 NA 的一行数据框，保留列名
    return(as.data.frame(matrix(NA, nrow = 1, ncol = length(all_columns),
                                dimnames = list(NULL, all_columns))))
  }
  # 填补缺失列
  missing_cols <- setdiff(all_columns, names(x))
  for (col in missing_cols) {
    x[[col]] <- NA
  }
  # 按照统一顺序排列
  x <- x[all_columns]
  return(x)
})
# 合并为一个数据框
final_result <- bind_rows(results_list_fixed)
dim(final_result)
#去掉那些采样失败的滑窗模型
cleaned_result <- final_result %>%
  filter(!is.na(`Gethp:climate`)) 
Offs_SurRandpre_Coo <- cleaned_result
dim(Offs_SurRandpre_Coo) #检查有效数据量


#################成体存活窗口
climbio<-read.csv("P_SUR_Age_TQ_Ex.csv")
climbio <- climbio %>% filter(Status != "H")
climbio$year <- as.character(climbio$Year)  # 确保年份格式一致
climbio $Experience1[is.na(climbio$Experience1)] <- median(climbio $Experience1, na.rm = TRUE)
climbio $TQ3[is.na(climbio$TQ3)] <- median(climbio $TQ3, na.rm = TRUE)
climbio <- climbio %>%mutate(Experience1 = ifelse(Experience1 >= 1, 1, Experience1),Recieved = as.numeric(Recieved))
dim(climbio)
head(climbio)
table(climbio$Recieved)
##################成体存活气候窗口
##温度
climlong<-read.csv("cleaned_climlongday.csv")
Adult_SurWintmp<- slidingwin(xvar = list(Tmp= climlong$tmpvalue),
                             cdate = climlong$Date,
                             bdate = climbio$Date,
                             baseline = glm(SUR  ~  Recieved + scale(TQ3) + Sex + scale(Experience1), data = climbio,family = binomial),
                             cinterval = "month",
                             range = c(10, 0),
                             type = "absolute", refday = c(30, 05),
                             stat = "mean",
                             func = "lin")
##随机化模型
Adult_SurRandtmp <- randwin_modified(repeats = 1000, 
                                     xvar = list(Tmp= climlong$tmpvalue),
                                     cdate = climlong$Date,
                                     bdate = climbio$Date,
                                     baseline = glm(SUR  ~  Recieved + scale(TQ3) + Sex + scale(Experience1), data = climbio,family = binomial),
                                     cinterval = "month",
                                     range = c(10, 0),
                                     type = "absolute", refday = c(30, 05),
                                     stat = "mean",
                                     func = "lin") 

#####################成体存活合作繁殖效应窗口
#温度
climbio$climate <- 1
Adult_SurWintmp_Coo<- slidingwin(xvar = list(Tmp= climlong$tmpvalue),
                                 cdate = climlong$Date,
                                 bdate = climbio$Date,
                                 baseline = glm(SUR  ~  Recieved + scale(TQ3) + Sex + scale(Experience1) + Recieved:climate, data = climbio,family = binomial),
                                 cinterval = "month",
                                 range = c(10, 0),
                                 type = "absolute", refday = c(30, 05),
                                 stat = "mean",
                                 func = "lin")
#随机化
Adult_SurRandtmp_Coo<- randwin(repeats = 5, 
                               xvar = list(Tmp= climlong$tmpvalue),
                               cdate = climlong$Date,
                               bdate = climbio$Date,
                               baseline = glm(SUR  ~  Recieved + scale(TQ3) + Sex + scale(Experience1) + Recieved:climate, data = climbio,family = binomial),
                               cinterval = "month",
                               range = c(10, 0),
                               type = "absolute", refday = c(30, 05),
                               stat = "mean",
                               func = "lin")  #没有估计交互项

# 获取两组列名
cols1 <- colnames(Adult_SurWintmp_Coo[[1]]$Dataset)
cols2 <- colnames(Adult_SurRandtmp_Coo[[1]])

#考虑到采样成功的概率，进行大量采样，保存每一次采样结果
results_list <- list()  
# 列表中的所有可能列名
all_columns <- union(cols1, cols2)

for (i in 1:12000) {
  # 使用 tryCatch 捕获每次重复的错误 
  result <- tryCatch({
    randwin(repeats = 1, 
            xvar = list(Tmp= climlong$tmpvalue),
            cdate = climlong$Date,
            bdate = climbio$Date,
            baseline = glm(SUR  ~  Recieved + scale(TQ3) + Sex + scale(Experience1) + Recieved:climate, data = climbio,family = binomial),
            cinterval = "month",
            range = c(10, 0),
            type = "absolute", refday = c(30, 05),
            stat = "mean",
            func = "lin")
  }, error = function(e) {
    message("Error in randomization number ", i, ": ", e)
    return(NULL)  # 可以用 NULL 标记出错
  })
  
  # 确保结果存在
  if (!is.null(result)) {
    # 手动添加缺失的列并填充NA
    missing_cols <- setdiff(all_columns, names(result[[1]]))
    for (col in missing_cols) {
      result[[1]][[col]] <- NA
    }
    results_list[[i]] <- result[[1]]  # 只保存第一部分的数据
  }
}

# 保留所有结果，包括失败的，统一为同样结构
results_list_fixed <- lapply(results_list, function(x) {
  if (is.null(x)) {
    # 构造全 NA 的一行数据框，保留列名
    return(as.data.frame(matrix(NA, nrow = 1, ncol = length(all_columns),
                                dimnames = list(NULL, all_columns))))
  }
  # 填补缺失列
  missing_cols <- setdiff(all_columns, names(x))
  for (col in missing_cols) {
    x[[col]] <- NA
  }
  # 按照统一顺序排列
  x <- x[all_columns]
  return(x)
})
# 合并为一个数据框
final_result <- bind_rows(results_list_fixed)
dim(final_result)
head(final_result)
#去掉那些采样失败的滑窗模型
cleaned_result <- final_result %>%
  filter(!is.na(`Recieved:climate`)) 
Adult_SurRandtmp_Coo<- cleaned_result
dim(Adult_SurRandtmp_Coo) #检查有效数据量1494个

#################成体存活窗口
climbio<-read.csv("P_SUR_Age_TQ_Ex.csv")
climbio <- climbio %>% filter(Status != "H")
climbio$year <- as.character(climbio$Year)  # 确保年份格式一致
climbio $Experience1[is.na(climbio$Experience1)] <- median(climbio $Experience1, na.rm = TRUE)
climbio $TQ3[is.na(climbio$TQ3)] <- median(climbio $TQ3, na.rm = TRUE)
climbio <- climbio %>%mutate(Experience1 = ifelse(Experience1 >= 1, 1, Experience1),Recieved = as.numeric(Recieved))
##################成体存活气候窗口
##降水
climlong<-read.csv("cleaned_dataday_pre.csv")
Adult_SurWinpre<- slidingwin(xvar = list(Pre= climlong$Pre),
                             cdate = climlong$Date,
                             bdate = climbio$Date,
                             baseline = glm(SUR  ~  Recieved + scale(TQ3) + Sex + scale(Experience1), data = climbio,family = binomial),
                             cinterval = "month",
                             range = c(15, 0),
                             type = "absolute", refday = c(30, 05),
                             stat = "mean",
                             func = "lin")
##随机化模型
Adult_SurRandpre <- randwin_modified(repeats = 1000, 
                                     xvar = list(Pre= climlong$Pre),
                                     cdate = climlong$Date,
                                     bdate = climbio$Date,
                                     baseline = glm(SUR  ~  Recieved + scale(TQ3) + Sex + scale(Experience1), data = climbio,family = binomial),
                                     cinterval = "month",
                                     range = c(15, 0),
                                     type = "absolute", refday = c(30, 05),
                                     stat = "mean",
                                     func = "lin") 

#####################成体存活合作繁殖效应窗口
#降水
climbio$climate <- 1
Adult_SurWinpre_Coo<- slidingwin(xvar = list(Pre= climlong$Pre),
                                 cdate = climlong$Date,
                                 bdate = climbio$Date,
                                 baseline = glm(SUR  ~  Recieved + scale(TQ3) + Sex + scale(Experience1) + Recieved:climate, data = climbio,family = binomial),
                                 cinterval = "month",
                                 range = c(15, 0),
                                 type = "absolute", refday = c(30, 05),
                                 stat = "mean",
                                 func = "lin")

#随机化
Adult_SurRandpre_Coo<- randwin_modified(repeats = 2, 
                                        xvar = list(Pre= climlong$Pre),
                                        cdate = climlong$Date,
                                        bdate = climbio$Date,
                                        baseline = glm(SUR  ~  Recieved + scale(TQ3) + Sex + scale(Experience1) + Recieved:climate, data = climbio,family = binomial),
                                        cinterval = "month",
                                        range = c(15, 0),
                                        type = "absolute", refday = c(30, 05),
                                        stat = "mean",
                                        func = "lin")  #没有估计交互项

# 获取两组列名
cols1 <- colnames(Adult_SurWinpre_Coo[[1]]$Dataset)
cols2 <- colnames(Adult_SurRandpre_Coo[[1]])

#考虑到采样成功的概率，进行大量采样，保存每一次采样结果
results_list <- list()  
# 列表中的所有可能列名
all_columns <- union(cols1, cols2)

for (i in 1:12000) {
  # 使用 tryCatch 捕获每次重复的错误 
  result <- tryCatch({
    randwin_modified(repeats = 1, 
                     xvar = list(Pre= climlong$Pre),
                     cdate = climlong$Date,
                     bdate = climbio$Date,
                     baseline = glm(SUR  ~  Recieved + scale(TQ3) + Sex + scale(Experience1) + Recieved:climate, data = climbio,family = binomial),
                     cinterval = "month",
                     range = c(15, 0),
                     type = "absolute", refday = c(30, 05),
                     stat = "mean",
                     func = "lin")
  }, error = function(e) {
    message("Error in randomization number ", i, ": ", e)
    return(NULL)  # 可以用 NULL 标记出错
  })
  
  # 确保结果存在
  if (!is.null(result)) {
    # 手动添加缺失的列并填充NA
    missing_cols <- setdiff(all_columns, names(result[[1]]))
    for (col in missing_cols) {
      result[[1]][[col]] <- NA
    }
    results_list[[i]] <- result[[1]]  # 只保存第一部分的数据
  }
}

# 保留所有结果，包括失败的，统一为同样结构
results_list_fixed <- lapply(results_list, function(x) {
  if (is.null(x)) {
    # 构造全 NA 的一行数据框，保留列名
    return(as.data.frame(matrix(NA, nrow = 1, ncol = length(all_columns),
                                dimnames = list(NULL, all_columns))))
  }
  # 填补缺失列
  missing_cols <- setdiff(all_columns, names(x))
  for (col in missing_cols) {
    x[[col]] <- NA
  }
  # 按照统一顺序排列
  x <- x[all_columns]
  return(x)
})
# 合并为一个数据框
final_result <- bind_rows(results_list_fixed)
dim(final_result)
head(final_result)
#去掉那些采样失败的滑窗模型
cleaned_result <- final_result %>%
  filter(!is.na(`Recieved:climate`)) 
Adult_SurRandpre_Coo<- cleaned_result
dim(Adult_SurRandpre_Coo) #有效数据量为0。
# saveRDS(final_result, file = "Adult_SurRandpre_Coo_final_result_list.rds")

####可视化模型的窗口
files_rds <- c(
  "Adult_SurRandpre.rds",
  "Adult_SurRandpre_Coo.rds",
  "Adult_SurRandpre_Coo_final_result_list.rds",
  "Adult_SurRandtmp.rds",
  "Adult_SurRandtmp_Coo.rds",
  "Adult_SurWinpre.rds",
  "Adult_SurWinpre_Coo.rds",
  "Adult_SurWintmp.rds",
  "Adult_SurWintmp_Coo.rds",
  "Broodsize_Pre.rds",
  "Broodsize_Pre_Coo.rds",
  "Broodsize_RandPre.rds",
  "Broodsize_RandPre_Coo.rds",
  "Broodsize_RandTmp.rds",
  "Broodsize_RandTmp_Coo.rds",
  "Broodsize_Tmp.rds",
  "Broodsize_Tmp_Coo.rds",
  "Offs_SurRandpre.rds",
  "Offs_SurRandpre_Coo.rds",
  "Offs_SurRandtmp.rds",
  "Offs_SurRandtmp_Coo.rds",
  "Offs_SurWinpre.rds",
  "Offs_SurWinpre_Coo.rds",
  "Offs_SurWintmp.rds",
  "Offs_SurWintmp_Coo.rds"
)
#读取其中之一进行查看
Adult_SurWintmp_Coo <- readRDS("Adult_SurWintmp_Coo.rds")

head(Adult_SurWintmp_Coo[[1]]$Dataset) #查看模型排序
medwin(Adult_SurWintmp_Coo[[1]]$Dataset) #总结所有可能的模型
BroodsizeOutput<-Adult_SurWintmp_Coo[[1]]$Dataset
#可视化窗口
plotdelta(dataset = BroodsizeOutput)
plotweights(dataset = BroodsizeOutput)
plotwin(dataset = BroodsizeOutput)


Offs_SurWintmp_Coo <- readRDS("Offs_SurWintmp_Coo.rds")

head(Offs_SurWintmp_Coo[[1]]$Dataset) #查看模型排序
medwin(Offs_SurWintmp_Coo[[1]]$Dataset) #总结所有可能的模型
BroodsizeOutput<-Offs_SurWintmp_Coo[[1]]$Dataset
#可视化窗口
plotdelta(dataset = BroodsizeOutput)
plotweights(dataset = BroodsizeOutput)
plotwin(dataset = BroodsizeOutput)

Broodsize_Pre_Coo <- readRDS("Broodsize_Pre_Coo.rds")

head(Broodsize_Pre_Coo[[1]]$Dataset) #查看模型排序
medwin(Broodsize_Pre_Coo[[1]]$Dataset) #总结所有可能的模型
BroodsizeOutput<-Broodsize_Pre_Coo[[1]]$Dataset
#可视化窗口
plotdelta(dataset = BroodsizeOutput)
plotweights(dataset = BroodsizeOutput)
plotwin(dataset = BroodsizeOutput)








