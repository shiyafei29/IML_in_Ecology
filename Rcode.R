
#读入数据
data <- read.csv("New_Alpine_Meadow_Plant_Diversity_Data.csv")
View(data)

#####################探索性数据分析-20250504##############################
# 计算描述性统计结果
summary_stats <- data.frame(
  Mean = sapply(data, mean, na.rm = TRUE),
  Median = sapply(data, median, na.rm = TRUE),
  SD = sapply(data, sd, na.rm = TRUE),
  Min = sapply(data, min, na.rm = TRUE),
  Max = sapply(data, max, na.rm = TRUE),
  N = sapply(data, function(x) sum(!is.na(x)))
)
print(round(summary_stats, 2))


# 加载必要包
library(corrplot)
library(ggplot2)

# 读取数据（请根据实际路径修改）
data <- read.csv("New_Alpine_Meadow_Plant_Diversity_Data.csv")
# 计算皮尔森相关性矩阵
cor_matrix <- cor(data, use = "complete.obs", method = "pearson")
# 打印相关性矩阵（保留两位小数）
print(round(cor_matrix, 2))
# 自定义双色调色板
flat_colors <- colorRampPalette(c("#e74c3c", "white", "#3498db"))(200)
# 绘图
corrplot(cor_matrix,
         method = "color",
         type = "full",
         col = flat_colors,
         addCoef.col = "black",
         tl.col = "black",
         tl.cex = 1.0,
         tl.srt = 45,
         tl.font = 3,              # 斜体标签
         number.cex = 0.9,
         mar = c(0, 0, 2, 0))
###############线性回归模型-20250504############################
# 安装和加载必要包
library(lm.beta)
library(ggplot2)
library(patchwork)

# 拟合线性回归模型
model <- lm(PlantDiversity ~ Elevation + Temperature + SoilMoisture +
              SoilBulk + SoilOrganicMatter, data = data)

# 原始回归系数提取
coef_df1 <- as.data.frame(summary(model)$coefficients)
coef_df1 <- coef_df1[-1, ]  # 去掉截距
coef_df1$Variable <- rownames(coef_df1)
rownames(coef_df1) <- NULL
names(coef_df1) <- c("Estimate", "Std.Error", "t.value", "P.Value", "Variable")

# 添加显著性标记
coef_df1$Sig <- cut(coef_df1$P.Value,
                    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                    labels = c("***", "**", "*", ""))

# 原始回归系数处理
coef_df1 <- as.data.frame(summary(model)$coefficients)[-1, ]
coef_df1$Variable <- rownames(summary(model)$coefficients)[-1]
names(coef_df1) <- c("Estimate", "Std.Error", "t.value", "P.Value", "Variable")
coef_df1$Sig <- cut(coef_df1$P.Value,
                    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                    labels = c("***", "**", "*", ""))
coef_df1$Signif <- coef_df1$P.Value < 0.05

# 图1：原始系数图（无图例）
ggplot(coef_df1, aes(y = reorder(Variable, Estimate), x = Estimate)) +
  geom_point(aes(color = Signif), size = 4) +
  geom_errorbarh(aes(xmin = Estimate - 1.96 * Std.Error, xmax = Estimate + 1.96 * Std.Error), height = 0.2) +
  geom_text(aes(label = Sig), hjust = -0.5, size = 5) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = c("TRUE" = "#2c3e50", "FALSE" = "gray70"), guide = "none") +
  labs(title = "Raw Coefficients", x = "Weight Estimate", y = "Variables") +
  theme_minimal(base_size = 13) +
  theme(panel.border = element_rect(fill = NA, color = "black"))

# 标准化回归系数处理
model_beta <- lm.beta(model)
coef_df2 <- as.data.frame(summary(model_beta)$coefficients)[-1, ]
coef_df2$Variable <- rownames(summary(model_beta)$coefficients)[-1]
coef_df2$Beta <- model_beta$standardized.coefficients[-1]
coef_df2$Lower <- coef_df2$Beta - 1.96 * coef_df2[, "Std. Error"]
coef_df2$Upper <- coef_df2$Beta + 1.96 * coef_df2[, "Std. Error"]
coef_df2$Sig <- cut(coef_df2$`Pr(>|t|)`,
                    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                    labels = c("***", "**", "*", ""))
coef_df2$Signif <- coef_df2$`Pr(>|t|)` < 0.05

# 图2：标准化系数图（无图例）
ggplot(coef_df2, aes(y = reorder(Variable, Beta), x = Beta)) +
  geom_point(aes(color = Signif), size = 4) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
  geom_text(aes(label = Sig), hjust = -0.5, size = 5) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = c("TRUE" = "#e67e22", "FALSE" = "gray70"), guide = "none") +
  labs(title = "Standardized Coefficients", x = "Standardized Beta", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.border = element_rect(fill = NA, color = "black"))



###########层次分隔
library(rdacca.hp)
# 设置响应变量和解释变量
Y <- data[, "PlantDiversity"]
X <- data[, c("Elevation", "Temperature", "SoilMoisture", "SoilBulk", "SoilOrganicMatter")]

# 正确调用（不要写 resp=, pred=, data=）
result <- rdacca.hp(Y, X, type = "adjR2", method = "RDA")

# 查看结果
print(result)

plot(result, plot.perc = TRUE)

# 提取变量名称和独立贡献（百分比）
plot_data <- data.frame(
  Variable = rownames(result$Hier.part),
  Contribution = result$Hier.part[, 1] * 100  # 乘100得到百分比
)

# 加载 ggplot2
library(ggplot2)

# 计算一个合理的顶部边界
ymax <- max(plot_data$Contribution) * 1.1
ggplot(plot_data, aes(y = reorder(Variable, Contribution), x = Contribution)) +
  geom_col(fill = "#2a9d8f", width = 0.6) +
  geom_text(aes(label = round(Contribution, 1)), hjust = -0.2, size = 3, family = "Arial") +
  labs(
    x = "Contribution (%)",
    y = "Variables"
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, max(plot_data$Contribution) * 1.1)
  ) +
  theme_minimal(base_family = "Arial", base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(color = "black")
  )



###############特征交互-20250505



##############广义可加模型GAM##################

# 加载必要的库
library(mgcv)

# 读取数据
data <- read.csv("Alpine_Meadow_Plant_Diversity_Data.csv")

# 拟合广义可加模型（GAM）
# 使用平滑函数 's()' 作用于每个自变量
gam_model <- gam(PlantDiversity ~ s(Elevation, k = 10) + s(Temperature, k = 10) +
                   s(SoilMoisture, k = 10) + s(SoilBulk, k = 10) + 
                   s(SoilOrganicMatter, k = 10), data = data)

# 查看模型摘要
summary(gam_model)

# 绘制模型，查看平滑函数的效果
plot(gam_model, pages = 1)




###########rdacca.hp包变差分解######################

# 加载rdacca.hp包
library(rdacca.hp)

# 定义响应变量和预测变量矩阵
response_var <- data$PlantDiversity
predictor_vars <- data[, c("Temperature", "Elevation", "SoilMoisture", "SoilBulk", "SoilOrganicMatter")]

# 加载数据并将每个变量包装为数据框
groups <- list(
  Temperature = data.frame(data$Temperature),
  Elevation = data.frame(data$Elevation),
  SoilMoisture = data.frame(data$SoilMoisture),
  SoilBulk = data.frame(data$SoilBulk),
  SoilOrganicMatter = data.frame(data$SoilOrganicMatter)
)

# 将响应变量包装为数据框
response_var <- data.frame(data$PlantDiversity)

# 执行方差分解分析
var_decomp_results <- rdacca.hp(response_var, groups)

# 查看结果
print(var_decomp_results)

# 可视化方差贡献
plot(var_decomp_results)



#############回归树-20250504####################

library(rpart)
library(rpart.plot)
data <- read.csv("New_Alpine_Meadow_Plant_Diversity_Data.csv")
# 限制最大深度为 5 层，允许较复杂结构但避免过拟合
tree_model <- rpart(
  PlantDiversity ~ Elevation + Temperature + SoilMoisture + SoilBulk + SoilOrganicMatter,
  data = data,
  method = "anova",
  control = rpart.control(maxdepth = 3, cp = 0.001, minsplit = 5)
)
View(data)
# 可视化
rpart.plot(tree_model, type = 2, extra = 101,
           fallen.leaves = TRUE,
           main = "Regression Tree for Plant Diversity (Max Depth = 3)")





#############随机森林模型-202050504###########
library(randomForest)
# 使用因变量 PlantDiversity 和其他自变量构建随机森林回归模型
set.seed(123) # 设置种子以确保结果可重复
#随机森林回归模型
rf_model <- randomForest(PlantDiversity ~ ., data = data, 
                         ntree = 500, importance = TRUE)

#############自变量重要性排序####################
# 加载必要的包
library(randomForest)

data <- read.csv("New_Alpine_Meadow_Plant_Diversity_Data.csv")
# 构建一个随机森林模型，预测变量为 mpg，使用其他变量作为自变量
set.seed(123)  # 设置随机种子以便结果可复现
rf_model <- randomForest(PlantDiversity ~ ., data = data, importance = TRUE)

# 查看特征的重要性
importance(rf_model)

# 将特征重要性进行排序并可视化(基于基尼重要性/Gini Importance)
importance_ordered <- importance(rf_model)[order(importance(rf_model)[,1],
                                                 decreasing = TRUE), ]
print(importance_ordered)

# 绘制特征重要性图
varImpPlot(rf_model, main = "Variable Importance in Random Forest")

# 查看特征重要性(基于置换重要性/Permutation Importance)
importance(rf_model, type = 2)
# 将特征重要性进行排序并可视化(置换重要性/Permutation Importance)
importance_ordered <- importance(rf_model)[order(importance(rf_model)[,1],
                                                 decreasing = TRUE), ]# 绘制特征重要性图
varImpPlot(rf_model, type = 2, main = "Variable Importance in Random Forest")
?varImpPlot










library(randomForest)

# 读取数据
data <- read.csv("New_Alpine_Meadow_Plant_Diversity_Data.csv")

# 构建随机森林模型
set.seed(123)
rf_model <- randomForest(PlantDiversity ~ ., data = data, importance = TRUE)

# ----------------------------
# 特征重要性排序并打印
# ----------------------------
# Permutation Importance（置换重要性）
perm_importance <- importance(rf_model, type = 2)
perm_ordered <- perm_importance[order(perm_importance[, 1], decreasing = TRUE), ]
print("Permutation Importance:")
print(perm_ordered)

# Gini Importance（基尼不纯度减少）
gini_importance <- importance(rf_model, type = 1)
gini_ordered <- gini_importance[order(gini_importance[, 1], decreasing = TRUE), ]
print("Gini Importance:")
print(gini_ordered)



# ----------------------------
# 组图显示（1行2列）
# ----------------------------

par(mfrow = c(1, 2))  # 设置画布为 1 行 2 列

varImpPlot(rf_model, type = 1, main = "Permutation Importance")          # 左图
varImpPlot(rf_model, type = 2, main = "Gini Importance")   # 右图

par(mfrow = c(1, 1))  # 恢复为单图模式







############部分依赖图PDF(基于随机森林-20250504)##############

# 加载包
library(randomForest)
library(pdp)

# 使用因变量 PlantDiversity 和其他自变量构建随机森林回归模型
set.seed(123) # 设置种子以确保结果可重复
#随机森林回归模型
rf_model <- randomForest(PlantDiversity ~ ., data = data, 
                         ntree = 500, importance = TRUE)

# 打印模型信息
print(rf_model)

# 查看变量重要性
importance(rf_model)
varImpPlot(rf_model)

# 使用pdp绘制部分依赖图
# 以Elevation为例，可以更改为其他自变量，如Temperature
partial_dependence <- partial(rf_model, pred.var = "Elevation", 
                              grid.resolution = 50)
# 绘制部分依赖图
plotPartial(partial_dependence)
# 也可以为其他变量绘制部分依赖图，例如Temperature
partial_dependence_temp <- partial(rf_model, pred.var = "Temperature", grid.resolution = 50)
plotPartial(partial_dependence_temp)



# 要绘图的变量
vars <- c("Elevation", "Temperature", "SoilMoisture", "SoilBulk", "SoilOrganicMatter")
# 用于存放每个图的列表
pdp_plot_list <- list()

# 使用 for 循环，逐个绘图
for (var in vars) {
  pd <- partial(rf_model, pred.var = var, grid.resolution = 50)
  
  p <- ggplot(pd, aes_string(x = var, y = "yhat")) +
    geom_line(size = 1, color = "#2c7bb6") +
    labs(
      title = paste("PDP:", var),
      x = var,
      y = "Predicted Plant Diversity"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  pdp_plot_list[[var]] <- p
}

# 补一个空白图让图形对齐 2x3
pdp_plot_list[["blank"]] <- plot_spacer()

# 拼图展示
final_plot <- wrap_plots(pdp_plot_list, ncol = 3) +
  plot_annotation(
    title = "Partial Dependence Plots for Key Environmental Variables",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

# 展示
print(final_plot)


############局部累计效应ALE(基于随机森林-20250504)##############

# 加载包
library(iml)
library(randomForest)

# 使用随机森林构建回归模型
set.seed(123)
rf_model <- randomForest(PlantDiversity ~ ., data = data, ntree = 500)


# 3. 创建 Predictor 对象
X <- data[, c("Elevation", "Temperature", "SoilMoisture", "SoilBulk", "SoilOrganicMatter")]
predictor_rf <- Predictor$new(rf_model, data = X, y = data$PlantDiversity)

# 4. 循环生成每个变量的 ALE 图
vars <- colnames(X)
ale_plot_list <- list()

for (v in vars) {
  ale_obj <- FeatureEffect$new(predictor_rf, feature = v, method = "ale")
  
# 基础图
p <- ale_obj$plot() +
    ggtitle(paste("ALE:", v)) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  #橘色、加粗线条
  if (!is.null(p$layers[[1]]$aes_params)) {
    p$layers[[1]]$aes_params$colour <- "#e67e22"  # 橘色
    p$layers[[1]]$aes_params$size <- 1          # 加粗
  }
  
  ale_plot_list[[v]] <- p
}

# 5. 添加空白图，保证拼图对齐（2x3）
ale_plot_list[["blank"]] <- patchwork::plot_spacer()

# 6. 拼图展示
final_ale_plot <- wrap_plots(ale_plot_list, ncol = 3) +
  plot_annotation(
    title = "Accumulated Local Effects (ALE) for Environmental Variables",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

# 7. 显示图像
print(final_ale_plot)





############特征交互-20250505############
set.seed(123) # 设置种子以确保结果可重复
#随机森林回归模型
rf_model <- randomForest(PlantDiversity ~ ., data = data, 
                         ntree = 500, importance = TRUE)
library(iml)
#####(1) 变量间的总体交互
# 构建 predictor（你应已有）
X <- data[, c("Elevation", "Temperature", "SoilMoisture", "SoilBulk", "SoilOrganicMatter")]
predictor_rf <- Predictor$new(rf_model, data = X, y = data$PlantDiversity)
# 计算每个变量与其他变量的平均交互强度
interaction_obj <- Interaction$new(predictor_rf)
# 查看结果
print(interaction_obj$results)
# 绘图
plot1 <- plot(interaction_obj) +
  ggtitle("Overall Interaction Strength") +
  theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.minor = element_blank()
  ) +
  scale_color_manual(values = c("black")) +
  scale_x_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  )

######(2)两两变量间的交互响应：以Elevation 与其他变量的交互强度为例
interaction_elev <- Interaction$new(predictor_rf, feature = "Elevation")
print(interaction_elev$results)
# 绘图（可视化 Elevation 与其他变量的 H-statistic）
plot2 <- plot(interaction_elev) +
  ggtitle("Interaction with Elevation") +
  theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.minor = element_blank()
  ) +
  scale_color_manual(values = c("black")) +
  scale_x_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  )
#组图
final_plot <- plot1 + plot2 +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(
    title = "Feature Interaction Summary",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = 16, hjust = 0.5))
  )
#显示
print(final_plot)


######2D-PDP和2D-ALE
#加载包
library(randomForest)
library(pdp)
library(iml)
library(ggplot2)
library(patchwork)
library(yaImpute)
#计算Elevation与Temperature的PDP 
pdp_2d <- partial(rf_model,
                  pred.var = c("Elevation", "Temperature"),
                  grid.resolution = 20,
                  progress = "none")
#绘制PDP
p1 <- ggplot(pdp_2d, aes(x = Elevation, y = Temperature, fill = yhat)) +
  geom_tile() +
  #geom_contour(aes(z = yhat), color = "white", alpha = 0.5) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = mean(pdp_2d$yhat, na.rm = TRUE),
    name = "Predicted\nDiversity"
  ) +
  labs(
    title = "2D PDP: Elevation × Temperature",
    x = "Elevation", y = "Temperature"
  ) +
  theme_minimal(base_size = 12)

#计算ALE
ale_effect <- FeatureEffect$new(
  predictor_rf,
  feature = c("Elevation", "Temperature"),
  method = "ale",
  grid.size = 20  # 或更高
)
p2 <-plot(ale_effect)
# 用 %+% 和 + 来修改颜色
p2 <- p2 + scale_fill_gradient2(
  low = "blue", mid = "white", high = "red", midpoint = 0,
  name = "ALE"
)

#组图
final_plot <- p1  + p2 + plot_layout(guides = "collect")
print(final_plot)

##############全局代理模型-20250505##############
#1.加载库
library(randomForest)
library(ggplot2)
library(broom)
library(dplyr)
#2. 构建随机森林模型（黑箱模型）
library(randomForest)
# 使用因变量 PlantDiversity 和其他自变量构建随机森林回归模型
set.seed(123) # 设置种子以确保结果可重复
#随机森林回归模型
rf_model <- randomForest(PlantDiversity ~ ., data = data, 
                         ntree = 500, importance = TRUE)
#3. 获取随机森林预测结果（作为 surrogate 模型的目标变量）
data$y_rf <- predict(rf_model, newdata = data)
# 4. 构建线性代理模型（用 RF 的预测值做因变量）
lm_surrogate <- lm(y_rf ~ Elevation + Temperature + SoilMoisture + SoilBulk + SoilOrganicMatter, data = data)
# 5. 模型摘要
summary(lm_surrogate)  # 查看系数显著性和 R²
# 6. 可视化：变量估计值与置信区间
coef_df <- tidy(lm_surrogate, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    signif_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE            ~ ""
    )
  )
ggplot(coef_df, aes(x = reorder(term, estimate), y = estimate)) +
  geom_point(size = 4, color = "#e74c3c", stroke = 1.1, shape = 21, fill = "#e74c3c") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_text(aes(label = signif_label), hjust = -0.6, size = 5, color = "black") +
  coord_flip() +
  labs(
    title = "Global Surrogate (Linear Model)",
    x = "Variables", y = "Estimated Effect"
  ) +
  theme_bw(base_size = 13) +  # ✅ 替换 minimal
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin = margin(5.5, 30, 5.5, 5.5)  # 右边加大防止星号被裁剪
  )



##############局部代理模型LIME-20250505#########################
# 安装并加载必要的包
library(lime)
library(randomForest)

# 假设你已经有一个随机森林模型 'rf_model'，以及数据 'data_frame'
# 使用一个随机森林模型并进行训练
set.seed(123)
rf_model <- randomForest(PlantDiversity ~ ., data = data, ntree = 500)

# 定义model_type方法 - 指定模型类型为回归
model_type.randomForest <- function(x, ...) {
  return("regression")  # 如果是分类问题，返回 "classification"
}

# 定义predict_model方法 - 用于预测模型的输出
predict_model.randomForest <- function(x, newdata, ...) {
  results <- predict(x, newdata)
  return(as.data.frame(results))
}

# 创建解释器对象
explainer <- lime(data, rf_model)

# 选择一个特定的实例进行解释
instance_to_explain <- data[30,]  # 选择数据框的第一行

# 使用LIME解释该数据点
explanation <- explain(instance_to_explain, explainer, n_labels = 1, n_features = 6)
library(dplyr)
explanation_filtered <- explanation %>%
  filter(!grepl("y_rf", feature))
# 打印解释结果
print(explanation)
# 通过可视化查看特征贡献
plot_features(explanation_filtered)

###########SHAP局部解释-202050505###################
# 加载包
library(iml)
library(randomForest)

# 使用因变量 PlantDiversity 和其他自变量构建随机森林回归模型
set.seed(123)
rf_model <- randomForest(PlantDiversity ~ ., data = data, ntree = 500)
# 创建Predictor对象
predictor <- Predictor$new(rf_model, data = data)
# 计算SHAP值
shapley <- Shapley$new(predictor, x.interest = data[30, , drop = FALSE])# 解释第30个实例
# 打印结果
print(shapley)
# 绘制SHAP值
shapley$plot()


###ggplot2绘图
# 明确移除响应变量，只用特征数据
X <- data[, setdiff(names(data), "PlantDiversity")]
# 正确创建 predictor
predictor <- Predictor$new(rf_model, data = X, y = data$PlantDiversity)
# 选择一个观测点
x_interest <- X[30, , drop = FALSE]
# 计算 SHAP 值
shapley <- Shapley$new(predictor, x.interest = x_interest)
library(ggplot2)
library(dplyr)
# 提取 SHAP 数据并处理 feature + value 标签
shap_df <- shapley$results %>%
  filter(feature != "y_rf") %>%  # 移除 y_rf 行
  mutate(
    value_num = as.numeric(sub(".*=", "", feature.value)),  # 提取等号后面的值并转为数值
    feature_value = paste0(feature, " = ", round(value_num, 2))  # 拼接成特征标签
  ) %>%
  arrange(phi) %>%
  mutate(feature_value = factor(feature_value, levels = feature_value))  # 控制 y 轴顺序

# 获取预测值和 baseline（平均预测）
prediction <- as.numeric(shapley$y.hat)
baseline <- as.numeric(shapley$y.hat.mean)

# 绘图
ggplot(shap_df, aes(x = phi, y = feature_value, fill = phi > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#e74c3c", "FALSE" = "#3498db")) +
  labs(
    title = "SHAP Explanation",
    subtitle = paste0("Prediction: 38.27", round(prediction, 2),
                      "  |  Baseline: 35.01", round(baseline, 2)),
    x = "SHAP value (phi)",
    y = "Feature = Value"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")





