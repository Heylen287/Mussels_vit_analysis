library(readxl)
library(dplyr)
library(Momocs)
library(ggplot2)
library(vegan)
library(ggforce)

file_path <- "E:/мидийное/mussels morphology analysis/Морфометрия Витязь.xlsx"
df_morph <- read_excel(file_path, sheet = 4)

df_morph  <- df_morph [!is.na(df_morph $ID), ]
df_morph 


#outlines

folder <- paste("data/", list.files("data"), sep = "")

coords <- import_txt(txt.paths = folder)

myt <- Out(coords)
myt 
names(myt)

df_names_myt <- data.frame(ID = (names(myt)))

coo_check(myt)

str(myt)

panel(myt, names = FALSE)

df_gen <- merge(df_names_myt, df_morph, sort = F)

df_gen %<>% 
  filter(as.character(ID) %in% names(myt$coo))

myt$fac <- df_gen
nrow(df_gen)
pca_results <- myt %>%                  
  efourier(6, norm = TRUE) %>%  
  PCA()

# Смотрим на график (пока без группировки)
plot_PCA(pca_results, points = TRUE, morphospace = TRUE)

myt %T>%                  
  stack() %>%                  
  efourier(6, norm=TRUE) %>%  
  PCA()  %T>%                 
  plot_PCA(points = T, morphospace = T, chull= T, chullfilled =T, labelpoints = T, labelgroups = T) ->
  pca_results

str(pca_results)

summary(pca_results)

df_out <- as.data.frame(pca_results$x[,1:2]) #главные компоненты PCA по аутлайнам

df_out$ID <- row.names(df_out)

df_out$ID <- gsub(x = df_out$ID, pattern = "v", replacement = "V")



df_morph %>% 
  mutate(category = case_when(`Structure_MG` < 0.05 ~ "MT_pure", 
                              `Structure_MG` < 0.25 & `Structure_MG` > 0.05 ~ "MT_like_hybrids",
                              `Structure_MG` <= 0.75 & `Structure_MG` >= 0.25 ~ "F1",
                              `Structure_MG` > 0.75 & `Structure_MG` < 0.95 ~ "MG_like_hybrids",
                              `Structure_MG` > 0.95 ~ "MG_pure")) ->
  df_morph
df_morph$category <- factor(df_morph$category, levels = c("MT_pure", "MT_like_hybrids", "F1", "MG_like_hybrids", "MG_pure"))

df_all <- merge(df_morph, df_out, by = "ID",) #таблица с с данными STRUCTURE, классифицированными генотипами, замерами по (McDonald et al., 1991) и главными компонентами PCA по аутлайнам

theme_set(theme_classic())

df_all %>% 
  ggplot(aes(x = Structure_MG, y = PC1)) +
  geom_point() +
  geom_smooth(method = "gam")

df_all %>% 
  ggplot(aes(x = category, y = PC1)) +
  geom_boxplot() +
  geom_hline(yintercept = 0)

df_all %>% 
  ggplot(aes(x = Structure_MG, y = PC2)) +
  geom_point() +
  geom_smooth(method = "gam")

df_all %>% 
  ggplot(aes(x = category, y = PC2)) +
  geom_boxplot() +
  geom_hline(yintercept = 0)

df_all %>% 
  ggplot(aes(PC1, PC2, colour = Morphotype)) +
  geom_point()
  
df_all %>% 
  ggplot(aes(x = Morphotype, y = PC2)) +
  geom_boxplot() +
  geom_hline(yintercept = 0)

library(ggplot2)
library(interp)        # интерполяция
library(metR)          # подписи изолиний
library(ggnewscale)    # отдельные цветовые шкалы для точек и линий

grid <- with(df_all, interp(x = PC1, y = PC2, z = Structure_MG,
                        nx = 120, ny = 120,
                        duplicate = "mean",
                        extrap = FALSE))
grid_df <- as.data.frame(grid)   # столбцы x, y, z

ggplot(df_all, aes(PC1, PC2)) +
  geom_point(size = 2) +        # облако точек
  new_scale_color() +                                # сброс шкалы
  geom_contour(data = grid_df,
               aes(x, y, z = z, color = after_stat(level)),
               bins = 8, linewidth = 0.6) +
  geom_text_contour(data = grid_df, aes(x, y, z = z),
                    stroke = 0.2, size = 3) +
  scale_color_viridis_c(name = "z") +
  coord_fixed() +
  theme_minimal()

head(df_all)

df_all %>% 
  mutate(MG = df_morph$MG, ME = df_morph$ME, MT = df_morph$MT, Site = df_morph$site, Label = row.names(df_morph), Genotype = df_morph$Genotype) %>% 
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(color = factor(Genotype), size = MG, shape = Site)) +
  scale_color_manual(values = c("darkgreen", "blue", "orange", "red")) +
  scale_shape_manual(values = c(15, 16, 17, 18)) -> Pl3

#PCA по всему
morph_all <- morph_log

morph_all$PC1 <- df_all$PC1
morph_all$PC2 <- df_all$PC2
morph_all

morph_all.rda <- rda(morph_all, scale = F)

plot(morph_all.rda, type = "t")

as.data.frame(scores(morph_all.rda, choices = 1:3)$sites) %>% 
  mutate(MG = df_all$MG, ME = df_all$ME, MT = df_all$MT, Site = df_all$site, Label = row.names(df_all), Genotype = df_all$Genotype) %>% 
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(color = factor(Genotype), size = MG, shape = Site)) +
  scale_color_manual(values = c("darkgreen", "blue", "orange", "red")) +
  scale_shape_manual(values = c(15, 16, 17, 18)) -> Pl4

library(gridExtra)
grid.arrange(Pl1, Pl2, Pl3, Pl4)

#Кажется, что замеры по (McDonald et al., 1991) работают лучше, чем аутлайны и замеры+аутлайны



#' Если я правильно понимаю, то наша задача понять какие из признаков более всего коррелируют с генетическим статусом мидии. Последний выражается в индексах STRUCTURE (MG, ME, MT).
#' Стало быть это задача для какого-то регрессионного анализа. Типа как признак такой-то зависит от MG, например.  
#' В данном случае самая подходящая форма регрессионного анализа - это ординация в ограниченных осях, то есть RDA, CCA или dbRDA. Чтобы не мучить людей лишними знаниями предлагаю остановиться на RDA.
#' Задача этого анализа выявить те признаки, которые более всего связаны с Генотипами. 

morph_log$PC1_outline <- df_all$PC1
morph_log$PC2_outline <- df_all$PC2


all.rda <- rda(morph_log ~ MG + MT +  ht + wid, data = df_morph, scale = T) 

#' Беру MT чисто из соображений, что они, как кажется, более дальние родственники для MG, чем ME (или я не прав?)
#' Поскольку влияние размеров в датасете morph_log не удалено, то в модель включаем что-то что характеризует размеры. Пусть это будет wid и ht.

vif.cca(all.rda)

anova(all.rda)

anova(all.rda, by = "axis")
anova(all.rda, by = "margin")


plot(all.rda, type = "t", display = c("species", "cn"))

#'Из этой картинки видно, что все метрические признаки, несмотря на то, что они были нормированы по длине (ведь так?) имеют связь с размером. Это не супер интересный результат. Важно то, что RDA1 -это однозначно размер. 
#'Поэтому изображаем в осях RDA2-RDA3


plot(all.rda, choice = 2:3, type = "t", display = c("species", "cn"))

#'Вторая ось - это как раз что-то про MT vs MG.
#' Ну а дальше pal, hp, lpr, aam - это то, чего много у MT, а wpr, lig (ну и то, что рядом) это то, чего много у MG.
#' Да, важно, что PC1_outline  и PC2_outline похоже не коррелируют с генотипами. То есть форма раковины лишь ситуативно выделяет gallo, а в целом не работает.

# Вот например для троссулусного маркера
qplot(x = df_morph$MT, y = morph_log$hp) + 
  geom_smooth(method = "lm")

qplot(x = df_morph$MG, y = morph_log$hp) + 
  geom_smooth(method = "lm")


####
qplot(x = df_morph$MG, y = morph_log$PC1_outline) + 
  geom_smooth(method = "lm")
qplot(x = df_morph$MT, y = morph_log$PC1_outline) + 
  geom_smooth(method = "lm")

###



# А это для галло
qplot(x = df_morph$MT, y = morph_log$wpr) + 
  geom_smooth(method = "lm")

qplot(x = df_morph$MG, y = morph_log$wpr) + 
  geom_smooth(method = "lm")


# Далее я предлагаю оставит в анализе только те маркеры, которые сильно коррелируют с RDA2 и уже по ним построить PCA для того, чтобы понять как ординируются особи из разных мест. 


