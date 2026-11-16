#это тот код с обработкой данных по классической морфометрии который я откопала на диске, Ваш файл morph я пока нормально не осознала так что здесь он не учтен
library(readxl)
library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(mgcv)

set.seed(123)
theme_set(theme_bw(base_size = 14))

file_path <- "data/Морфометрия Витязь.xlsx"
dat <- read_excel(file_path, sheet = 4)
dat <- dat[!is.na(dat$ID), ]

traits <- c("aam", "dpr", "hp", "ht", "lar", "lig", "lpr", "pad", "padp",
            "padv", "pal", "ppad", "ular", "vpr", "wid", "wpr", "tkn")

stopifnot(!anyNA(dat[, traits]))
table(dat$Morphotype)   

dat %>%
  mutate(category = case_when(Structure_MG < 0.05 ~ "MT_pure",
                              Structure_MG < 0.25 ~ "MT_like_hybrids",
                              Structure_MG <= 0.75 ~ "F1",
                              Structure_MG < 0.95 ~ "MG_like_hybrids",
                              TRUE ~ "MG_pure")) -> dat
dat$category <- factor(dat$category, levels = c("MT_pure", "MT_like_hybrids", "F1",
                                                "MG_like_hybrids", "MG_pure"))
table(dat$category)

dat %>%
  filter(!is.na(Z)) %>%
  ggplot(aes(x = category, y = Z)) +
  geom_boxplot()

dat %>% mutate(Z_zero = ifelse(Z == 0, 1, 0)) -> dat

mod_T <- gam(Z_zero ~ s(Structure_MG), family = "binomial", data = dat)
summary(mod_T)

Mydata <- data.frame(Structure_MG = seq(0, 1, 0.01))
pr <- predict(mod_T, newdata = Mydata, type = "link", se.fit = TRUE)
Mydata$fit <- plogis(pr$fit)
Mydata$lo  <- plogis(pr$fit - 1.96 * pr$se.fit)
Mydata$hi  <- plogis(pr$fit + 1.96 * pr$se.fit)

p_gam <- ggplot(Mydata, aes(x = Structure_MG, y = fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70", alpha = 0.5) +
  geom_line(linewidth = 1) +
  geom_jitter(data = filter(dat, !is.na(Z_zero)), aes(y = Z_zero),
              height = 0.03, width = 0, alpha = 0.5) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(x = "Вклад аллелей M. galloprovincialis",
       y = "Вероятность Z = 0")
p_gam

# блок трансформации данных + пца
morph <- as.data.frame(dat[, traits])
rownames(morph) <- dat$ID

env <- data.frame(MG         = dat$Structure_MG,
                  logL       = log(dat$L),
                  L          = dat$L,
                  Morphotype = factor(dat$Morphotype),
                  row.names  = dat$ID)
summary(env)

morph_log  <- log(morph)
morph_log2 <- morph_log / env$logL

fill_mg <- scale_fill_gradient(low = "#00FF7F", high = "#530096",
                               name = "Вклад аллелей\nM. galloprovincialis")

# PCA остатка
pca_res <- rda(morph_log2 ~ logL, data = env)
pca_res

pc_scores <- as.data.frame(scores(pca_res, choices = 1:4)[2])
pc_traits <- data.frame(pc_scores, env)

lab_res <- function(i) paste0("PC", i, " остатка (",
                              round(100 * pca_res$CA$eig[i] / pca_res$tot.chi, 1),
                              " % общей вариации)")

plot_pca_res <- function(i, j) {
  st <- as.data.frame(scores(pca_res, display = "sites",   choices = c(i, j) + 1, scaling = 2))
  sp <- as.data.frame(scores(pca_res, display = "species", choices = c(i, j) + 1, scaling = 2))
  colnames(st) <- colnames(sp) <- c("x", "y")
  st$MG <- env$MG
  
  # стрелки растянуты под облако точек
  span <- max(abs(st[, 1:2]))
  sp[, 1:2] <- sp[, 1:2] * 0.9 * span / max(abs(sp[, 1:2]))
  sp$label <- rownames(sp)
  
  # подписи стрелок
  len <- sqrt(sp$x^2 + sp$y^2)
  sh  <- ifelse(len > 0.35 * span, 0.055 * span, 0.02 * span)
  sp$lab_x <- sp$x + sh * sp$x / len
  sp$lab_y <- sp$y + sh * sp$y / len
  
  ggplot() +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(data = st, aes(x = x, y = y, fill = MG),
               shape = 21, size = 3, colour = "black") +
    geom_segment(data = sp, aes(x = 0, y = 0, xend = x, yend = y),
                 arrow = arrow(length = unit(0.18, "cm")), colour = "grey30") +
    geom_text_repel(data = sp, aes(x = lab_x, y = lab_y, label = label),
                    colour = "grey20", size = 4, max.overlaps = Inf,
                    box.padding = 0.1, point.padding = 0, force = 0.3,
                    segment.colour = NA, seed = 1) +
    fill_mg +
    labs(x = lab_res(i), y = lab_res(j))
}

p_res1 <- plot_pca_res(1, 2)
p_res1

p_res1b <- plot_pca_res(2, 3)
p_res1b

p_res1c <- ggplot(pc_traits, aes(x = sites.PC1, y = MG, colour = Morphotype)) +
  geom_point(size = 3) +
  labs(x = lab_res(1), y = "Вклад аллелей M. galloprovincialis\nв генотипы особей")
p_res1c


# еще про морфотипы
p_res2 <- ggplot(pc_traits, aes(x = sites.PC1)) +
  geom_density(aes(fill = Morphotype), alpha = 0.5) +
  labs(x = lab_res(1))
p_res2

p_res3 <- ggplot(pc_traits, aes(x = sites.PC2)) +
  geom_density(aes(fill = Morphotype), alpha = 0.5) +
  labs(x = lab_res(2))
p_res3

# MG против осей остатка, панели по морфотипам
connect_MG_PC1 <- ggplot(pc_traits, aes(x = MG, y = sites.PC1)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~Morphotype) +
  labs(y = lab_res(1))
connect_MG_PC1

connect_MG_PC2 <- ggplot(pc_traits, aes(x = MG, y = sites.PC2)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~Morphotype) +
  labs(y = lab_res(2))
connect_MG_PC2

# двумерная плотность по морфотипам
plot_density_mg <- function(d, y, size_var = NULL) {
  yr <- range(d[[y]])
  p <- ggplot(d, aes(x = MG, y = .data[[y]])) +
    stat_density_2d(aes(fill = Morphotype, alpha = after_stat(level)),
                    geom = "polygon", show.legend = c(size = FALSE)) +
    geom_density_2d(aes(colour = Morphotype), linewidth = 0.4,
                    show.legend = c(size = FALSE)) +
    scale_alpha_continuous(range = c(0.15, 0.6), guide = "none") +
    scale_x_continuous(limits = c(-0.5, 1.5)) +
    scale_y_continuous(limits = yr + c(-1, 1) * diff(yr)) +
    coord_cartesian(xlim = c(-0.05, 1.05), ylim = yr + c(-0.05, 0.05) * diff(yr))
  if (is.null(size_var)) p + geom_point(alpha = 0.6)
  else p + geom_point(aes(size = .data[[size_var]]), alpha = 0.6)
}

p_dens1 <- plot_density_mg(pc_traits, "sites.PC1") + labs(y = lab_res(1))
p_dens1
