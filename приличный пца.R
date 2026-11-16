library(readxl)
library(dplyr)
library(magrittr)
library(Momocs)
library(ggplot2)

# ============ 1. Данные и контуры ============
file_path <- "E:/мидийное/mussels morphology analysis/Морфометрия Витязь.xlsx"
df_morph  <- read_excel(file_path, sheet = 4)
df_morph  <- df_morph[!is.na(df_morph$ID), ]

folder <- paste("data/", list.files("data"), sep = "")
coords <- import_txt(txt.paths = folder)
myt    <- Out(coords)

df_names_myt <- data.frame(ID = names(myt))
df_gen <- merge(df_names_myt, df_morph, sort = F)
df_gen %<>% filter(as.character(ID) %in% names(myt$coo))
myt$fac <- df_gen

pca_results <- myt %>%
  efourier(6, norm = TRUE) %>%
  PCA()

df_out <- as.data.frame(pca_results$x[, 1:2])
df_out$ID <- row.names(df_out)
df_out$ID <- gsub(x = df_out$ID, pattern = "v", replacement = "V")

df_morph %>%
  mutate(category = case_when(Structure_MG <  0.05 ~ "MT_pure",
                              Structure_MG <  0.25 ~ "MT_like_hybrids",
                              Structure_MG <= 0.75 ~ "F1",
                              Structure_MG <= 0.95 ~ "MG_like_hybrids",
                              Structure_MG >  0.95 ~ "MG_pure")) -> df_morph
df_morph$category <- factor(df_morph$category,
                            levels = c("MT_pure", "MT_like_hybrids", "F1",
                                       "MG_like_hybrids", "MG_pure"))

df_all <- merge(df_morph, df_out, by = "ID")

# ============ 2. Параметры графика ============
n_cols       <- 6      # контуров в сетке по горизонтали
n_rows       <- 5      # контуров в сетке по вертикали
shape_size   <- 0.9    # размер контура относительно ячейки сетки (0.9 = почти вплотную)
plot_aspect  <- 0.75   # высота панели / ширина панели (больше = выше график)
pad          <- 0.08   # запас вокруг облака точек
out_dir      <- "."

pal_cat <- c(MT_pure         = "#A00E00FF",
             MT_like_hybrids = "#D04E00FF",
             F1              = "#F6C200FF",
             MG_like_hybrids = "#0086A8FF",
             MG_pure         = "#132B69FF")

lab_cat <- c(MT_pure         = "MT_pure",
             MT_like_hybrids = "MT_like_hybrids",
             F1              = "F1",
             MG_like_hybrids = "MG_like_hybrids",
             MG_pure         = "MG_pure")

stopifnot(all(c("PC1", "PC2", "category") %in% names(df_all)))
stopifnot(!anyNA(df_all$category))
stopifnot(identical(pca_results$method[1], "efourier"))

# ============ 3. Подписи осей ============
var_expl <- 100 * pca_results$sdev^2 / sum(pca_results$sdev^2)
lab_x <- sprintf("PC1 (%.1f%%)", var_expl[1])
lab_y <- sprintf("PC2 (%.1f%%)", var_expl[2])

# ============ 4. Контуры морфопространства ============
efourier_xy <- function(coe, nb_pts = 80) {
  nb_h <- length(coe) / 4
  an <- coe[1:nb_h];                    bn <- coe[(nb_h + 1):(2 * nb_h)]
  cn <- coe[(2 * nb_h + 1):(3 * nb_h)]; dn <- coe[(3 * nb_h + 1):(4 * nb_h)]
  t  <- seq(0, 2 * pi, length.out = nb_pts + 1)
  x <- y <- numeric(length(t))
  for (h in seq_len(nb_h)) {
    x <- x + an[h] * cos(h * t) + bn[h] * sin(h * t)
    y <- y + cn[h] * cos(h * t) + dn[h] * sin(h * t)
  }
  cbind(x = x, y = y)
}

rx <- range(df_all$PC1); ry <- range(df_all$PC2)
x_lim <- rx + c(-1, 1) * pad * diff(rx)
y_lim <- ry + c(-1, 1) * pad * diff(ry)
Sx <- diff(x_lim); Sy <- diff(y_lim)

# центры ячеек сетки: ячейки равномерно покрывают всю панель
grid_pos <- expand.grid(gx = x_lim[1] + (seq_len(n_cols) - 0.5) * Sx / n_cols,
                        gy = y_lim[1] + (seq_len(n_rows) - 0.5) * Sy / n_rows)

rot    <- pca_results$rotation[, 1:2]
mshape <- if (!is.null(pca_results$mshape)) pca_results$mshape else pca_results$center

raw_shapes <- lapply(seq_len(nrow(grid_pos)), function(i) {
  coe <- mshape + as.numeric(rot %*% c(grid_pos$gx[i], grid_pos$gy[i]))
  xy  <- efourier_xy(coe)
  
  # Центрируем контур
  xy[, 1] <- xy[, 1] - mean(range(xy[, 1]))
  xy[, 2] <- xy[, 2] - mean(range(xy[, 2]))
  
  # ОТРАЖЕНИЕ ПО ВЕРТИКАЛИ (по оси Y)
  xy[, 2] <- -xy[, 2] 
  
  xy
})
# единый масштаб для всех контуров: самый большой вписывается в ячейку
# (размеры в долях панели: ширина = 1, высота = plot_aspect)
w_max <- max(sapply(raw_shapes, function(s) diff(range(s[, 1]))))
h_max <- max(sapply(raw_shapes, function(s) diff(range(s[, 2]))))
scale_vis <- shape_size * min((1 / n_cols) / w_max, (plot_aspect / n_rows) / h_max)

# перевод из долей панели в единицы PC1 и PC2 (контуры не искажаются)
ux <- Sx
uy <- Sy / plot_aspect

shapes_df <- bind_rows(lapply(seq_along(raw_shapes), function(i) {
  data.frame(id = i,
             x = raw_shapes[[i]][, 1] * scale_vis * ux + grid_pos$gx[i],
             y = raw_shapes[[i]][, 2] * scale_vis * uy + grid_pos$gy[i])
}))

# ============ 5. Графики ============
base_plot <- function() {
  ggplot() +
    geom_polygon(data = shapes_df, aes(x, y, group = id),
                 fill = NA, colour = "grey70", linewidth = 0.6) +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.4) +
    coord_cartesian(xlim = x_lim, ylim = y_lim, expand = FALSE) +
    labs(x = lab_x, y = lab_y) +
    theme_minimal(base_size = 16) +
    theme(aspect.ratio    = plot_aspect,
          panel.grid      = element_blank(),
          axis.line       = element_line(colour = "grey30", linewidth = 0.4),
          axis.ticks      = element_line(colour = "grey30", linewidth = 0.4),
          axis.title      = element_text(colour = "grey20"),
          legend.position = "bottom",
          legend.title    = element_blank())
}

# Вариант 1: одноцветные точки
p_plain <- base_plot() +
  geom_point(data = df_all, aes(PC1, PC2),
             colour = "#68228B", size = 2.6, alpha = 0.7)

# Вариант 2: точки по category
p_cat <- base_plot() +
  geom_point(data = df_all, aes(PC1, PC2, colour = category),
             size = 2.6, alpha = 0.9) +
  scale_colour_manual(values = pal_cat, labels = lab_cat, drop = FALSE) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                               override.aes = list(size = 3.5, alpha = 0.8)))

print(p_plain)
print(p_cat)

ggsave(file.path(out_dir, "PCA_outlines_plain.png"),    p_plain, width = 8, height = 7, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "PCA_outlines_category.png"), p_cat,   width = 8, height = 7, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "PCA_outlines_plain.pdf"),    p_plain, width = 8, height = 7)
ggsave(file.path(out_dir, "PCA_outlines_category.pdf"), p_cat,   width = 8, height = 7)


# ============ Боксплоты PC1/PC2 по категориям (без легенды) ============
library(ggplot2)

pal_cat <- c(MT_pure         = "#A00E00FF",
             MT_like_hybrids = "#D04E00FF",
             F1              = "#F6C200FF",
             MG_like_hybrids = "#0086A8FF",
             MG_pure         = "#132B69FF")

lab_short <- c(MT_pure         = "M. trossulus",
               MT_like_hybrids = "MT-like",
               F1              = "F1",
               MG_like_hybrids = "MG-like",
               MG_pure         = "M. galloprovincialis")

var_expl <- 100 * pca_results$sdev^2 / sum(pca_results$sdev^2)
lab_pc1  <- sprintf("PC1 (%.1f%%)", var_expl[1])
lab_pc2  <- sprintf("PC2 (%.1f%%)", var_expl[2])

theme_box <- theme_minimal(base_size = 16) +
  theme(panel.grid        = element_blank(),
        axis.line         = element_line(colour = "grey30", linewidth = 0.4),
        axis.ticks        = element_line(colour = "grey30", linewidth = 0.4),
        axis.title        = element_text(colour = "grey20"),
        axis.text         = element_text(colour = "grey20"),
        legend.position   = "none",
        plot.margin       = margin(10, 15, 10, 10))

make_box <- function(data, yvar, ylab) {
  ggplot(data, aes(x = category, y = .data[[yvar]], fill = category)) +
    geom_hline(yintercept = 0, colour = "grey45",
               linewidth = 0.4, linetype = "dashed") +
    geom_boxplot(width = 0.62,
                 outlier.shape = 21, outlier.size = 2,
                 outlier.fill = "white", outlier.colour = "grey20",
                 colour = "grey20", linewidth = 0.5, alpha = 0.9) +
    scale_fill_manual(values = pal_cat, drop = FALSE) +
    scale_x_discrete(labels = lab_short) +
    labs(x = NULL, y = ylab) +
    theme_box
}

p_box_pc1 <- make_box(df_all, "PC1", lab_pc1)
p_box_pc2 <- make_box(df_all, "PC2", lab_pc2)

print(p_box_pc1)
print(p_box_pc2)

ggsave("boxplot_PC1_category.png", p_box_pc1, width = 8, height = 5, dpi = 300, bg = "white")
ggsave("boxplot_PC2_category.png", p_box_pc2, width = 8, height = 5, dpi = 300, bg = "white")
ggsave("boxplot_PC1_category.pdf", p_box_pc1, width = 8, height = 5)
ggsave("boxplot_PC2_category.pdf", p_box_pc2, width = 8, height = 5)
