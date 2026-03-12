# Clear the workspace
rm(list = ls())

# Load required packages
library(ggplot2)
library(patchwork)

# ============================================================
# Helper functions for equilateral mixture triangle geometry
# ============================================================

sqrt3 <- sqrt(3)

# Convert barycentric coordinates (carb, protein, fat) to 2D Cartesian
# coordinates for an equilateral mixture triangle
bary_to_xy <- function(carb, protein, fat) {
  x <- protein + 0.5 * fat
  y <- (sqrt3 / 2) * fat
  data.frame(x, y)
}

# Define triangle vertices
v1 <- bary_to_xy(1, 0, 0)
v2 <- bary_to_xy(0, 1, 0)
v3 <- bary_to_xy(0, 0, 1)
tri2d <- rbind(v1, v2, v3, v1)

# Build internal grid lines at specified mixture levels
mkseg <- function(p1, p2) {
  data.frame(x = p1$x, y = p1$y, xend = p2$x, yend = p2$y)
}

levs <- c(0.2, 0.4, 0.6, 0.8)

grid_carb <- do.call(rbind, lapply(levs, function(t) {
  mkseg(bary_to_xy(t, 0, 1 - t), bary_to_xy(t, 1 - t, 0))
}))

grid_prot <- do.call(rbind, lapply(levs, function(t) {
  mkseg(bary_to_xy(0, t, 1 - t), bary_to_xy(1 - t, t, 0))
}))

grid_fat <- do.call(rbind, lapply(levs, function(t) {
  mkseg(bary_to_xy(0, 1 - t, t), bary_to_xy(1 - t, 0, t))
}))

# ============================================================
# Define food compositions and target points
# ============================================================

averagingRatio <- 0.075

# Food items in barycentric coordinates
carb    <- c(0.10, 0.30, 0.40, 0.30, 0.15)
protein <- c(0.30, 0.30, 0.50, 0.50, 0.80)
fat     <- c(0.60, 0.40, 0.10, 0.20, 0.05)

foods <- cbind(carb = carb, protein = protein, fat = fat)
foods_df <- as.data.frame(foods)

# Project foods into EMT coordinates
foods_xy <- bary_to_xy(foods_df$carb, foods_df$protein, foods_df$fat)
foods_df <- cbind(
  food = factor(c("i", "ii", "iii", "iv", "v"), levels = c("i", "ii", "iii", "iv", "v")),
  foods_df,
  foods_xy
)

# Create a separate label data frame so individual labels can be nudged
foods_labels_df <- foods_df
foods_labels_df$x[foods_labels_df$food == "iv"] <- foods_labels_df$x[foods_labels_df$food == "iv"] + 0.055
foods_labels_df$y[foods_labels_df$food == "iv"] <- foods_labels_df$y[foods_labels_df$food == "iv"] - 0.045

# Define intake target and its EMT coordinates
I <- c(carb = 0.2, protein = 0.4, fat = 0.4)
I_xy <- bary_to_xy(I["carb"], I["protein"], I["fat"])

# Compute the center of mass of the food items
zbar <- colMeans(foods)
zbar_xy <- bary_to_xy(zbar[1], zbar[2], zbar[3])

# Define the starting intake state
O0 <- c(carb = 0.6, protein = 0.2, fat = 0.2)

# Compute the convex hull of the available foods in EMT coordinates
h_idx <- chull(foods_xy$x, foods_xy$y)
h_idx <- c(h_idx, h_idx[1])

hull_df <- data.frame(
  x = foods_xy$x[h_idx],
  y = foods_xy$y[h_idx]
)

# ============================================================
# Decision rules for choosing food items
# ============================================================

# Angular alignment score between current state, a food item, and target
d_score <- function(O, Fi, I) {
  num <- sum((O - I) * (Fi - I))
  den <- sqrt(sum((O - I)^2)) * sqrt(sum((Fi - I)^2))
  if (den == 0) return(0)
  -num / den
}

# Exponential moving-average update toward a chosen food item
ema_step <- function(O, Fi, alpha = averagingRatio) {
  O_new <- (1 - alpha) * O + alpha * Fi
  O_new / sum(O_new)
}

# Sample a food item from a softmax distribution over scores
softmax_pick <- function(scores, k, w = NULL) {
  if (is.null(w)) w <- rep(1, length(scores))
  logits <- log(pmax(w, 1e-12)) + k * scores
  p <- exp(logits - max(logits))
  p <- p / sum(p)
  sample.int(length(scores), 1, prob = p)
}

# Weighted angular alignment score emphasizing selected nutrients
d_score_weighted <- function(O, Fi, I, W) {
  CI <- W * (O - I)
  num <- sum(CI * (Fi - I))
  den <- sqrt(sum(CI^2)) * sqrt(sum((Fi - I)^2))
  if (den == 0) return(0)
  -num / den
}

# ============================================================
# Simulation functions
# ============================================================

# Simulate a trajectory for a fixed number of steps under a food-choice rule
simulate_fixed <- function(O_start, chooser_fn, foods, steps = 50) {
  O <- O_start
  path <- matrix(NA_real_, nrow = steps + 1, ncol = 3)
  path[1, ] <- O
  chosen <- integer(steps)
  
  for (t in 1:steps) {
    j <- chooser_fn(O, foods)
    chosen[t] <- j
    O <- ema_step(O, foods[j, ], alpha = averagingRatio)
    path[t + 1, ] <- O
  }
  
  list(path = path, chosen = chosen)
}

# ============================================================
# Food-choice strategies
# ============================================================

# Choose uniformly at random among foods
chooser_random <- function() {
  function(O, foods) {
    sample.int(nrow(foods), 1)
  }
}

# Choose the food that minimizes squared distance to the target after one step
chooser_greedy <- function() {
  function(O, foods) {
    dists <- apply(foods, 1, function(Fi) {
      cand <- ema_step(O, Fi, alpha = averagingRatio)
      sum((cand - I)^2)
    })
    which.min(dists)
  }
}

# Prefer a particular food item while still weighting by directional alignment
chooser_food_bias <- function(prefer_idx = 3, k = 1,
                              pref_weights = c(0.025, 0.025, 0.9, 0.025, 0.025)) {
  function(O, foods) {
    d <- apply(foods, 1, function(Fi) d_score(O, Fi, I))
    softmax_pick(d, k = k, w = pref_weights)
  }
}

# Prefer foods that align with a nutrient-weighted version of the target direction
chooser_protein_bias <- function(k = 2, W = c(0.05, 0.90, 0.05)) {
  function(O, foods) {
    d <- apply(foods, 1, function(Fi) d_score_weighted(O, Fi, I, W))
    softmax_pick(d, k = k)
  }
}

# ============================================================
# Run simulations
# ============================================================

set.seed(42)

sim_random   <- simulate_fixed(O0, chooser_random(),       foods, steps = 50)
sim_greedy   <- simulate_fixed(O0, chooser_greedy(),       foods, steps = 50)
sim_foodbias <- simulate_fixed(O0, chooser_food_bias(),    foods, steps = 50)
sim_protbias <- simulate_fixed(O0, chooser_protein_bias(), foods, steps = 50)

# Convert simulated trajectories to EMT coordinates
path_df <- function(sim) {
  df <- as.data.frame(sim$path)
  names(df) <- c("carb", "protein", "fat")
  xy <- bary_to_xy(df$carb, df$protein, df$fat)
  cbind(xy, t = seq_len(nrow(xy)))
}

# Convert chosen food indices to a data frame for time-series plotting
choices_df <- function(sim) {
  data.frame(
    t = seq_along(sim$chosen),
    food_idx = sim$chosen
  )
}

p_random_path   <- path_df(sim_random)
p_greedy_path   <- path_df(sim_greedy)
p_foodbias_path <- path_df(sim_foodbias)
p_protbias_path <- path_df(sim_protbias)

c_random <- choices_df(sim_random)
c_greedy <- choices_df(sim_greedy)
c_food   <- choices_df(sim_foodbias)
c_prot   <- choices_df(sim_protbias)

# Map numeric food indices to roman numeral labels
idx_to_label <- c("i", "ii", "iii", "iv", "v")

force_levels <- function(df) {
  df$food_label <- factor(idx_to_label[df$food_idx], levels = idx_to_label)
  df
}

c_random <- force_levels(c_random)
c_greedy <- force_levels(c_greedy)
c_food   <- force_levels(c_food)
c_prot   <- force_levels(c_prot)

# ============================================================
# Plot colors and axis annotations
# ============================================================

BLACK <- "black"

# Protein = green (left), fat = blue (right), carb = red (bottom)
colProtein <- "#00A000"
colFat     <- "#0055FF"
colCarb    <- "#D00000"

# Colors used for the internal grid lines
axis_cols_inside <- c(
  carb = colProtein,
  protein = colCarb,
  fat = colFat
)

# Tick labels for each side of the triangle
tick_vals <- levs

ticks_bottom <- data.frame(
  bary_to_xy(1 - tick_vals, tick_vals, 0),
  lab = paste0(100 * tick_vals, "%"),
  col = colCarb,
  nudge_x = 0,
  nudge_y = -0.04,
  hjust = 0.5,
  vjust = 1
)

ticks_left <- data.frame(
  bary_to_xy(1 - tick_vals, 0, tick_vals),
  lab = paste0(100 * tick_vals, "%"),
  col = colProtein,
  nudge_x = -0.04,
  nudge_y = 0,
  hjust = 1,
  vjust = 0.5
)

ticks_right <- data.frame(
  bary_to_xy(0, 1 - tick_vals, tick_vals),
  lab = paste0(100 * tick_vals, "%"),
  col = colFat,
  nudge_x = 0.04,
  nudge_y = 0,
  hjust = 0,
  vjust = 0.5
)

# Axis title anchor points
prot_mid <- bary_to_xy(0.5, 0.0, 0.5)
fat_mid  <- bary_to_xy(0.0, 0.5, 0.5)
carb_mid <- bary_to_xy(0.5, 0.5, 0.0)

# Axis title offsets
PROT_DELTA <- 0.10
FAT_DELTA  <- 0.10

PROT_NUDGE_X <- -(sqrt3 / 2) * PROT_DELTA
PROT_NUDGE_Y <- 0.5 * PROT_DELTA

FAT_NUDGE_X <- (sqrt3 / 2) * FAT_DELTA
FAT_NUDGE_Y <- 0.5 * FAT_DELTA

CARB_NUDGE_Y <- -0.12

# Axis title data frame
axis_titles <- rbind(
  data.frame(
    x = prot_mid$x,
    y = prot_mid$y,
    lab = "Protein %",
    col = colProtein,
    hjust = 0.5,
    vjust = 0.5,
    nudge_x = PROT_NUDGE_X,
    nudge_y = PROT_NUDGE_Y,
    angle = 60
  ),
  data.frame(
    x = fat_mid$x,
    y = fat_mid$y,
    lab = "Fat %",
    col = colFat,
    hjust = 0.5,
    vjust = 0.5,
    nudge_x = FAT_NUDGE_X,
    nudge_y = FAT_NUDGE_Y,
    angle = -60
  ),
  data.frame(
    x = carb_mid$x,
    y = carb_mid$y,
    lab = "Carb %",
    col = colCarb,
    hjust = 0.5,
    vjust = 1,
    nudge_x = 0,
    nudge_y = CARB_NUDGE_Y,
    angle = 0
  )
)

# ============================================================
# Plot scaffolds
# ============================================================

tri_xmin <- -0.15
tri_xmax <- 1.15
tri_ymin <- -0.15
tri_ymax <- sqrt3 / 2 + 0.15

tri_x_range <- tri_xmax - tri_xmin
tri_y_range <- tri_ymax - tri_ymin
tri_height_width_ratio <- tri_y_range / tri_x_range

# Base EMT plot with triangle, grid, convex hull, and axis annotations
emt_base <- function() {
  ggplot() +
    geom_path(data = tri2d, aes(x, y), linewidth = 1.2, colour = BLACK) +
    geom_segment(data = grid_carb, aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 0.9, colour = axis_cols_inside["carb"], alpha = 1) +
    geom_segment(data = grid_prot, aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 0.9, colour = axis_cols_inside["protein"], alpha = 1) +
    geom_segment(data = grid_fat, aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 0.9, colour = axis_cols_inside["fat"], alpha = 1) +
    geom_polygon(data = hull_df, aes(x, y), fill = "grey70", alpha = 0.30, colour = NA) +
    geom_path(data = hull_df, aes(x, y), linewidth = 1.0, colour = "grey35") +
    geom_text(
      data = ticks_bottom,
      aes(x = x, y = y, label = lab),
      nudge_x = ticks_bottom$nudge_x,
      nudge_y = ticks_bottom$nudge_y,
      colour = ticks_bottom$col,
      fontface = "bold",
      size = 4.5,
      hjust = ticks_bottom$hjust,
      vjust = ticks_bottom$vjust,
      angle = 0
    ) +
    geom_text(
      data = ticks_left,
      aes(x = x, y = y, label = lab),
      nudge_x = ticks_left$nudge_x,
      nudge_y = ticks_left$nudge_y,
      colour = ticks_left$col,
      fontface = "bold",
      size = 4.5,
      hjust = ticks_left$hjust,
      vjust = ticks_left$vjust,
      angle = 60
    ) +
    geom_text(
      data = ticks_right,
      aes(x = x, y = y, label = lab),
      nudge_x = ticks_right$nudge_x,
      nudge_y = ticks_right$nudge_y,
      colour = ticks_right$col,
      fontface = "bold",
      size = 4.5,
      hjust = ticks_right$hjust,
      vjust = ticks_right$vjust,
      angle = -60
    ) +
    geom_text(
      data = axis_titles,
      aes(x = x, y = y, label = lab, angle = angle),
      nudge_x = axis_titles$nudge_x,
      nudge_y = axis_titles$nudge_y,
      colour = axis_titles$col,
      fontface = "bold",
      size = 6
    ) +
    coord_fixed(
      xlim = c(tri_xmin, tri_xmax),
      ylim = c(tri_ymin, tri_ymax),
      ratio = 1,
      expand = FALSE
    ) +
    theme_bw(base_size = 16) +
    theme(
      axis.title   = element_blank(),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      panel.grid   = element_blank(),
      panel.border = element_blank(),
      plot.title   = element_blank()
    )
}

# Points and labels overlaid on the EMT plot
emt_overlay_points <- function() {
  list(
    geom_point(data = foods_df, aes(x, y), shape = 16, size = 3.8, colour = BLACK),
    geom_text(data = foods_labels_df, aes(x, y, label = food),
              nudge_y = 0.055, fontface = "bold", size = 7.5),
    geom_point(aes(x = zbar_xy$x, y = zbar_xy$y),
               shape = 16, size = 4.6, colour = "steelblue4"),
    geom_point(aes(x = I_xy$x, y = I_xy$y),
               shape = 15, size = 4.8, colour = "goldenrod2")
  )
}

# Base time-series plot showing food choice across steps
ts_base <- function(df) {
  df$food_numeric <- as.numeric(df$food_label)
  
  x_min <- min(df$t)
  x_max <- max(df$t)
  x_range <- x_max - x_min
  
  y_min_raw <- min(df$food_numeric)
  y_max_raw <- max(df$food_numeric)
  y_min <- y_min_raw - 0.2
  y_max <- y_max_raw + 0.2
  y_range <- y_max - y_min
  
  ratio_ts <- tri_height_width_ratio * (x_range / y_range)
  
  ggplot(df, aes(x = t, y = food_numeric, group = 1)) +
    geom_step(linewidth = 1.1) +
    scale_y_continuous(
      breaks = 1:length(idx_to_label),
      labels = idx_to_label,
      expand = c(0, 0)
    ) +
    labs(x = "Step", y = "Food Item") +
    coord_fixed(
      ratio = ratio_ts,
      xlim = c(x_min, x_max),
      ylim = c(y_min, y_max),
      expand = FALSE
    ) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_blank(),
      axis.title = element_text(size = 20),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16, face = "bold")
    )
}

# ============================================================
# Assemble the four figure rows
# ============================================================

traj_col <- "#FF7F0E"

row_greedy <- (
  emt_base() +
    geom_path(data = p_greedy_path, aes(x, y), linewidth = 1.4, colour = traj_col, alpha = 0.98) +
    emt_overlay_points()
) | ts_base(c_greedy)
row_greedy <- row_greedy + plot_layout(widths = c(2, 3))

row_random <- (
  emt_base() +
    geom_path(data = p_random_path, aes(x, y), linewidth = 1.4, colour = traj_col, alpha = 0.98) +
    emt_overlay_points()
) | ts_base(c_random)
row_random <- row_random + plot_layout(widths = c(2, 3))

row_food <- (
  emt_base() +
    geom_path(data = p_foodbias_path, aes(x, y), linewidth = 1.4, colour = traj_col, alpha = 0.98) +
    emt_overlay_points()
) | ts_base(c_food)
row_food <- row_food + plot_layout(widths = c(2, 3))

row_prot <- (
  emt_base() +
    geom_path(data = p_protbias_path, aes(x, y), linewidth = 1.4, colour = traj_col, alpha = 0.98) +
    emt_overlay_points()
) | ts_base(c_prot)
row_prot <- row_prot + plot_layout(widths = c(2, 3))

# Common title styling
title_theme <- theme(
  plot.title = element_text(size = 26, face = "bold", hjust = 0.0)
)

# Add row titles
row_A <- row_greedy + plot_annotation(
  title = "A) Closest Distance Optimization",
  theme = title_theme
)

row_B <- row_random + plot_annotation(
  title = "B) Random Choice",
  theme = title_theme
)

row_C <- row_food + plot_annotation(
  title = "C) Food Item Bias",
  theme = title_theme
)

row_D <- row_prot + plot_annotation(
  title = "D) Macronutrient Bias",
  theme = title_theme
)

# ============================================================
# Display plots
# ============================================================

print(row_A)
print(row_B)
print(row_C)
print(row_D)