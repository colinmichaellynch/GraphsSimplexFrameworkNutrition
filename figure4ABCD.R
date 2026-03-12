# Clear the workspace
rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(fields)
  library(scales)
  library(grid)
})

# Define shared plotting constants
BLACK <- "black"
sqrt3 <- sqrt(3)

# ============================================================
# Helper functions for equilateral mixture triangle geometry
# ============================================================

# Convert barycentric coordinates (x1, x2, x3) to 2D Cartesian
# coordinates for an equilateral triangle
bary_to_xy <- function(x1, x2, x3) {
  x <- x2 + 0.5 * x3
  y <- (sqrt3 / 2) * x3
  data.frame(x, y)
}

# Convert 2D Cartesian coordinates back to barycentric coordinates
xy_to_bary <- function(x, y) {
  x3 <- (2 / sqrt3) * y
  x2 <- x - 0.5 * x3
  x1 <- 1 - x2 - x3
  cbind(x1 = x1, x2 = x2, x3 = x3)
}

# Determine whether a point lies inside the simplex
inside_simplex <- function(x, y, tol = 1e-9) {
  b <- xy_to_bary(x, y)
  apply(b >= -tol, 1, all)
}

# Define triangle vertices and centroid
v1 <- bary_to_xy(1, 0, 0)
v2 <- bary_to_xy(0, 1, 0)
v3 <- bary_to_xy(0, 0, 1)
tri2d <- rbind(v1, v2, v3, v1)
center <- bary_to_xy(1 / 3, 1 / 3, 1 / 3)

# ============================================================
# Build triangle scaffold, grid, and axis annotations
# ============================================================

# Build internal grid lines for specified mixture levels
build_grid_2d_levels <- function(levels = c(0.2, 0.4, 0.6, 0.8)) {
  mk <- function(p1, p2) {
    data.frame(x = p1$x, y = p1$y, xend = p2$x, yend = p2$y)
  }
  
  out <- list()
  for (t in levels) {
    out[[length(out) + 1]] <- mk(bary_to_xy(t, 0, 1 - t), bary_to_xy(t, 1 - t, 0))
    out[[length(out) + 1]] <- mk(bary_to_xy(0, t, 1 - t), bary_to_xy(1 - t, t, 0))
    out[[length(out) + 1]] <- mk(bary_to_xy(0, 1 - t, t), bary_to_xy(1 - t, 0, t))
  }
  
  do.call(rbind, out)
}

grid2d <- build_grid_2d_levels(c(0.2, 0.4, 0.6, 0.8))

# Compute the outward unit normal vector for a triangle edge
edge_outward_normal <- function(vA, vB) {
  e <- c(vB$x - vA$x, vB$y - vA$y)
  n1 <- c(-e[2], e[1])
  n1 <- n1 / sqrt(sum(n1^2))
  n2 <- -n1
  
  mid <- c((vA$x + vB$x) / 2, (vA$y + vB$y) / 2)
  to_center <- c(center$x - mid[1], center$y - mid[2])
  
  if (sum(n1 * to_center) < 0) n1 else n2
}

# Build percentage tick labels positioned outside a triangle edge
build_edge_tick_labels <- function(vA, vB, values, offset = 0.06) {
  n <- edge_outward_normal(vA, vB)
  
  do.call(rbind, lapply(values, function(t) {
    P <- c((1 - t) * vA$x + t * vB$x, (1 - t) * vA$y + t * vB$y)
    data.frame(
      x = P[1] + n[1] * offset,
      y = P[2] + n[2] * offset,
      txt = paste0(round(100 * t), "%")
    )
  }))
}

# Build an axis label centered outside a triangle edge
build_axis_label <- function(sym, vA, vB, offset = 0.10) {
  n <- edge_outward_normal(vA, vB)
  mid <- c((vA$x + vB$x) / 2, (vA$y + vB$y) / 2)
  
  data.frame(
    x = mid[1] + n[1] * offset,
    y = mid[2] + n[2] * offset,
    txt = sym
  )
}

ticks <- c(0.2, 0.4, 0.6, 0.8)

# Axis mapping: r1 on bottom, r2 on left, r3 on right
labs_r1 <- build_edge_tick_labels(v1, v2, ticks, offset = 0.06)
labs_r2 <- build_edge_tick_labels(v1, v3, ticks, offset = 0.06)
labs_r3 <- build_edge_tick_labels(v2, v3, ticks, offset = 0.06)

axis_labs <- rbind(
  build_axis_label("r[1]", v1, v2, offset = 0.10),
  build_axis_label("r[2]", v1, v3, offset = 0.10),
  build_axis_label("r[3]", v2, v3, offset = 0.10)
)

base_limits <- list(
  xlim = c(-0.15, 1.15),
  ylim = c(-0.15, sqrt3 / 2 + 0.15)
)

# ============================================================
# Generate a lattice of mixture compositions over the simplex
# ============================================================

# Build a dense triangular lattice over the simplex
simplex_grid <- function(m = 60) {
  L <- lapply(0:m, function(i) {
    lapply(0:(m - i), function(j) {
      k <- m - i - j
      data.frame(x1 = i / m, x2 = j / m, x3 = k / m)
    })
  })
  
  df <- do.call(rbind, unlist(L, recursive = FALSE))
  cbind(df, bary_to_xy(df$x1, df$x2, df$x3))
}

S <- simplex_grid(m = 60)

# ============================================================
# Construct synthetic response surfaces over mixture space
# ============================================================

# Offspring response surface
S$Offspring <- with(S, {
  p1 <- 1.25 * exp(-((x1 - 0.58)^2 / 0.025 + (x2 - 0.26)^2 / 0.035 + (x3 - 0.16)^2 / 0.035))
  p2 <- 0.28 * exp(-((x1 - 0.42)^2 / 0.055 + (x2 - 0.46)^2 / 0.055 + (x3 - 0.12)^2 / 0.055))
  base <- 0.07 + 0.09 * x1 - 0.04 * x3 - 0.08 * (x2 - 0.34)^2 + 0.03 * (x1 * x2)
  p1 + p2 + base
})

# Normalize offspring response for later use
o_rng <- range(S$Offspring)
S$O_norm <- (S$Offspring - o_rng[1]) / (diff(o_rng) + 1e-12)

# Body size response surface
S$BodySize <- with(S, {
  p1 <- 0.98 * exp(-((x1 - 0.22)^2 / 0.012 + (x2 - 0.68)^2 / 0.034 + (x3 - 0.10)^2 / 0.034))
  ridge <- 0.20 * exp(-((x1 - (0.72 * x2 + 0.07))^2) / 0.008) * exp(-((x3 - 0.22)^2) / 0.07)
  base <- 0.05 + 0.08 * x2 - 0.035 * (x1 - 0.25)^2 + 0.02 * (x2 * x3)
  p1 + ridge + base
})

# Lifespan response surface
S$Lifespan <- with(S, {
  p1 <- 1.15 * exp(-((x1 - 0.14)^2 / 0.050 + (x2 - 0.24)^2 / 0.050 + (x3 - 0.64)^2 / 0.016))
  p2 <- 0.16 * exp(-((x1 - 0.12)^2 / 0.015 + (x2 - 0.63)^2 / 0.040 + (x3 - 0.25)^2 / 0.040))
  base <- 0.07 + 0.10 * x3 - 0.05 * (x1 - 0.20)^2 + 0.02 * (x1 * x3)
  dip_left <- 0.35 * exp(-((x1 - 0.025)^2 / 0.1 + (x2 - 0.98)^2 / 0.1 + (x3 - 0.55)^2 / 0.1))
  (p1 + p2 + base) - 0.44 * S$O_norm - dip_left
})

# Assay performance response surface
S$AssayPerf <- with(S, {
  p1 <- 1.02 * exp(-((x1 - 0.36)^2 / 0.024 + (x2 - 0.40)^2 / 0.024 + (x3 - 0.24)^2 / 0.024))
  synergy <- 0.30 * (x1 * x2 + x2 * x3 + x1 * x3)
  base <- 0.05 - 0.028 * (x1 - 0.63)^2 - 0.028 * (x2 - 0.22)^2 + 0.02 * (x1 * x3)
  p1 + synergy + base
})

# ============================================================
# Compute Derringer-Suich desirability scores
# ============================================================

# Larger-is-better desirability function
d_max <- function(y, L, T, r = 1) {
  d <- ifelse(
    y <= L, 0,
    ifelse(y >= T, 1, ((y - L) / (T - L))^r)
  )
  pmin(pmax(d, 0), 1)
}

# Return lower and upper quantiles for desirability thresholds
qfun <- function(v) {
  quantile(v, probs = c(0.10, 0.90), na.rm = TRUE)
}

# Compute desirability for each response
S$d_Off <- d_max(S$Offspring, qfun(S$Offspring)[1], qfun(S$Offspring)[2], r = 0.9)
S$d_Bod <- d_max(S$BodySize,  qfun(S$BodySize)[1],  qfun(S$BodySize)[2],  r = 0.9)
S$d_Lif <- d_max(S$Lifespan,  qfun(S$Lifespan)[1],  qfun(S$Lifespan)[2],  r = 0.9)
S$d_Ass <- d_max(S$AssayPerf, qfun(S$AssayPerf)[1], qfun(S$AssayPerf)[2], r = 0.9)

# Weighted geometric mean of individual desirability scores
w <- c(Off = 1.1, Bod = 0.9, Lif = 1.1, Ass = 0.9)
S$Desirability <- (S$d_Off^w["Off"] * S$d_Bod^w["Bod"] * S$d_Lif^w["Lif"] * S$d_Ass^w["Ass"])^(1 / sum(w))

# ============================================================
# Build a rectangular prediction grid for TPS smoothing
# ============================================================

# Construct a rectangular grid and mask it to the simplex
make_rect_grid <- function(nx = 300, ny = 300) {
  xlo <- -0.01
  xhi <- 1.01
  ylo <- -0.01
  yhi <- sqrt3 / 2 + 0.01
  
  xo <- seq(xlo, xhi, length.out = nx)
  yo <- seq(ylo, yhi, length.out = ny)
  grd <- expand.grid(x = xo, y = yo)
  keep <- inside_simplex(grd$x, grd$y, tol = 1e-9)
  
  list(grid = grd, keep = keep)
}

RECT <- make_rect_grid(nx = 300, ny = 300)

# ============================================================
# Thin-plate spline smoothing and prediction
# ============================================================

# Estimate a range parameter for fastTps from pairwise distances
.est_aRange <- function(X, max_n = 1200) {
  n <- nrow(X)
  
  if (n > max_n) {
    set.seed(123)
    X <- X[sample(n, max_n), , drop = FALSE]
  }
  
  d <- fields::rdist(X)
  stats::median(d[upper.tri(d)], na.rm = TRUE)
}

# Fit a TPS model and predict it on the rectangular grid
tps_to_rect <- function(df, zcol, RECT, prefer_fast = TRUE) {
  XY <- cbind(df$x, df$y)
  z <- df[[zcol]]
  
  dup <- !duplicated(XY)
  XYu <- XY[dup, , drop = FALSE]
  zu <- z[dup]
  
  fit <- NULL
  if (prefer_fast && "fastTps" %in% getNamespaceExports("fields")) {
    ar <- .est_aRange(XYu)
    fit <- tryCatch(
      fields::fastTps(XYu, zu, aRange = ar),
      error = function(e) NULL
    )
  }
  
  if (is.null(fit)) {
    fit <- fields::Tps(XYu, zu)
  }
  
  grd <- RECT$grid
  grd$z <- as.vector(predict(fit, grd))
  grd[RECT$keep, , drop = FALSE]
}

# ============================================================
# Plot helper for smoothed triangular surfaces
# ============================================================

# Plot a smoothed surface over the mixture triangle
tri_surface <- function(df, zcol, title, palette = "Viridis",
                        bins = 22, alpha_fill = 0.70, pal_values = NULL) {
  if (is.null(pal_values)) {
    pal <- grDevices::hcl.colors(bins, palette = palette, rev = FALSE)
  } else {
    pal <- grDevices::colorRampPalette(pal_values)(bins)
  }
  
  Z <- tps_to_rect(df, zcol, RECT)
  
  ggplot() +
    stat_contour_filled(
      data = Z,
      inherit.aes = FALSE,
      aes(x = x, y = y, z = z, fill = after_stat(level)),
      bins = bins,
      alpha = alpha_fill
    ) +
    scale_fill_manual(values = pal, drop = FALSE) +
    geom_path(data = tri2d, aes(x, y), linewidth = 1, colour = BLACK) +
    geom_segment(
      data = grid2d,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.5,
      alpha = 0.9,
      colour = BLACK
    ) +
    geom_text(data = labs_r1, aes(x, y, label = txt), size = 3.6, colour = BLACK) +
    geom_text(data = labs_r2, aes(x, y, label = txt), size = 3.6, colour = BLACK) +
    geom_text(data = labs_r3, aes(x, y, label = txt), size = 3.6, colour = BLACK) +
    geom_text(data = axis_labs, aes(x, y, label = txt), parse = TRUE, size = 5, colour = BLACK) +
    coord_fixed(xlim = base_limits$xlim, ylim = base_limits$ylim, expand = FALSE) +
    theme_bw(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.position = "none"
    ) +
    labs(title = title)
}

# ============================================================
# Plot individual response surfaces
# ============================================================

p_off <- tri_surface(S, "Offspring", "A) # Offspring", palette = "YlOrRd", bins = 22, alpha_fill = 0.72)
p_bod <- tri_surface(S, "BodySize", "B) Body Size", palette = "BluGrn", bins = 22, alpha_fill = 0.72)
p_lif <- tri_surface(S, "Lifespan", "C) Lifespan", palette = "Plasma", bins = 22, alpha_fill = 0.72)
p_ass <- tri_surface(S, "AssayPerf", "D) Assay Performance", palette = "Viridis", bins = 22, alpha_fill = 0.72)

# Arrange the four response surfaces into a 2x2 layout
p_four <- (p_off | p_bod) / (p_lif | p_ass)

# ============================================================
# Plot composite desirability surface
# ============================================================

# Custom color palette for the composite desirability surface
des_pal_base <- c(
  "#0B1026", "#0E3A5A", "#1E6C8C", "#4BA3C3", "#CFE9E7",
  "#F2E9C9", "#EEB44F", "#E1792C", "#C13B2A", "#8A2BE2"
)

# Identify the point of maximal desirability on the smoothed surface
Z_des <- tps_to_rect(S, "Desirability", RECT)
imax <- which.max(Z_des$z)
star_df <- data.frame(x = Z_des$x[imax], y = Z_des$y[imax])
gold <- "#41DC8E"

# Plot composite desirability and mark the optimum
p_des <- tri_surface(
  S,
  "Desirability",
  "E) Composite Metric",
  bins = 26,
  alpha_fill = 0.78,
  pal_values = des_pal_base
) +
  geom_point(
    data = star_df,
    aes(x, y),
    shape = 8,
    size = 5.2,
    stroke = 1.2,
    colour = gold
  )

# ============================================================
# Identify the Pareto frontier for offspring and lifespan
# ============================================================

# Determine whether a point is dominated in the two-objective space
is_dominated <- function(A, B, Aall, Ball) {
  any((Aall >= A & Ball >= B) & (Aall > A | Ball > B))
}

# Extract nondominated points
nd_idx <- which(!mapply(
  is_dominated,
  S$Offspring, S$Lifespan,
  MoreArgs = list(Aall = S$Offspring, Ball = S$Lifespan)
))

Pareto <- S[nd_idx, , drop = FALSE]
Pareto_os <- Pareto[order(Pareto$Offspring, Pareto$Lifespan), ]

# Define colors for Pareto visualizations
pareto_col <- "#7B4CEB"
offspring_col <- "#D95F0E"
lifespan_col <- "#E91E63"

# Build a convex hull around the Pareto set in mixture space
if (nrow(Pareto) >= 3) {
  hull_idx <- grDevices::chull(Pareto$x, Pareto$y)
  Pareto_hull <- Pareto[hull_idx, , drop = FALSE]
} else {
  Pareto_hull <- Pareto
}

# Split Pareto points into low- and high-offspring halves for annotation
med_off <- median(Pareto$Offspring, na.rm = TRUE)
P_lo <- Pareto[Pareto$Offspring <= med_off, , drop = FALSE]
P_hi <- Pareto[Pareto$Offspring >= med_off, , drop = FALSE]

cent_lo <- data.frame(x = mean(P_lo$x), y = mean(P_lo$y))
cent_hi <- data.frame(x = mean(P_hi$x), y = mean(P_hi$y))

# ============================================================
# Plot Pareto region on the mixture triangle
# ============================================================

p_frontier_on_triangle <- ggplot() +
  geom_polygon(
    data = Pareto_hull,
    aes(x = x, y = y),
    fill = alpha(pareto_col, 0.28),
    colour = pareto_col,
    linewidth = 1.1
  ) +
  geom_path(data = tri2d, aes(x, y), linewidth = 1, colour = BLACK) +
  geom_segment(
    data = grid2d,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5,
    alpha = 0.9,
    colour = BLACK
  ) +
  geom_label(
    data = cent_hi,
    aes(x, y),
    label = "↑ # Offspring",
    fontface = "bold",
    size = 4.3,
    colour = offspring_col,
    fill = "#FFFFFF",
    alpha = 1,
    label.size = 0.6,
    label.r = unit(2, "pt"),
    label.padding = unit(3, "pt")
  ) +
  geom_label(
    data = cent_lo,
    aes(x, y),
    label = "↑ Lifespan",
    fontface = "bold",
    size = 4.3,
    colour = lifespan_col,
    fill = "#FFFFFF",
    alpha = 1,
    label.size = 0.6,
    label.r = unit(2, "pt"),
    label.padding = unit(3, "pt")
  ) +
  geom_text(data = labs_r1, aes(x, y, label = txt), size = 3.6, colour = BLACK) +
  geom_text(data = labs_r2, aes(x, y, label = txt), size = 3.6, colour = BLACK) +
  geom_text(data = labs_r3, aes(x, y, label = txt), size = 3.6, colour = BLACK) +
  geom_text(data = axis_labs, aes(x, y, label = txt), parse = TRUE, size = 5, colour = BLACK) +
  coord_fixed(xlim = base_limits$xlim, ylim = base_limits$ylim, expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  ) +
  labs(title = "F)")

# ============================================================
# Plot Pareto frontier in objective space
# ============================================================

S_all <- transform(S, grp = "Other points")
Pareto_lab <- transform(Pareto_os, grp = "Pareto frontier")

p_pf <- ggplot() +
  geom_point(
    data = S_all,
    aes(x = Offspring, y = Lifespan, colour = grp),
    alpha = 0.35,
    size = 0.8
  ) +
  geom_path(
    data = Pareto_lab,
    aes(x = Offspring, y = Lifespan, colour = grp),
    linewidth = 1.2
  ) +
  geom_point(
    data = Pareto_lab,
    aes(x = Offspring, y = Lifespan, colour = grp),
    size = 1.8
  ) +
  scale_colour_manual(
    name = NULL,
    values = c("Other points" = "gray70", "Pareto frontier" = pareto_col)
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.key = element_rect(fill = "white"),
    legend.margin = margin(4, 6, 4, 6),
    legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "G)", x = "# Offspring", y = "Lifespan") +
  coord_equal()

# ============================================================
# Display plots
# ============================================================

p_four
p_des
p_frontier_on_triangle
p_pf
