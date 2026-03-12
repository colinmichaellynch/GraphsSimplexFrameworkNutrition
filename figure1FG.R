# Clear the workspace
rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# Define shared plotting constants
set.seed(7)
BLACK <- "black"
sqrt3 <- sqrt(3)

# ============================================================
# Helper functions for equilateral mixture triangle geometry
# ============================================================

# Convert barycentric coordinates (z1, z2, z3) to 2D Cartesian
# coordinates for an equilateral triangle
bary_to_xy <- function(z1, z2, z3) {
  x <- z2 + 0.5 * z3
  y <- (sqrt3 / 2) * z3
  data.frame(x, y)
}

# Define triangle vertices and centroid
v1 <- bary_to_xy(1, 0, 0)
v2 <- bary_to_xy(0, 1, 0)
v3 <- bary_to_xy(0, 0, 1)
tri2d <- rbind(v1, v2, v3, v1)
center <- bary_to_xy(1 / 3, 1 / 3, 1 / 3)

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

# Build percentage labels positioned outside a triangle edge
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

# ============================================================
# Right-angled mixture triangle (RMT) scaffold
# ============================================================

# Triangle outline in right-angled coordinates
rmt_tri <- data.frame(
  z1 = c(0, 1, 0, 0),
  z2 = c(0, 0, 1, 0)
)

# Build internal grid lines for the RMT
build_rmt_grid <- function(levels = c(0.2, 0.4, 0.6, 0.8)) {
  out <- list()
  
  for (t in levels) {
    out[[length(out) + 1]] <- data.frame(x = t, y = 0,     xend = t,     yend = 1 - t)
    out[[length(out) + 1]] <- data.frame(x = 0, y = t,     xend = 1 - t, yend = t)
    out[[length(out) + 1]] <- data.frame(x = 0, y = 1 - t, xend = 1 - t, yend = 0)
  }
  
  do.call(rbind, out)
}

rmt_grid <- build_rmt_grid(ticks)

# Percentage labels for the bottom and left edges of the RMT
rmt_edge_ticks <- function(values, which = c("z1", "z2"), offset = 0.06) {
  which <- match.arg(which)
  
  if (which == "z1") {
    data.frame(x = values, y = 0 - offset, txt = paste0(round(100 * values), "%"))
  } else {
    data.frame(x = 0 - offset, y = values, txt = paste0(round(100 * values), "%"))
  }
}

rmt_labs_z1 <- rmt_edge_ticks(ticks, "z1", offset = 0.06)
rmt_labs_z2 <- rmt_edge_ticks(ticks, "z2", offset = 0.06)

# Percentage labels for the diagonal r3 direction
rmt_bisector_ticks_z3 <- function(values, advance = 0.015) {
  do.call(rbind, lapply(values, function(t) {
    pos <- (1 - t) / 2
    data.frame(
      x = pos + advance,
      y = pos + advance,
      txt = paste0(round(100 * t), "%")
    )
  }))
}

rmt_labs_z3 <- rmt_bisector_ticks_z3(ticks, advance = 0.015)

# Axis labels for the RMT
rmt_axis_labels <- {
  mid_bottom <- c(0.5, 0.0)
  mid_left   <- c(0.0, 0.5)
  mid_hyp    <- c(0.5, 0.5)
  
  n_bottom <- c(0, -1)
  n_left   <- c(-1, 0)
  n_hyp    <- c(1, 1) / sqrt(2)
  
  off <- 0.10
  
  data.frame(
    x = c(
      mid_bottom[1] + off * n_bottom[1],
      mid_left[1]   + off * n_left[1],
      mid_hyp[1]    + off * n_hyp[1]
    ),
    y = c(
      mid_bottom[2] + off * n_bottom[2],
      mid_left[2]   + off * n_left[2],
      mid_hyp[2]    + off * n_hyp[2]
    ),
    txt = c("r[1]", "r[2]", "r[3]")
  )
}

# ============================================================
# Functions for generating anisotropic clusters in simplex space
# ============================================================

# Normalize a vector to unit length
unit <- function(v) {
  v / sqrt(sum(v^2))
}

# Construct an orthonormal basis for a plane in simplex space
make_plane_basis <- function(axis_dir) {
  n <- c(1, 1, 1)
  u1 <- unit(axis_dir)
  
  w0 <- c(1, -1, 0)
  w <- w0 - sum(w0 * n) / sum(n * n) * n
  w <- w - sum(w * u1) * u1
  u2 <- unit(w)
  
  list(u1 = u1, u2 = u2)
}

# Sample a cluster of points centered at mu, elongated along axis_dir,
# and constrained to remain inside the simplex
sample_group <- function(mu, axis_dir, n = 450, sigma_hi = 0.06, sigma_lo = 0.02, max_tries = 30000) {
  B <- make_plane_basis(axis_dir)
  u1 <- B$u1
  u2 <- B$u2
  
  out <- matrix(NA_real_, nrow = n, ncol = 3)
  got <- 0
  tries <- 0
  
  while (got < n && tries < max_tries) {
    a <- rnorm(1, 0, sigma_hi)
    b <- rnorm(1, 0, sigma_lo)
    z <- mu + a * u1 + b * u2
    
    if (all(z > 0) && all(z < 1)) {
      z <- z / sum(z)
      if (all(z > 0) && all(z < 1)) {
        got <- got + 1
        out[got, ] <- z
      }
    }
    
    tries <- tries + 1
  }
  
  if (got < n) {
    warning(sprintf(
      "Only obtained %d / %d points; consider smaller variances or moving the mean inward.",
      got, n
    ))
    out <- out[1:got, , drop = FALSE]
  }
  
  as.data.frame(out) |> setNames(c("z1", "z2", "z3"))
}

# ============================================================
# Generate three groups with different dominant directions
# ============================================================

# Group means in simplex space
muA <- c(0.55, 0.30, 0.15)
muB <- c(0.20, 0.60, 0.20)
muC <- c(0.25, 0.25, 0.50)

# Dispersion parameters
sigma_hi <- 0.06
sigma_lo <- 0.02
squish_boost <- 1

# Dominant elongation directions for each group
axis_A <- c(1, 0, -1)
axis_B <- c(0, 1, -1)
axis_C <- c(-1 / 2, -1 / 2, 1)

# Sample the three groups
A <- sample_group(muA, axis_A, n = 450, sigma_hi = sigma_hi, sigma_lo = sigma_lo) |>
  mutate(group = "A: hi || r[1]")

B <- sample_group(muB, axis_B, n = 450, sigma_hi = sigma_hi, sigma_lo = sigma_lo) |>
  mutate(group = "B: hi || r[2]")

C <- sample_group(muC, axis_C, n = 450, sigma_hi = sigma_hi * squish_boost, sigma_lo = sigma_lo) |>
  mutate(group = "C: hi along r[3]")

dat <- bind_rows(A, B, C)

# ============================================================
# Build EMT coordinates with r1 and r2 swapped
# ============================================================

# Swap the first two simplex components before projection into EMT space
tmp <- dat |> transmute(group, r1s = z2, r2s = z1, r3s = z3)
dat_emt <- bind_cols(tmp, bary_to_xy(tmp$r1s, tmp$r2s, tmp$r3s))

# Build the EMT grid after swapping r1 and r2
build_grid_2d_levels_swapped <- function(levels = c(0.2, 0.4, 0.6, 0.8)) {
  mk <- function(p1, p2) {
    data.frame(x = p1$x, y = p1$y, xend = p2$x, yend = p2$y)
  }
  
  out <- list()
  for (t in levels) {
    out[[length(out) + 1]] <- mk(bary_to_xy(0,     t,     1 - t), bary_to_xy(1 - t, t,     0))
    out[[length(out) + 1]] <- mk(bary_to_xy(t,     0,     1 - t), bary_to_xy(t,     1 - t, 0))
    out[[length(out) + 1]] <- mk(bary_to_xy(0,     1 - t, t),     bary_to_xy(1 - t, 0,     t))
  }
  
  do.call(rbind, out)
}

grid2d_emt <- build_grid_2d_levels_swapped(ticks)

# Axis labels for the EMT after swapping r1 and r2
axis_labs_emt <- rbind(
  build_axis_label("r[1]", v1, v2, offset = 0.10),
  build_axis_label("r[2]", v3, v1, offset = 0.10),
  build_axis_label("r[3]", v2, v3, offset = 0.10)
)

# Percentage labels along each EMT edge
labs_bottom_emt <- build_edge_tick_labels(v1, v2, ticks, offset = 0.06)
labs_left_emt   <- build_edge_tick_labels(v3, v1, ticks, offset = 0.06)
labs_right_emt  <- build_edge_tick_labels(v2, v3, ticks, offset = 0.06)

# ============================================================
# Plot the equilateral mixture triangle
# ============================================================

p_emt <- ggplot() +
  geom_path(data = tri2d, aes(x, y), linewidth = 1, colour = BLACK) +
  geom_segment(
    data = grid2d_emt,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5,
    alpha = 0.8,
    colour = BLACK
  ) +
  geom_text(data = labs_bottom_emt, aes(x, y, label = txt), size = 3.4, colour = BLACK) +
  geom_text(data = labs_left_emt,   aes(x, y, label = txt), size = 3.4, colour = BLACK) +
  geom_text(data = labs_right_emt,  aes(x, y, label = txt), size = 3.4, colour = BLACK) +
  geom_text(data = axis_labs_emt,   aes(x, y, label = txt), parse = TRUE, size = 5, colour = BLACK) +
  geom_point(data = dat_emt, aes(x, y, color = group), alpha = 0.55, size = 1.5, show.legend = FALSE) +
  coord_fixed(
    xlim = c(-0.15, 1.15),
    ylim = c(-0.15, sqrt3 / 2 + 0.15),
    expand = FALSE
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title      = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_blank(),
    legend.position = "none",
    plot.title      = element_blank()
  )

# ============================================================
# Plot the right-angled mixture triangle
# ============================================================

p_rmt <- ggplot() +
  geom_path(data = rmt_tri, aes(z1, z2), linewidth = 1, colour = BLACK) +
  geom_segment(
    data = rmt_grid,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5,
    alpha = 0.8,
    colour = BLACK
  ) +
  geom_point(data = dat, aes(z1, z2, color = group), alpha = 0.55, size = 1.5) +
  geom_text(data = rmt_labs_z1, aes(x, y, label = txt), size = 3.4, colour = BLACK) +
  geom_text(data = rmt_labs_z2, aes(x, y, label = txt), size = 3.4, colour = BLACK) +
  geom_text(
    data = rmt_labs_z3,
    aes(x, y, label = txt),
    size = 3.6,
    colour = BLACK,
    angle = 45,
    hjust = 0,
    vjust = 0.5
  ) +
  geom_text(data = rmt_axis_labels, aes(x, y, label = txt), parse = TRUE, size = 5, colour = BLACK) +
  coord_equal(
    xlim = c(-0.20, 1.08),
    ylim = c(-0.12, 1.08),
    expand = FALSE
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title      = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_blank(),
    legend.position = "none",
    plot.title      = element_blank(),
    plot.margin     = margin(10, 20, 10, 45)
  )

# ============================================================
# Display plots side by side
# ============================================================

p_emt + p_rmt

