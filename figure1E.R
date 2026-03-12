# Clear the workspace
rm(list = ls())

# Load required packages
set.seed(7)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)
})

# Define shared plotting constants
BLACK <- "black"

# ============================================================
# Helper functions for EMT geometry
# ============================================================

# Convert barycentric coordinates (z1, z2, z3) to 2D Cartesian
# coordinates for an equilateral mixture triangle
bary_to_xy <- function(z1, z2, z3) {
  x <- z2 + 0.5 * z3
  y <- (sqrt(3) / 2) * z3
  data.frame(x, y)
}

# Define EMT vertices and triangle outline
v1 <- bary_to_xy(1, 0, 0)
v2 <- bary_to_xy(0, 1, 0)
v3 <- bary_to_xy(0, 0, 1)
tri2d <- rbind(v1, v2, v3, v1)

# Build the internal EMT grid
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

grid2d <- build_grid_2d_levels()

# ============================================================
# Build EMT axis labels and tick labels
# ============================================================

# Define EMT centroid
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

# Map EMT axes to triangle edges
edges <- list(
  z1 = list(A = v2, B = v3),
  z2 = list(A = v3, B = v1),
  z3 = list(A = v1, B = v2)
)

# Build an axis label centered outside a triangle edge
build_axis_label <- function(sym, vA, vB, offset = 0.16) {
  n <- edge_outward_normal(vA, vB)
  mid <- c((vA$x + vB$x) / 2, (vA$y + vB$y) / 2)
  
  data.frame(
    x = mid[1] + n[1] * offset,
    y = mid[2] + n[2] * offset,
    txt = sym
  )
}

ticks <- c(0.2, 0.4, 0.6, 0.8)

# EMT axis labels
axis_labs <- rbind(
  build_axis_label("r[3]", edges$z1$A, edges$z1$B, 0.16),
  build_axis_label("r[2]", edges$z2$A, edges$z2$B, 0.16),
  build_axis_label("r[1]", edges$z3$A, edges$z3$B, 0.16)
)

# Interpolate a point along a triangle edge
point_on_edge <- function(vA, vB, s) {
  data.frame(
    x = (1 - s) * vA$x + s * vB$x,
    y = (1 - s) * vA$y + s * vB$y
  )
}

# Build EMT tick labels outside a triangle edge
make_edge_ticks <- function(vA, vB, offset = 0.08) {
  n <- edge_outward_normal(vA, vB)
  
  do.call(rbind, lapply(ticks, function(s) {
    p <- point_on_edge(vA, vB, s)
    data.frame(
      x = p$x + n[1] * offset,
      y = p$y + n[2] * offset,
      txt = sprintf("%.1f", s)
    )
  }))
}

# EMT tick labels shown in Step 1
emt_ticks_all <- rbind(
  make_edge_ticks(edges$z3$A, edges$z3$B, 0.08),
  make_edge_ticks(edges$z2$A, edges$z2$B, 0.08),
  make_edge_ticks(edges$z1$A, edges$z1$B, 0.08)
) %>%
  mutate(step = factor("Step 1", levels = paste0("Step ", 1:6)))

# ============================================================
# Build the right-angled mixture triangle scaffold
# ============================================================

# Triangle outline in RMT coordinates
rmt_tri <- data.frame(
  z1 = c(0, 1, 0, 0),
  z2 = c(0, 0, 1, 0)
)

# Build the internal RMT grid
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

# Build RMT axis label positions
mid_bottom <- c(0.5, 0.0)
mid_left   <- c(0.0, 0.5)
mid_hyp    <- c(0.5, 0.5)

n_bottom <- c(0, -1)
n_left   <- c(-1, 0)
n_hyp    <- c(1, 1) / sqrt(2)

off <- 0.10

rmt_axis_labels <- data.frame(
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

# ============================================================
# Transform the RMT into EMT-scaled coordinates
# ============================================================

x_off <- 1.3
y_scale <- sqrt(3) / 2

# Shift and scale RMT coordinates so they can be interpolated from EMT
shift_scale <- function(df, kind = c("grid", "tri", "axislab")) {
  kind <- match.arg(kind)
  
  if (kind == "grid") {
    transform(df, xs = x + x_off, ys = y * y_scale, xend = xend + x_off, yend = yend * y_scale)
  } else if (kind == "tri") {
    transform(df, xs = z1 + x_off, ys = z2 * y_scale)
  } else {
    transform(df, xs = x + x_off, ys = y * y_scale)
  }
}

rmt_tri_ss  <- shift_scale(rmt_tri, "tri")
rmt_grid_ss <- shift_scale(rmt_grid, "grid")
rmt_ax_ss   <- shift_scale(rmt_axis_labels, "axislab")

# ============================================================
# Interpolate the morph from EMT to RMT
# ============================================================

# Interpolate triangle outlines between EMT and transformed RMT
interp_tri <- function(t) {
  X <- as.matrix(tri2d[, c("x", "y")])
  Y <- as.matrix(rmt_tri_ss[, c("xs", "ys")])
  Z <- (1 - t) * X + t * Y
  data.frame(x = Z[, 1], y = Z[, 2], t = t)
}

# Interpolate internal grids between EMT and transformed RMT
interp_grid <- function(t) {
  Xs <- as.matrix(grid2d[, c("x", "y", "xend", "yend")])
  Ys <- as.matrix(rmt_grid_ss[, c("xs", "ys", "xend", "yend")])
  Zs <- (1 - t) * Xs + t * Ys
  data.frame(x = Zs[, 1], y = Zs[, 2], xend = Zs[, 3], yend = Zs[, 4], t = t)
}

# Define interpolation steps
t_seq <- c(0, 0.25, 0.5, 0.75, 1)

tri_frames <- do.call(rbind, lapply(t_seq, interp_tri))
grid_frames <- do.call(rbind, lapply(t_seq, interp_grid))

tri_frames$step  <- factor(paste0("Step ", match(tri_frames$t, t_seq)), levels = paste0("Step ", 1:6))
grid_frames$step <- factor(paste0("Step ", match(grid_frames$t, t_seq)), levels = paste0("Step ", 1:6))

# ============================================================
# Define axis labels for EMT end state and RMT end state
# ============================================================

emt_ax_end <- axis_labs %>%
  mutate(step = factor("Step 1", levels = paste0("Step ", 1:6)))

rmt_ax_end <- data.frame(x = rmt_ax_ss$xs, y = rmt_ax_ss$ys, txt = rmt_ax_ss$txt) %>%
  mutate(step = factor("Step 5", levels = paste0("Step ", 1:6)))

# Reposition the r2 label for Step 5
ax_step5_z2_top <- data.frame(
  x = x_off - 0.02,
  y = y_scale + 0.06,
  txt = "r[2]",
  step = factor("Step 5", levels = paste0("Step ", 1:6))
)

rmt_ax_end <- rmt_ax_end %>%
  dplyr::filter(!(step == "Step 5" & txt == "r[2]")) %>%
  bind_rows(ax_step5_z2_top)

# ============================================================
# Define the final Step 6 RMT view
# ============================================================

tri_step6 <- transform(rmt_tri_ss, x = xs, y = ys)[, c("x", "y")]
grid_step6 <- transform(rmt_grid_ss, x = xs, y = ys, xend = xend, yend = yend)[, c("x", "y", "xend", "yend")]

ax_step6 <- data.frame(
  x = c(
    rmt_ax_ss$xs[rmt_axis_labels$txt == "r[1]"],
    x_off - 0.02,
    rmt_ax_ss$xs[rmt_axis_labels$txt == "r[3]"]
  ),
  y = c(
    rmt_ax_ss$ys[rmt_axis_labels$txt == "r[1]"],
    y_scale + 0.06,
    rmt_ax_ss$ys[rmt_axis_labels$txt == "r[3]"]
  ),
  txt = c("r[1]", "r[2]", "r[3]")
)

tri_step6$step  <- factor("Step 6", levels = paste0("Step ", 1:6))
grid_step6$step <- factor("Step 6", levels = paste0("Step ", 1:6))
ax_step6$step   <- factor("Step 6", levels = paste0("Step ", 1:6))

# ============================================================
# Define auxiliary annotations for Steps 5 and 6
# ============================================================

# Step 5 annotation settings
left_tick_dx  <- 0.10
left_arrow_dx <- 0.22

# Step 6 annotation settings
right_tick_dx  <- -0.10
right_arrow_dx <- -0.22
right_text_dx  <- -0.22

tick_y <- ticks * y_scale

# Step 5: labels and arrow for the "before" scale
before_labels_5 <- data.frame(
  x = x_off - left_tick_dx,
  y = tick_y,
  txt = sprintf("%.1f", 1 - ticks),
  step = factor("Step 5", levels = paste0("Step ", 1:6))
)

arrow_before_5 <- data.frame(
  x = x_off - left_arrow_dx,
  y = 0.80 * y_scale,
  xend = x_off - left_arrow_dx,
  yend = 0.18 * y_scale,
  step = factor("Step 5", levels = paste0("Step ", 1:6))
)

label_before_5 <- transform(
  arrow_before_5,
  lx = x,
  ly = (y + yend) / 2,
  lbl = "before"
)

# Step 6: labels and arrow for the "after" scale
after_labels_6 <- data.frame(
  x = x_off + right_tick_dx,
  y = tick_y,
  txt = sprintf("%.1f", ticks),
  step = factor("Step 6", levels = paste0("Step ", 1:6))
)

arrow_after_6 <- data.frame(
  x = x_off + right_arrow_dx,
  y = 0.18 * y_scale,
  xend = x_off + right_arrow_dx,
  yend = 0.80 * y_scale,
  step = factor("Step 6", levels = paste0("Step ", 1:6))
)

label_after_6 <- data.frame(
  lx = x_off + right_text_dx,
  ly = (0.18 * y_scale + 0.80 * y_scale) / 2,
  lbl = "after",
  step = factor("Step 6", levels = paste0("Step ", 1:6))
)

# ============================================================
# Combine all frame data
# ============================================================

tri_frames <- bind_rows(tri_frames, tri_step6)
grid_frames <- bind_rows(grid_frames, grid_step6)
axis_end <- bind_rows(emt_ax_end, rmt_ax_end, ax_step6)

# ============================================================
# Build the faceted morph plot
# ============================================================

p <- ggplot() +
  geom_segment(
    data = grid_frames,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5,
    alpha = 0.85,
    colour = BLACK
  ) +
  geom_path(
    data = tri_frames,
    aes(x, y, group = step),
    linewidth = 1,
    colour = BLACK
  ) +
  geom_text(
    data = axis_end,
    aes(x, y, label = txt, group = step),
    parse = TRUE,
    size = 4.8
  ) +
  geom_text(
    data = emt_ticks_all,
    aes(x, y, label = txt),
    size = 3.6,
    colour = "black"
  ) +
  geom_text(
    data = before_labels_5,
    aes(x, y, label = txt),
    size = 3.6,
    colour = "#B00020"
  ) +
  geom_segment(
    data = arrow_before_5,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(6, "pt"), type = "closed"),
    linewidth = 0.6,
    colour = "#B00020"
  ) +
  geom_label(
    data = label_before_5,
    aes(x = lx, y = ly, label = lbl),
    angle = 270,
    label.size = 0.6,
    fill = "white",
    size = 3.8,
    label.r = unit(0.8, "pt"),
    label.padding = unit(2.2, "pt"),
    colour = "#B00020"
  ) +
  geom_text(
    data = after_labels_6,
    aes(x, y, label = txt),
    size = 3.6,
    colour = "#005BBB"
  ) +
  geom_segment(
    data = arrow_after_6,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(6, "pt"), type = "closed"),
    linewidth = 0.6,
    colour = "#005BBB"
  ) +
  geom_label(
    data = label_after_6,
    aes(x = lx, y = ly, label = lbl),
    angle = 90,
    label.size = 0.6,
    fill = "white",
    size = 3.8,
    label.r = unit(0.8, "pt"),
    label.padding = unit(2.2, "pt"),
    colour = "#005BBB"
  ) +
  facet_wrap(~step, nrow = 2) +
  coord_fixed(
    xlim = c(-0.25, x_off + 1.2),
    ylim = c(-0.22, y_scale + 0.15),
    expand = FALSE,
    clip = "off"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.title = element_blank(),
    strip.text = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(16, 12, 28, 12)
  )

# Display the plot
p