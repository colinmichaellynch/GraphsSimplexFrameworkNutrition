# Clear the workspace
rm(list = ls())

# Define the nutrient composition of the available foods
foods <- rbind(
  c(0.75, 0.15, 0.10),
  c(0.20, 0.85, 0.05),
  c(0.15, 0.30, 0.55),
  c(0.35, 0.50, 0.15),
  c(0.40, 0.15, 0.45)
)
colnames(foods) <- c("z1", "z2", "z3")

# Load required package
set.seed(42)
library(ggplot2)

# Define shared plotting constants
BLACK <- "black"
sqrt3 <- sqrt(3)

# ============================================================
# Helper functions for equilateral mixture triangle geometry
# ============================================================

# Convert barycentric coordinates (z1, z2, z3) to 2D Cartesian
# coordinates for an equilateral triangle
bary_to_xy <- function(z1, z2, z3) {
  x <- z2 + 0.5 * z3
  y <- (sqrt(3) / 2) * z3
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

# ============================================================
# Define colors and simplex grid
# ============================================================

# Colors for the internal grid lines corresponding to each nutrient axis
axis_cols_inside <- c(
  z1 = "#00A000",
  z2 = "#D00000",
  z3 = "#0055FF"
)

# Colors for the axis labels and tick labels
axis_cols_labels <- c(
  z1 = "#0055FF",
  z2 = "#00A000",
  z3 = "#D00000"
)

# Mixture levels used for internal grid lines
levs <- c(0.2, 0.4, 0.6, 0.8)

# Helper to build a line segment between two projected points
mkseg <- function(p1, p2) {
  data.frame(x = p1$x, y = p1$y, xend = p2$x, yend = p2$y)
}

# Internal grid lines for each simplex coordinate
grid2d_z1 <- do.call(rbind, lapply(levs, function(t) {
  mkseg(bary_to_xy(t, 0, 1 - t), bary_to_xy(t, 1 - t, 0))
}))

grid2d_z2 <- do.call(rbind, lapply(levs, function(t) {
  mkseg(bary_to_xy(0, t, 1 - t), bary_to_xy(1 - t, t, 0))
}))

grid2d_z3 <- do.call(rbind, lapply(levs, function(t) {
  mkseg(bary_to_xy(0, 1 - t, t), bary_to_xy(1 - t, 0, t))
}))

# Define the triangle edges associated with each nutrient axis
edges <- list(
  z1 = list(A = v2, B = v3),
  z2 = list(A = v3, B = v1),
  z3 = list(A = v1, B = v2)
)

# ============================================================
# Functions for external tick labels and axis labels
# ============================================================

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

# Build a label centered outside a triangle edge
build_axis_label <- function(label, vA, vB, offset = 0.10) {
  n <- edge_outward_normal(vA, vB)
  mid <- c((vA$x + vB$x) / 2, (vA$y + vB$y) / 2)
  
  data.frame(
    x = mid[1] + n[1] * offset,
    y = mid[2] + n[2] * offset,
    txt = label
  )
}

ticks <- c(0.2, 0.4, 0.6, 0.8)

# Tick labels along each edge
labs_z1 <- build_edge_tick_labels(edges$z1$A, edges$z1$B, ticks, offset = 0.06)
labs_z2 <- build_edge_tick_labels(edges$z2$A, edges$z2$B, ticks, offset = 0.06)
labs_z3 <- build_edge_tick_labels(edges$z3$A, edges$z3$B, ticks, offset = 0.06)

# Nutrient axis labels
axis_lab_z1 <- build_axis_label("Fat %",     edges$z1$A, edges$z1$B, offset = 0.10)
axis_lab_z2 <- build_axis_label("Protein %", edges$z2$A, edges$z2$B, offset = 0.10)
axis_lab_z3 <- build_axis_label("Carb %",    edges$z3$A, edges$z3$B, offset = 0.10)

# ============================================================
# Project foods into 2D and compute summary quantities
# ============================================================

# Assign equal weights to all foods
p <- rep(1 / nrow(foods), nrow(foods))

# Project food compositions into Cartesian coordinates
foods_xy <- bary_to_xy(foods[, "z1"], foods[, "z2"], foods[, "z3"])
foods_df <- cbind(as.data.frame(foods), foods_xy)

# Compute the center of mass of the foods
zbar <- colSums(foods * p)
zbar_xy <- bary_to_xy(zbar[1], zbar[2], zbar[3])

# Construct a target point as a convex combination of foods, then nudge it
# slightly within the simplex
w <- runif(nrow(foods))
w <- w / sum(w)
target <- as.numeric(t(w) %*% foods)

dx <- 0.08
dy <- 0.045
dz3 <- (2 / sqrt3) * dy
dz2 <- dx - 0.5 * dz3
dz1 <- -(dz2 + dz3)

target <- target + c(dz1, dz2, dz3)
target <- pmax(1e-6, target)
target <- target / sum(target)
target_xy <- bary_to_xy(target[1], target[2], target[3])

# ============================================================
# Compute a 95% confidence ellipse around the food center
# ============================================================

# Approximate the covariance of the sample mean in projected 2D space
n_ind <- 30
mx <- sum(p * foods_df$x)
my <- sum(p * foods_df$y)

cx  <- sum(p * (foods_df$x - mx)^2)
cy  <- sum(p * (foods_df$y - my)^2)
cxy <- sum(p * (foods_df$x - mx) * (foods_df$y - my))

Sigma_one  <- matrix(c(cx, cxy, cxy, cy), 2, 2, byrow = TRUE)
Sigma_mean <- Sigma_one / n_ind

level <- 0.95
chi2_95 <- qchisq(level, df = 2)

# Construct ellipse coordinates from the covariance matrix
Ssym <- (Sigma_mean + t(Sigma_mean)) / 2
eg <- eigen(Ssym, symmetric = TRUE)
L <- eg$vectors %*% diag(sqrt(pmax(eg$values, 0)))

theta <- seq(0, 2 * pi, length.out = 360)
unitC <- rbind(cos(theta), sin(theta))
E <- L %*% (sqrt(chi2_95) * unitC)

ell_df <- data.frame(
  x = mx + E[1, ],
  y = my + E[2, ]
)

# ============================================================
# Compute the convex hull of the foods
# ============================================================

h_idx <- chull(foods_df$x, foods_df$y)
h_idx <- c(h_idx, h_idx[1])
hull_df <- foods_df[h_idx, c("x", "y")]

# ============================================================
# Functions for deterministic feeding trajectories
# ============================================================

# Move the current intake state partway toward a selected food
ema_step <- function(cur, food, alpha = 0.20) {
  nxt <- (1 - alpha) * cur + alpha * food
  s <- sum(nxt)
  if (s != 0) nxt <- nxt / s
  nxt
}

# Greedy policy that chooses the food minimizing squared distance to a goal
greedy_to_goal_policy <- function(goal, alpha = 0.20) {
  function(cur, foods) {
    dists <- apply(foods, 1, function(fj) {
      cand <- (1 - alpha) * cur + alpha * fj
      sum((cand - goal)^2)
    })
    which.min(dists)
  }
}

# Simulate a feeding trajectory under a specified policy
simulate_traj <- function(start_bary, steps, policy_fun, alpha = 0.20) {
  cur <- start_bary
  out <- matrix(NA_real_, nrow = steps + 1, ncol = 3)
  out[1, ] <- cur
  
  for (t in 1:steps) {
    j <- policy_fun(cur, foods)
    cur <- ema_step(cur, foods[j, ], alpha)
    out[t + 1, ] <- cur
  }
  
  out <- as.data.frame(out)
  names(out) <- c("z1", "z2", "z3")
  xy <- bary_to_xy(out$z1, out$z2, out$z3)
  
  cbind(xy, t = seq_len(nrow(out)))
}

# ============================================================
# Simulate two example trajectories
# ============================================================

# Define starting points and number of steps
start_A <- c(0.94, 0.04, 0.02)
start_B <- c(0.04, 0.02, 0.94)
steps <- 12

# Simulate one trajectory toward the food center of mass and one toward the target
trajA <- simulate_traj(start_A, steps, greedy_to_goal_policy(zbar, alpha = 0.20), alpha = 0.20)
trajB <- simulate_traj(start_B, steps, greedy_to_goal_policy(target, alpha = 0.20), alpha = 0.20)

traj_col <- "#FF7F0E"

# ============================================================
# Build the equilateral mixture triangle plot
# ============================================================

p_emt <- ggplot() +
  geom_path(data = tri2d, aes(x, y), linewidth = 1.2, colour = BLACK) +
  geom_segment(
    data = grid2d_z1,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.9,
    alpha = 1,
    colour = axis_cols_inside["z1"]
  ) +
  geom_segment(
    data = grid2d_z2,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.9,
    alpha = 1,
    colour = axis_cols_inside["z2"]
  ) +
  geom_segment(
    data = grid2d_z3,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.9,
    alpha = 1,
    colour = axis_cols_inside["z3"]
  ) +
  geom_text(data = labs_z1, aes(x, y, label = txt), size = 4.2, colour = axis_cols_labels["z1"]) +
  geom_text(data = labs_z2, aes(x, y, label = txt), size = 4.2, colour = axis_cols_labels["z2"]) +
  geom_text(data = labs_z3, aes(x, y, label = txt), size = 4.2, colour = axis_cols_labels["z3"]) +
  geom_text(
    data = axis_lab_z1,
    aes(x, y, label = txt),
    size = 5.6,
    fontface = "bold",
    colour = axis_cols_labels["z1"]
  ) +
  geom_text(
    data = axis_lab_z2,
    aes(x, y, label = txt),
    size = 5.6,
    fontface = "bold",
    colour = axis_cols_labels["z2"]
  ) +
  geom_text(
    data = axis_lab_z3,
    aes(x, y, label = txt),
    size = 5.6,
    fontface = "bold",
    colour = axis_cols_labels["z3"]
  ) +
  geom_polygon(data = hull_df, aes(x, y), fill = "grey70", alpha = 0.30, colour = NA) +
  geom_path(data = hull_df, aes(x, y), linewidth = 1.1, colour = "grey35") +
  geom_path(data = ell_df, aes(x, y), linewidth = 1.05, linetype = "dotted", colour = "steelblue4") +
  geom_path(data = trajA, aes(x, y), linewidth = 1.2, colour = traj_col, linetype = "solid",  alpha = 0.98) +
  geom_path(data = trajB, aes(x, y), linewidth = 1.2, colour = traj_col, linetype = "dashed", alpha = 0.98) +
  geom_point(data = foods_df, aes(x, y), shape = 16, size = 3.4, colour = BLACK) +
  geom_point(aes(x = zbar_xy$x,   y = zbar_xy$y),   shape = 16, size = 3.8, colour = "steelblue4") +
  geom_point(aes(x = target_xy$x, y = target_xy$y), shape = 16, size = 4.0, colour = "goldenrod2") +
  coord_fixed(
    xlim = c(-0.15, 1.15),
    ylim = c(-0.15, sqrt3 / 2 + 0.15),
    expand = FALSE
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    panel.grid   = element_blank(),
    panel.border = element_blank()
  )

# Display the plot
p_emt
