# Clear the workspace
rm(list = ls())

# Load required packages
library(ggplot2)
library(plotly)

# Detach ggtern if it is loaded, to avoid guide conflicts
if ("package:ggtern" %in% search()) {
  detach("package:ggtern", unload = TRUE, character.only = TRUE)
}

# Define shared plotting constants
BLACK <- "black"
WHITE <- "#FFFFFF"
sqrt3 <- sqrt(3)

# ============================================================
# Helper functions for simplex geometry
# ============================================================

# Convert barycentric coordinates (x1, x2, x3) to 2D Cartesian
# coordinates for an equilateral ternary triangle
bary_to_xy <- function(x1, x2, x3) {
  x <- x2 + 0.5 * x3
  y <- (sqrt(3) / 2) * x3
  data.frame(x, y)
}

# Build 2D ternary grid lines for specified mixture levels
build_grid_2d_levels <- function(levels = c(0.2, 0.4, 0.6, 0.8)) {
  mk <- function(p1, p2) {
    data.frame(x = p1$x, y = p1$y, xend = p2$x, yend = p2$y)
  }
  
  out <- list()
  for (t in levels) {
    out[[length(out) + 1]] <- mk(bary_to_xy(t, 0, 1 - t), bary_to_xy(t, 1 - t, 0)) # x1 = t
    out[[length(out) + 1]] <- mk(bary_to_xy(0, t, 1 - t), bary_to_xy(1 - t, t, 0)) # x2 = t
    out[[length(out) + 1]] <- mk(bary_to_xy(0, 1 - t, t), bary_to_xy(1 - t, 0, t)) # x3 = t
  }
  do.call(rbind, out)
}

# ============================================================
# Plot 1: One-dimensional simplex
# ============================================================

# Define the diagonal simplex from (0,1) to (1,0)
t_vals <- c(0.2, 0.4, 0.6, 0.8)

# Unit normal vector used to place tick marks and labels on either side
norm_vec <- c(1, 1) / sqrt(2)

tick_len  <- 0.02
label_off <- 0.06

# Tick marks centered on the simplex line
ticks_df <- do.call(rbind, lapply(t_vals, function(t) {
  Px <- t
  Py <- 1 - t
  data.frame(
    x    = Px - norm_vec[1] * tick_len,
    y    = Py - norm_vec[2] * tick_len,
    xend = Px + norm_vec[1] * tick_len,
    yend = Py + norm_vec[2] * tick_len
  )
}))

# Labels for r1 on the upper-right side of the simplex
labels_z1 <- do.call(rbind, lapply(t_vals, function(t) {
  Px <- t
  Py <- 1 - t
  data.frame(
    x   = Px + norm_vec[1] * label_off,
    y   = Py + norm_vec[2] * label_off,
    lab = sprintf('r[1] == %d*"%%"', round(100 * t))
  )
}))

# Labels for r2 on the lower-left side of the simplex
labels_z2 <- do.call(rbind, lapply(t_vals, function(t) {
  Px <- t
  Py <- 1 - t
  data.frame(
    x   = Px - norm_vec[1] * label_off,
    y   = Py - norm_vec[2] * label_off,
    lab = sprintf('r[2] == %d*"%%"', round(100 * (1 - t)))
  )
}))

# Create the 1D simplex plot
p1 <- ggplot() +
  geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0), linewidth = 1.2, colour = BLACK) +
  geom_point(aes(x = 0, y = 1), size = 3, colour = BLACK) +
  geom_point(aes(x = 1, y = 0), size = 3, colour = BLACK) +
  geom_segment(
    data = ticks_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.8,
    colour = BLACK
  ) +
  geom_text(data = labels_z1, aes(x, y, label = lab), parse = TRUE, hjust = 0, vjust = 0.5, size = 4) +
  geom_text(data = labels_z2, aes(x, y, label = lab), parse = TRUE, hjust = 1, vjust = 0.5, size = 4) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    x = expression(x[1]),
    y = expression(x[2]),
    title = "A)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title   = element_text(size = 24, face = "bold"),
    axis.title   = element_text(size = 20),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    panel.grid   = element_blank(),
    panel.border = element_rect(colour = BLACK, fill = NA, linewidth = 0.8)
  )

# ============================================================
# Plot 2: Two-dimensional simplex embedded in 3D space
# ============================================================

# Vertices of the triangular simplex
tri_verts <- data.frame(
  x = c(1, 0, 0),
  y = c(0, 1, 0),
  z = c(0, 0, 1)
)

# Triangle outline
tri_outline <- data.frame(
  x = c(1, 0, 0, 1),
  y = c(0, 1, 0, 0),
  z = c(0, 0, 1, 0)
)

# Construct ternary grid lines in 3D barycentric coordinates
bary_grid_levels_3d <- function(levels) {
  segs <- list()
  
  for (t in levels) {
    segs[[length(segs) + 1]] <- data.frame(
      x = c(t, t, NA), y = c(0, 1 - t, NA), z = c(1 - t, 0, NA)
    )
  }
  for (t in levels) {
    segs[[length(segs) + 1]] <- data.frame(
      x = c(0, 1 - t, NA), y = c(t, t, NA), z = c(1 - t, 0, NA)
    )
  }
  for (t in levels) {
    segs[[length(segs) + 1]] <- data.frame(
      x = c(0, 1 - t, NA), y = c(1 - t, 0, NA), z = c(t, t, NA)
    )
  }
  
  do.call(rbind, segs)
}

# Construct wireframe edges of the surrounding cube
cube_edges <- function() {
  e <- list()
  
  yz <- expand.grid(y = c(0, 1), z = c(0, 1))
  for (i in 1:nrow(yz)) {
    e[[length(e) + 1]] <- data.frame(x = c(0, 1, NA), y = rep(yz$y[i], 3), z = rep(yz$z[i], 3))
  }
  
  xz <- expand.grid(x = c(0, 1), z = c(0, 1))
  for (i in 1:nrow(xz)) {
    e[[length(e) + 1]] <- data.frame(x = rep(xz$x[i], 3), y = c(0, 1, NA), z = rep(xz$z[i], 3))
  }
  
  xy <- expand.grid(x = c(0, 1), y = c(0, 1))
  for (i in 1:nrow(xy)) {
    e[[length(e) + 1]] <- data.frame(x = rep(xy$x[i], 3), y = rep(xy$y[i], 3), z = c(0, 1, NA))
  }
  
  do.call(rbind, e)
}

grid3d <- bary_grid_levels_3d(c(0.2, 0.4, 0.6, 0.8))
cube_wire <- cube_edges()

# Create the 3D triangular simplex plot
p2 <- plot_ly() |>
  add_trace(
    type = "scatter3d", mode = "lines",
    x = cube_wire$x, y = cube_wire$y, z = cube_wire$z,
    line = list(width = 2, color = BLACK),
    hoverinfo = "none", showlegend = FALSE
  ) |>
  add_trace(
    type = "mesh3d",
    x = tri_verts$x, y = tri_verts$y, z = tri_verts$z,
    i = c(0), j = c(1), k = c(2),
    vertexcolor = rep(WHITE, nrow(tri_verts)),
    showscale = FALSE,
    lighting = list(ambient = 1, diffuse = 0, specular = 0),
    flatshading = TRUE,
    hoverinfo = "none"
  ) |>
  add_trace(
    type = "scatter3d", mode = "lines",
    x = grid3d$x, y = grid3d$y, z = grid3d$z,
    line = list(width = 2, color = BLACK),
    hoverinfo = "none", showlegend = FALSE
  ) |>
  add_trace(
    type = "scatter3d", mode = "lines",
    x = tri_outline$x, y = tri_outline$y, z = tri_outline$z,
    line = list(width = 4, color = BLACK),
    hoverinfo = "none", showlegend = FALSE
  ) |>
  layout(
    title = list(text = "B)", font = list(size = 24)),
    paper_bgcolor = "white",
    scene = list(
      xaxis = list(
        visible = TRUE, showgrid = FALSE, showticklabels = FALSE,
        ticks = "", zeroline = FALSE, showline = TRUE, linecolor = BLACK,
        range = c(0, 1), title = "x\u2081", titlefont = list(size = 20)
      ),
      yaxis = list(
        visible = TRUE, showgrid = FALSE, showticklabels = FALSE,
        ticks = "", zeroline = FALSE, showline = TRUE, linecolor = BLACK,
        range = c(0, 1), title = "x\u2082", titlefont = list(size = 20)
      ),
      zaxis = list(
        visible = TRUE, showgrid = FALSE, showticklabels = FALSE,
        ticks = "", zeroline = FALSE, showline = TRUE, linecolor = BLACK,
        range = c(0, 1), title = "x\u2083", titlefont = list(size = 20)
      ),
      bgcolor = "white",
      aspectmode = "data"
    ),
    margin = list(l = 0, r = 0, b = 0, t = 40)
  )

# ============================================================
# Plot 3: Two-dimensional simplex in equilateral ternary form
# ============================================================

# Triangle vertices in 2D
v1_2d <- bary_to_xy(1, 0, 0)
v2_2d <- bary_to_xy(0, 1, 0)
v3_2d <- bary_to_xy(0, 0, 1)
tri2d <- rbind(v1_2d, v2_2d, v3_2d, v1_2d)
center <- bary_to_xy(1 / 3, 1 / 3, 1 / 3)

# Compute outward unit normal for a given triangle edge
edge_outward_normal <- function(vA, vB) {
  e  <- c(vB$x - vA$x, vB$y - vA$y)
  n1 <- c(-e[2], e[1])
  n1 <- n1 / sqrt(sum(n1^2))
  n2 <- -n1
  
  mid <- c((vA$x + vB$x) / 2, (vA$y + vB$y) / 2)
  to_center <- c(center$x - mid[1], center$y - mid[2])
  
  if (sum(n1 * to_center) < 0) n1 else n2
}

# Build percentage labels outside a triangle edge
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

# Build axis labels outside the triangle
build_axis_label <- function(sym, vA, vB, offset = 0.10) {
  n <- edge_outward_normal(vA, vB)
  mid <- c((vA$x + vB$x) / 2, (vA$y + vB$y) / 2)
  
  data.frame(
    x = mid[1] + n[1] * offset,
    y = mid[2] + n[2] * offset,
    txt = sym
  )
}

grid2d <- build_grid_2d_levels(c(0.2, 0.4, 0.6, 0.8))
ticks <- c(0.2, 0.4, 0.6, 0.8)

# Assign each response axis to a triangle edge
edges <- list(
  z1 = list(A = v1_2d, B = v2_2d),
  z2 = list(A = v3_2d, B = v1_2d),
  z3 = list(A = v2_2d, B = v3_2d)
)

labs_z1 <- build_edge_tick_labels(edges$z1$A, edges$z1$B, ticks, offset = 0.06)
labs_z2 <- build_edge_tick_labels(edges$z2$A, edges$z2$B, ticks, offset = 0.06)
labs_z3 <- build_edge_tick_labels(edges$z3$A, edges$z3$B, ticks, offset = 0.06)

axis_labs <- rbind(
  build_axis_label("r[1]", edges$z1$A, edges$z1$B, offset = 0.10),
  build_axis_label("r[2]", edges$z2$A, edges$z2$B, offset = 0.10),
  build_axis_label("r[3]", edges$z3$A, edges$z3$B, offset = 0.10)
)

# Create the equilateral ternary plot
p3 <- ggplot() +
  geom_path(data = tri2d, aes(x, y), linewidth = 1, colour = BLACK) +
  geom_segment(
    data = grid2d,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5,
    alpha = 0.8,
    colour = BLACK
  ) +
  geom_text(data = labs_z1, aes(x, y, label = txt), size = 4, colour = BLACK) +
  geom_text(data = labs_z2, aes(x, y, label = txt), size = 4, colour = BLACK) +
  geom_text(data = labs_z3, aes(x, y, label = txt), size = 4, colour = BLACK) +
  geom_text(data = axis_labs, aes(x, y, label = txt), parse = TRUE, size = 7, colour = BLACK) +
  coord_fixed(
    xlim = c(-0.15, 1.15),
    ylim = c(-0.15, sqrt3 / 2 + 0.15),
    expand = FALSE
  ) +
  labs(title = "C)") +
  theme_bw(base_size = 16) +
  theme(
    plot.title   = element_text(size = 24, face = "bold"),
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    panel.grid   = element_blank(),
    panel.border = element_blank()
  )

# ============================================================
# Plot 4: Three-dimensional simplex as a regular tetrahedron
# ============================================================

# Define tetrahedron vertices
v1 <- c(0, 0, 0)
v2 <- c(1, 0, 0)
v3 <- c(0.5, sqrt(3) / 2, 0)
v4 <- c(0.5, sqrt(3) / 6, sqrt(6) / 3)

vx <- c(v1[1], v2[1], v3[1], v4[1])
vy <- c(v1[2], v2[2], v3[2], v4[2])
vz <- c(v1[3], v2[3], v3[3], v4[3])

# Face definitions for plotly mesh3d
i <- c(0, 0, 0, 1)
j <- c(1, 1, 2, 2)
k <- c(2, 3, 3, 3)

# Return the coordinates of a tetrahedron vertex
v_at <- function(k) c(vx[k], vy[k], vz[k])

# Construct triangular cross-sections corresponding to x_k = t
plane_triangle_edges <- function(k_idx, t) {
  others <- setdiff(1:4, k_idx)
  vk <- v_at(k_idx)
  
  p1 <- t * vk + (1 - t) * v_at(others[1])
  p2 <- t * vk + (1 - t) * v_at(others[2])
  p3 <- t * vk + (1 - t) * v_at(others[3])
  
  data.frame(
    x = c(p1[1], p2[1], NA_real_, p2[1], p3[1], NA_real_, p3[1], p1[1], NA_real_),
    y = c(p1[2], p2[2], NA_real_, p2[2], p3[2], NA_real_, p3[2], p1[2], NA_real_),
    z = c(p1[3], p2[3], NA_real_, p2[3], p3[3], NA_real_, p3[3], p1[3], NA_real_)
  )
}

# Build interior grid from triangular cross-sections
levels <- c(0.2, 0.4, 0.6, 0.8)
internal_wire_list <- list()
idx <- 0L

for (kk in 1:4) {
  for (tt in levels) {
    idx <- idx + 1L
    internal_wire_list[[idx]] <- plane_triangle_edges(kk, tt)
  }
}

internal_wire <- do.call(rbind, internal_wire_list)

# Build outer tetrahedron wireframe
edge_pairs <- rbind(
  c(1, 2), c(1, 3), c(1, 4),
  c(2, 3), c(2, 4),
  c(3, 4)
)

outer_wire_list <- lapply(seq_len(nrow(edge_pairs)), function(r) {
  a <- edge_pairs[r, 1]
  b <- edge_pairs[r, 2]
  data.frame(
    x = c(vx[a], vx[b], NA_real_),
    y = c(vy[a], vy[b], NA_real_),
    z = c(vz[a], vz[b], NA_real_)
  )
})

outer_wire <- do.call(rbind, outer_wire_list)

# Create the tetrahedral simplex plot
p4 <- plot_ly() |>
  add_trace(
    type = "mesh3d",
    x = vx, y = vy, z = vz,
    i = i, j = j, k = k,
    vertexcolor = rep(WHITE, 4),
    opacity = 0.32,
    showscale = FALSE,
    lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
    flatshading = TRUE,
    hoverinfo = "none"
  ) |>
  add_trace(
    type = "scatter3d", mode = "lines",
    x = outer_wire$x, y = outer_wire$y, z = outer_wire$z,
    line = list(width = 4, color = BLACK),
    hoverinfo = "none", showlegend = FALSE
  ) |>
  add_trace(
    type = "scatter3d", mode = "lines",
    x = internal_wire$x, y = internal_wire$y, z = internal_wire$z,
    line = list(width = 1.2, color = BLACK),
    opacity = 0.45,
    hoverinfo = "none", showlegend = FALSE
  ) |>
  layout(
    title = list(text = "D)", font = list(size = 24)),
    paper_bgcolor = "white",
    scene = list(
      xaxis = list(visible = FALSE),
      yaxis = list(visible = FALSE),
      zaxis = list(visible = FALSE),
      bgcolor = "white",
      aspectmode = "data"
    ),
    margin = list(l = 0, r = 0, b = 0, t = 40)
  )

# ============================================================
# Display plots
# ============================================================

p1
p2
p3
p4