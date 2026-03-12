# Equilateral Mixture Triangles in Nutritional Geometry

This repository contains the R scripts used to generate the figures for our paper on the use of **equilateral mixture triangles (EMTs)** in nutritional geometry.

## Overview

Animal nutritionists seek to understand how organisms regulate the intake and balance of multiple nutrients, but the design and analysis of these experiments are often shaped by how nutrient spaces are represented. In studies involving more than two nutrients, researchers have commonly used **right-angled mixture triangles (RMTs)**. Although convenient, RMTs can distort proportional relationships among nutrients and introduce visual artifacts that complicate interpretation.

In this paper, we reintroduce the **equilateral mixture triangle (EMT)** as a complementary and biologically meaningful representation of nutrient space. Because the EMT preserves mixture geometry, it provides a direct bridge between the **geometric framework for nutrition (GFN)** and the broader theory of **simplex-based mixture experiments** used in engineering and statistics.

Using this framework, we show how EMTs can be used to:

- represent multidimensional nutrient mixtures without the geometric distortions of RMTs,
- visualize how feeding trajectories are constrained by the **convex hull** of available foods,
- extend GFN ideas to higher-dimensional nutrient systems,
- distinguish **random choice** from **defense of an intake target**, and
- connect no-choice nutritional experiments to **response surface modeling** and **mixture design optimization**.

## Summary of the paper

The paper argues that EMTs provide a more geometrically faithful way to represent nutrient mixtures than RMTs when three or more nutrients are involved. This matters because the choice of coordinate system can shape visual interpretation, theoretical reasoning, and experimental design.

More specifically, the paper shows that EMTs can unify several parts of nutritional geometry that are often treated separately:

- **Choice experiments:** feeding trajectories are constrained by the convex hull of available foods, and this geometry can be used to test whether observed trajectories are consistent with random choice or with defense of an intake target.
- **Higher-dimensional nutrient systems:** the same geometric reasoning extends naturally beyond three nutrients.
- **No-choice experiments:** EMTs connect directly to simplex mixture designs, allowing response-surface methods to be used for estimating performance landscapes and identifying optimal nutrient ratios.

Together, these results position the EMT as a bridge between biological theory and the statistical design of mixture experiments.

## Repository contents

This repository contains **9 R scripts**, each of which generates one or more figure panels from the paper.

Scripts are named according to the figure and panel(s) they produce. For example:

- `Figure1ABCD.R`
- `Figure2AB.R`
- `Figure3C.R`

In general, the naming convention is:

```text
Figure[number][panel letters].R

So a script named `Figure1ABCD.R` produces panels A–D of Figure 1.

## Requirements

The scripts use base R plus several contributed packages. Depending on which figure you run, you may need some or all of the following:

- `ggplot2`
- `plotly`
- `patchwork`
- `dplyr`
- `fields`
- `AlgDesign`
- `dtw`

You can install them with:

```r id="install_pkgs_01"
install.packages(c(
  "ggplot2",
  "plotly",
  "patchwork",
  "dplyr",
  "fields",
  "AlgDesign",
  "dtw"
))

## How to run the scripts

Each script is designed to be run independently.

1. Clone or download this repository.
2. Open R or RStudio.
3. Set your working directory to the repository folder if needed.
4. Run the script corresponding to the figure you want to reproduce.

For example:

```r id="run_example_01"
source("Figure1ABCD.R")

Most scripts will display the figure directly in the plotting window. Some scripts also print numerical summaries or diagnostic output to the console.

## Notes on reproducibility

- Scripts are intended to be self-contained, so helper functions and data objects are typically defined within each file.
- Some figures rely on simulation or stochastic procedures. Where relevant, random seeds are set within the scripts so that outputs are reproducible.
- A small number of scripts print summary statistics in addition to producing figures.
- Because the scripts were written on a figure-by-figure basis, some helper functions appear in multiple files. This is intentional and helps keep each script portable and easy to run on its own.

## Purpose of the repository

This repository is intended to make the visual and computational components of the paper transparent and reproducible. Beyond reproducing the figures themselves, the scripts illustrate how EMTs can be used for:

- geometric reasoning about nutrient mixtures,
- comparison of alternative feeding rules,
- analysis of intake trajectories,
- simplex-lattice and simplex-centroid designs,
- constrained mixture optimization, and
- response surface visualization in nutritional geometry.

## Citation

If you use these scripts or build on this work, please cite the associated paper.

## Contact

For questions, suggestions, or reproducibility issues, please open an issue in this repository.
