# ==============================================================================
# Phylogenetic Tree Visualization of Haemaphysalis Mitochondrial Genomes
# ==============================================================================
#
# Description:
# This script demonstrates how to visualize a phylogenetic tree generated from
# mitochondrial genome sequences of Haemaphysalis ticks. 
#
# IMPORTANT:
# Before running the script, change the INPUT PATH, TREE FILE, and
# METADATA FILE below according to your own dataset.
#
# Required metadata columns:
#   Accession  : Accession number / sequence identifier matching tree tips
#   Species    : Species name
#   Genus      : Genus name
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. LOAD REQUIRED PACKAGES
# ------------------------------------------------------------------------------

library(ape)          # Phylogenetic tree manipulation
library(ggtree)       # Phylogenetic tree visualization
library(ggplot2)      # Plotting
library(dplyr)        # Data manipulation
library(readxl)       # Reading Excel metadata files
library(phytools)     # Phylogenetic utilities
library(svglite)      # SVG output
library(grid)         # Plot layout utilities


# ------------------------------------------------------------------------------
# 2. SET INPUT AND OUTPUT PATHS
# ------------------------------------------------------------------------------

# >>> CHANGE THIS PATH <<<
# Directory containing the tree and metadata files.
input_dir <- "PATH/TO/YOUR/INPUT/FILES"

# >>> CHANGE THESE FILE NAMES <<<
tree_file <- "your_phylogenetic_tree.contree"
metadata_file <- "your_metadata.xlsx"

# Output directory for generated figures
output_dir <- file.path(input_dir, "Rplots")

# Create output directory if it does not already exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 3. READ INPUT DATA
# ------------------------------------------------------------------------------

# Read phylogenetic tree in Newick format
tree <- read.tree(file.path(input_dir, tree_file))

# Read metadata table
metadata <- read_excel(file.path(input_dir, metadata_file))

# ------------------------------------------------------------------------------
# 4. CHECK AND PREPARE METADATA
# ------------------------------------------------------------------------------

# Required metadata columns
required_columns <- c("Accession", "Species", "Genus")

missing_columns <- setdiff(required_columns, colnames(metadata))

if (length(missing_columns) > 0) {
  stop(
    paste(
      "The following required metadata columns are missing:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# Convert relevant columns to character
metadata <- metadata %>%
  mutate(
    Accession = as.character(Accession),
    Species   = as.character(Species),
    Genus     = as.character(Genus)
  )

# ------------------------------------------------------------------------------
# 5. IDENTIFY STUDY TICK SAMPLES
# ------------------------------------------------------------------------------

# Tick samples are identified when "Tick" occurs in the accession/sample ID.
# Modify this condition if a different naming convention is used in your dataset.

metadata <- metadata %>%
  mutate(
    tick_mark = grepl("Tick", Accession, ignore.case = TRUE)
  )

# ------------------------------------------------------------------------------
# 6. SELECT THE TARGET GENUS
# ------------------------------------------------------------------------------

# >>> CHANGE THIS IF REQUIRED <<<
target_genus <- "Haemaphysalis"

metadata_filtered <- metadata %>%
  filter(Genus == target_genus)

# Match metadata accessions with tree tip labels
keep_tips <- intersect(
  tree$tip.label,
  metadata_filtered$Accession
)

if (length(keep_tips) == 0) {
  stop(
    paste0(
      "No matching tree tips were found for genus: ",
      target_genus,
      ". Check the Genus column and tree tip labels."
    )
  )
}

# Retain only matching sequences
tree_filtered <- keep.tip(tree, keep_tips)

# ------------------------------------------------------------------------------
# 7. MIDPOINT ROOT THE TREE
# ------------------------------------------------------------------------------

rooted_tree <- midpoint.root(tree_filtered)

# ------------------------------------------------------------------------------
# 8. PREPARE TREE DATA FOR ANNOTATION
# ------------------------------------------------------------------------------

base_tree <- ggtree(rooted_tree)

tree_data <- base_tree$data

# Combine tree information with metadata
annotated_data <- full_join(
  tree_data,
  metadata_filtered,
  by = c("label" = "Accession")
)

# Tip data
tip_data <- annotated_data %>%
  filter(isTip)

# Tick samples to be highlighted
tick_data <- annotated_data %>%
  filter(isTip, tick_mark)

# ------------------------------------------------------------------------------
# 9. DEFINE VISUALIZATION PARAMETERS
# ------------------------------------------------------------------------------

# Assign a single colour to the target genus.
# This can be modified if multiple genera are included.
genera <- sort(unique(na.omit(metadata_filtered$Genus)))

genus_colors <- setNames(
  rep("#00BA38", length(genera)),
  genera
)

# Tree dimensions
n_tips <- length(rooted_tree$tip.label)

tree_height <- max(
  tree_data$x,
  na.rm = TRUE
)

# Position of annotation columns
accession_x <- tree_height * 1.20
species_x   <- tree_height * 1.38

line_end_x <- accession_x - tree_height * 0.02

# Automatically adjust figure height according to number of tips
plot_height <- max(8.5, n_tips * 0.11)

plot_width <- 14

# ------------------------------------------------------------------------------
# 10. GENERATE PHYLOGENETIC TREE
# ------------------------------------------------------------------------------

p <- ggtree(
  rooted_tree,
  size = 0.35,
  color = "black"
) %<+% annotated_data +

  # --------------------------------------------------------------------------
  # Connector lines between tree tips and accession labels
  # --------------------------------------------------------------------------
  geom_segment(
    data = tip_data,
    aes(
      x = x,
      xend = line_end_x,
      y = y,
      yend = y,
      color = Genus
    ),
    linewidth = 0.20,
    inherit.aes = FALSE,
    na.rm = TRUE,
    show.legend = FALSE
  ) +

  # --------------------------------------------------------------------------
  # Highlight study tick samples using red dots
  # --------------------------------------------------------------------------
  geom_point(
    data = tick_data,
    aes(
      x = line_end_x,
      y = y
    ),
    color = "red",
    fill = "red",
    shape = 16,
    size = 1.5,
    inherit.aes = FALSE,
    na.rm = TRUE
  ) +

  # --------------------------------------------------------------------------
  # Accession / sample ID labels
  # --------------------------------------------------------------------------
  geom_text(
    data = tip_data,
    aes(
      x = accession_x,
      y = y,
      label = label,
      color = Genus
    ),
    hjust = 0,
    size = 2.6,
    fontface = "italic",
    inherit.aes = FALSE,
    na.rm = TRUE,
    show.legend = FALSE
  ) +

  # --------------------------------------------------------------------------
  # Species labels
  # --------------------------------------------------------------------------
  geom_text(
    data = tip_data,
    aes(
      x = species_x,
      y = y,
      label = Species,
      color = Genus
    ),
    hjust = 0,
    size = 2.6,
    fontface = "italic",
    inherit.aes = FALSE,
    na.rm = TRUE,
    show.legend = FALSE
  ) +

  # --------------------------------------------------------------------------
  # Invisible points used to generate the genus legend
  # --------------------------------------------------------------------------
  geom_point(
    data = tip_data,
    aes(
      x = x,
      y = y,
      color = Genus
    ),
    inherit.aes = FALSE,
    size = 0.01,
    alpha = 0,
    show.legend = TRUE
  ) +

  # --------------------------------------------------------------------------
  # Display bootstrap support values >60
  # --------------------------------------------------------------------------
  geom_text2(
    aes(
      subset =
        !isTip &
        !is.na(suppressWarnings(as.numeric(label))) &
        suppressWarnings(as.numeric(label)) > 60,
      label = label
    ),
    hjust = 0.5,
    vjust = -0.2,
    size = 2.5,
    color = "blue"
  ) +

  # --------------------------------------------------------------------------
  # Phylogenetic scale bar
  # --------------------------------------------------------------------------
  geom_treescale(
    x = 0,
    y = -1,
    width = 0.05,
    fontsize = 2.5,
    color = "black",
    offset = 0.6
  ) +

  # --------------------------------------------------------------------------
  # Genus colour legend
  # --------------------------------------------------------------------------
  scale_color_manual(
    values = genus_colors,
    na.value = "grey40",
    name = "Genus",
    breaks = genera
  ) +

  # Extend plotting area to accommodate labels
  xlim(
    0,
    species_x + tree_height * 0.42
  ) +

  # Allow labels to extend outside the tree panel
  coord_cartesian(
    clip = "off"
  ) +

  # Tree theme
  theme_tree2() +

  # --------------------------------------------------------------------------
  # Publication-oriented figure formatting
  # --------------------------------------------------------------------------
  theme(
    legend.position = c(0.93, 0.52),
    legend.justification = c(0, 0.5),

    legend.title = element_text(
      size = 8,
      face = "bold"
    ),

    legend.text = element_text(
      size = 6.5
    ),

    legend.key.size = unit(
      0.45,
      "lines"
    ),

    legend.key = element_rect(
      fill = "white",
      color = NA
    ),

    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.6
    ),

    legend.box.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.6
    ),

    legend.margin = margin(
      4, 5, 4, 5
    ),

    legend.box.margin = margin(
      2, 2, 2, 2
    ),

    legend.spacing.y = unit(
      0.08,
      "cm"
    ),

    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),

    # Additional right margin for accession and species labels
    plot.margin = margin(
      6, 170, 6, 6
    )
  ) +

  # Square legend symbols
  guides(
    color = guide_legend(
      override.aes = list(
        shape = 15,
        size = 3.5,
        alpha = 1,
        linetype = 0
      ),
      ncol = 1
    )
  )

# ------------------------------------------------------------------------------
# 11. EXPORT PUBLICATION-QUALITY FIGURES
# ------------------------------------------------------------------------------

# PDF
ggsave(
  filename = file.path(
    output_dir,
    paste0("phylogenetic_tree_", target_genus, ".pdf")
  ),
  plot = p,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 600,
  limitsize = FALSE
)

# SVG
ggsave(
  filename = file.path(
    output_dir,
    paste0("phylogenetic_tree_", target_genus, ".svg")
  ),
  plot = p,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 600,
  limitsize = FALSE
)

# PNG
ggsave(
  filename = file.path(
    output_dir,
    paste0("phylogenetic_tree_", target_genus, ".png")
  ),
  plot = p,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 600,
  limitsize = FALSE,
  bg = "white"
)

# ------------------------------------------------------------------------------
# 12. COMPLETION MESSAGE
# ------------------------------------------------------------------------------

message(
  paste0(
    "\nPhylogenetic tree generated successfully.\n",
    "Target genus: ", target_genus, "\n",
    "Number of tips: ", n_tips, "\n",
    "Output directory: ", output_dir, "\n",
    "Formats: PDF, SVG, PNG\n"
  )
)

# ==============================================================================
