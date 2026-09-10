# Tick Mitochondrial Genome Analysis

This folder contains the command-line workflows and R scripts used for
mitochondrial genome analysis of *Haemaphysalis* ticks.

### Contents

* **Command-line workflow** for reference-based mitochondrial genome assembly,
  consensus generation, and low/no-coverage masking from Oxford Nanopore reads.
* **R scripts** for phylogenetic analysis, tree visualization, and comparative
  analysis of mitochondrial genomes.

### Command-line Workflow

Tools used include `Porechop`, `minimap2`, `samtools`, `bcftools`, and `bedtools`.

### Phylogenetic Analysis

R scripts require a Newick-format phylogenetic tree and an Excel metadata file
containing `Accession`, `Species`, and `Genus`.

### Output

The workflows generate final mitochondrial genome sequences and
publication-quality phylogenetic trees in **PDF, SVG, and PNG** formats.
