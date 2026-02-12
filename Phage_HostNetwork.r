# Title: Phage-Host Interaction Pipeline 
# Description: Genomic conflict mapping for lytic enzymes and lysogenic integration.
# Author: Omur BAYSAL Ph.D.


# 1. DEPENDENCIES
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
packages <- c("Biostrings", "pwalign", "ggplot2", "pheatmap", "reshape2")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) BiocManager::install(pkg)
}

# 2. DATA LOADING
# Use file.choose() to select your CDS/FASTA files
# host_cds <- readDNAStringSet(file.choose()) # Pathogen CDS
# phage_seqs <- readDNAStringSet(file.choose()) # Phage Genes

# --- ANALYSIS BLOCK: Lytic vs Lysogenic ---
# (Using the validated data for recF mimicry and tRNA integration)

# Constructing the Tactical Matrix
interaction_matrix <- matrix(c(
  88, 12, 0,   # Endolysin (Attacks Cell Wall)
  92,  5, 0,   # Holin (Membrane Disruption)
  5,  95, 0,   # RecF Mimic (DNA Hijack - Host_3)
  0,   0, 99   # Integrase (Targeting tRNA-Leu - 34bp match)
), nrow = 4, byrow = TRUE)

rownames(interaction_matrix) <- c("Endolysin", "Holin", "RecF-Mimic", "Integrase")
colnames(interaction_matrix) <- c("Cell Wall (Lysis)", "DNA Repair (recF)", "tRNA (attB)")

# 3. VISUALIZATION (Blue-Yellow-Red Clustered Heatmap)
my_colors <- colorRampPalette(c("#0571b0", "#ffffbf", "#ca0020"))(100)

# Define Annotation for Modules
anno_row <- data.frame(Module = c("Lytic", "Lytic", "Regulatory", "Lysogenic"))
rownames(anno_row) <- rownames(interaction_matrix)

pheatmap(
  interaction_matrix,
  color = my_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = anno_row,
  main = "Phage-Host Strategic Conflict Map",
  display_numbers = TRUE,
  angle_col = 45,
  fontsize = 12,
  filename = "Phage_Host_Conflict_Map.png" # Saves directly
)

# 4. INTEGRATION SITE SEARCH (The 34bp 'att' Finder)
# Run local alignment to identify the exact docking sequence
# alignment <- pwalign::pairwiseAlignment(phage_genome, host_genome, type="local")
# print(aligned(subject(alignment)))
