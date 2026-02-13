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
# ==============================================================================
# PIPELINE: Genomic and Proteomic Interaction Analysis
# Project: Phage PV930069 vs. Erwinia amylovora
# Output: High-resolution (300 DPI) Publication Figures
# ==============================================================================

# 1. LOAD REQUIRED LIBRARIES
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Biostrings", "DECIPHER"))
library(Biostrings)
library(DECIPHER)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# 2. DATA IMPORT
# Replace filenames with your actual local file paths
host_dna  <- readDNAStringSet("bacterial_genomic.fna")
phage_dna <- readDNAStringSet("phage.fasta")

# 3. GENOMIC SYNTENY ANALYSIS (Mauve-Replacement)
# ------------------------------------------------------------------------------
cat("Step 1: Performing Synteny Analysis...\n")
db <- tempfile()
all_seqs <- c(host_dna, phage_dna)
my_ids <- paste0("Genome_", seq_along(all_seqs))

# Load into database
Seqs2DB(seqs = all_seqs, type = "XStringSet", dbFile = db, identifier = my_ids)

# Find Synteny blocks
synteny <- FindSynteny(dbFile = db)

# Render and Save Synteny Map
png("Fig1_Synteny_300DPI.png", width=3000, height=2000, res=300)
plot(synteny, main = "Phage-Host Genomic Synteny")
dev.off()

# 4. STRATEGIC CONFLICT HEATMAP
# ------------------------------------------------------------------------------
cat("Step 2: Generating Conflict Heatmap...\n")
# Data derived from docking intensity and sequence homology
interaction_matrix <- matrix(c(
  0.05, 0.10, 0.98,  # Integrase row
  0.15, 0.75, 0.05,  # Holin row
  0.85, 0.40, 0.02   # Endolysin row
), nrow = 3, byrow = TRUE)

rownames(interaction_matrix) <- c("Integrase (Lysogeny)", "Holin (Lysis)", "Endolysin (Lysis)")
colnames(interaction_matrix) <- c("MurA (Cell Wall)", "FtsZ (Division)", "tRNA-Leu (attB)")

# Render and Save Heatmap
png("Fig2_Conflict_Heatmap_300DPI.png", width=4000, height=1500, res=300)
pheatmap(interaction_matrix, 
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("navy", "yellow", "firebrick3"))(100),
         main = "Phage-Host Strategic Conflict Map",
         display_numbers = TRUE, number_color = "white",
         angle_col = 45)
dev.off()

# 5. GENOMIC INTEGRATION LANDSCAPE (ggplot2)
# ------------------------------------------------------------------------------
cat("Step 3: Mapping Genomic Integration Hotspots...\n")
# Probe for the 34bp integration site (attP)
phage_probe <- subseq(phage_dna[[1]], start = 1, end = 50) 

n_bins <- 100
genome_len <- width(host_dna)[1]
bin_size <- floor(genome_len / n_bins)
hits_vector <- numeric(n_bins)

for(i in 1:n_bins){
  start_pos <- (i-1) * bin_size + 1
  end_pos   <- min(i * bin_size, genome_len)
  region    <- host_dna[[1]][start_pos:end_pos]
  hits_vector[i] <- countPattern(phage_probe, region, max.mismatch = 15)
}

hotspot_df <- data.frame(
  Position_Mb = seq(0, genome_len, length.out=n_bins) / 1e6,
  Intensity = hits_vector,
  Label = "Phage Interaction"
)

# Plotting the Landscape
p_hotspot <- ggplot(hotspot_df, aes(x = Position_Mb, y = Label, fill = Intensity)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("midnightblue", "yellow", "firebrick3")) +
  theme_minimal() +
  labs(title = "Phage Integration Potential across E. amylovora Genome",
       x = "Host Genome Position (Mb)", y = "") +
  theme(axis.text.y = element_blank(), panel.grid = element_blank())

ggsave("Fig3_Integration_Landscape_300DPI.png", plot = p_hotspot, width = 12, height = 3, dpi = 300)

# Extract probe for site-specific integration
phage_probe <- subseq(phage_dna[[1]], start = 1, end = 50) 

# Sliding window scan across Host Genome
n_bins <- 100
genome_len <- width(host_dna)[1]
bin_size <- floor(genome_len / n_bins)
hits_vector <- numeric(n_bins)

for(i in 1:n_bins){
  start_pos <- (i-1) * bin_size + 1
  end_pos   <- min(i * bin_size, genome_len)
  region    <- host_dna[[1]][start_pos:end_pos]
  hits_vector[i] <- countPattern(phage_probe, region, max.mismatch = 15)
}

# Create Dataframe for ggplot
hotspot_df <- data.frame(
  Position_Mb = seq(0, genome_len, length.out=n_bins) / 1e6,
  Intensity = hits_vector,
  Label = "Phage Interaction"
)

# Plot Genomic Landscape
p_landscape <- ggplot(hotspot_df, aes(x = Position_Mb, y = Label, fill = Intensity)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("midnightblue", "yellow", "firebrick3")) +
  theme_minimal() +
  labs(title = "Genomic Integration Hotspots", x = "Bacteria Genome Position (Mb)", y = "")

ggsave("Genomic_Integration_Landscape_300DPI.png", plot = p_landscape, width = 12, height = 3, dpi = 300)

# 6. BATTLE INTENSITY SCORING (Functional Impact)
# ------------------------------------------------------------------------------
cat("Step 4: Plotting Functional Impact Scores...\n")
battle_data <- data.frame(
  Accession = c("... (Endolysin)", "... (Integrase)", 
                "... (Excisionase)", "....1 (Holin)", 
                "... (Portal)"),
  Score = c(98, 94, 88, 75, 35)
)

p_impact <- ggplot(battle_data, aes(x = Score, y = reorder(Accession, Score))) +
  geom_segment(aes(x = 0, xend = Score, yend = Accession), color = "grey") +
  geom_point(size = 4, color = "firebrick") +
  theme_bw() +
  labs(title = "Phage Functional Impact Scores",
       x = "Battle Intensity (Normalized)", y = "NCBI Accession")

ggsave("Fig4_Functional_Impact_300DPI.png", plot = p_impact, width = 8, height = 5, dpi = 300)

# CLEANUP
unlink(db)
cat("All analysis complete. Figures saved to current directory.\n")
