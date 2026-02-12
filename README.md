# Phage-Bacteria_Network-
Phage-Host Interaction Pipeline 
# Phage-Host Conflict Map: *X_bacterialpathogen* vs. Phage Genome

This repository contains the bioinformatic pipeline used to analyze the genomic interactions between a novel phage and its host, *Erwinia amylovora* (Strain EaSmR). 

## Project Overview
The analysis focuses on two primary survival strategies used by the phage:
1. **The Lytic Phase:** Identifying high-efficiency enzymes (Endolysins/Holins) that degrade the host cell wall.
2. **The Lysogenic Phase:** Predicting the site-specific integration of the phage genome into the host chromosome.

##  Key Findings
### 1. Molecular Mimicry (The "Red Hotspot")
Our analysis identified a significant mimicry event at **Host_3 (recF)**. The phage produces a protein mimicking the host's DNA repair machinery, likely to facilitate genome stability during the early stages of infection.

### 2. tRNA Integration (The 34bp "att" Site)
We discovered a perfect **34 bp sequence identity** between the phage and a host **tRNA-Leu gene**. 
- **Sequence:** `ATGCGTACGTTAGCTAGCCTAGCTAGCTAGTATA`
- **Mechanism:** This acts as the *attB* attachment site for site-specific recombination via the phage-encoded Integrase.



### 3. Structural Survival (Electron Microscopy Correlation)
Genomic data correlates with EM observations showing that **flagella remain intact** during lysis. The conflict map confirms that the phage lacks enzymes targeting the flagellin protein (**Host_1 / dnaA zone**), focusing purely on peptidoglycan degradation.


##  How to Run
1. Clone the repo: `git clone https://github.com/yourusername/phage-conflict-map.git`
2. Open `phage_host_analysis.R` in RStudio.
3. Ensure the following packages are installed: `pwalign`, `pheatmap`, `ggplot2`.
4. Provide your Pathogen CDS and Phage FASTA files when prompted.

##  Repository Structure
- `/data`: Raw FASTA files (excluded via .gitignore).
- `/scripts`: R scripts for alignment and clustering.
- `/results`: High-resolution Blue-Yellow-Red heatmaps and integration diagrams.
