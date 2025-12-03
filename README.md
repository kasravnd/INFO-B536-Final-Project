# Gibbs Sampling for DNA Motif Discovery

**Course:** B536 - Computational Biology  

## Project Overview

This project addresses the computational challenge of identifying transcription factor binding sites (TFBS) by implementing and evaluating the **Gibbs Sampling algorithm**. We compared this stochastic method against a deterministic **Greedy algorithm** to demonstrate the importance of avoiding local optima in genomic search spaces.

Our implementation was developed *de novo* in Python and rigorously tested on:
1.  **Synthetic Data:** Random DNA with mutated implanted motifs.
2.  **Biological Data:** 121,197 promoter sequences from the JASPAR database (CREB1).

## Key Results

Our comparative analysis yielded the following performance metrics:

| Algorithm | Success Rate | Median Score | Behavior |
|-----------|--------------|--------------|----------|
| **Greedy** | 2.0% | High (~26) | Trapped in local optima |
| **Gibbs** | **56.0%** | Low (~3) | Successfully explores solution space |

## Project Structure & File Descriptions

### Source Code
* **`motif_discovery.py`**: The main entry point for the project. It handles:
    * Loading real biological data (FASTA format).
    * Running the **Gibbs Sampling** algorithm with random restarts.
    * Running the **Greedy** algorithm (for comparison).
    * Generating the consensus motif and score metrics.
* **`visualize_comparison.py`**: The statistical benchmarking script. It:
    * Generates synthetic DNA with mutated implanted motifs.
    * Runs 50+ independent trials of Gibbs vs. Greedy.
    * Produces the *Success Rate* and *Score Distribution* charts.
* **`visualize_results.py`**: A helper module containing the plotting functions for Sequence Logos and PSSM Heatmaps.
* **`run_comparison.py`**: A lightweight script to run a quick, single-pass comparison between the algorithms without generating full charts.

### Data Sets
* **`MA0018.4.sites`**: The raw input dataset obtained from the [JASPAR Database](https://jaspar.elixir.no/). It contains 121,197 promoter sequences known to bind the **CREB1** transcription factor.

### Output Files & Visualizations
* **`results_summary.txt`**: A text report generated after running `motif_discovery.py`, containing the best motif sequence, score, and parameters used.
* **`found_motifs.csv`**: A raw data file listing the specific motif instance found in every single DNA sequence.
* **`final_sequence_logo.png`**: Visual representation of the biological motif information content (in bits).
* **`final_pssm_heatmap.png`**: A heatmap showing the probability matrix of the discovered motif.
* **`comparison_success_rate.png`**: Bar chart proving Gibbs (56%) outperforms Greedy (2%) on synthetic data.
* **`comparison_boxplot.png`**: Box plot showing the variance and stability of the algorithms.

### Configuration
* **`requirements.txt`**: A list of Python dependencies (`numpy`, `pandas`, `logomaker`, `seaborn`, `tqdm`) required to run the project.

## Installation

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/YOUR_USERNAME/motif-discovery.git](https://github.com/YOUR_USERNAME/motif-discovery.git)
   cd motif-discovery