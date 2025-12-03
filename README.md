# Gibbs Sampling for DNA Motif Discovery

**Authors:** Kasra Vand, Matthew Ng, Akanksha Tipparti, Upol Chowdhury  
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


## 🛠️ Installation

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/YOUR_USERNAME/motif-discovery.git](https://github.com/YOUR_USERNAME/motif-discovery.git)
   cd motif-discovery