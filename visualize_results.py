import logomaker
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def plot_sequence_logo(motifs, title="Motif Sequence Logo", save_path="motif_logo.png"):
    """
    Generates and SAVES a Sequence Logo.
    """
    # Create counts matrix
    matrix_counts = logomaker.alignment_to_matrix(motifs)
    
    # Transform to Information Content (Bits)
    matrix_info = logomaker.transform_matrix(matrix_counts, from_type='counts', to_type='information')

    # Create a fresh figure
    plt.figure(figsize=(10, 4))
    
    logo = logomaker.Logo(matrix_info, 
                          color_scheme='classic',
                          vpad=.1,
                          fade_probabilities=True,
                          stack_order='small_on_top')

    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    plt.title(title, fontsize=14)
    plt.ylabel("Information (bits)")
    plt.xlabel("Position")
    
    # SAVE instead of show
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close() # Close to free memory and avoid interference
    print(f"Saved Sequence Logo to: {save_path}")

def plot_pssm_heatmap(motifs, title="PSSM Heatmap", save_path="motif_heatmap.png"):
    """
    Generates and SAVES a Heatmap of the PSSM.
    """
    # Create counts matrix
    counts_df = logomaker.alignment_to_matrix(motifs)
    
    # Normalize to get Probabilities
    # Add a tiny epsilon to avoid division by zero errors if a row is empty
    prob_df = counts_df.div(counts_df.sum(axis=1), axis=0).fillna(0)
    
    # Transpose: Nucleotides on Y-axis, Positions on X-axis
    prob_df_T = prob_df.T 

    # Create a fresh figure
    plt.figure(figsize=(12, 5))
    
    sns.heatmap(prob_df_T, 
                annot=True,       # Show numbers
                fmt=".2f",        # 2 decimal places
                cmap="Blues",     # Color scheme
                cbar_kws={'label': 'Probability'})
    
    plt.title(title, fontsize=14)
    plt.xlabel("Position")
    plt.ylabel("Nucleotide")
    
    # SAVE instead of show
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close() # Close to free memory
    print(f"Saved PSSM Heatmap to: {save_path}")

# ==========================================
# 3. Convergence Plot (Gibbs vs Greedy)
# ==========================================
def plot_convergence(gibbs_scores, greedy_scores):
    """
    Plots the score improvement over iterations.
    You need to modify your algorithms to return a list of scores history.
    """
    plt.figure(figsize=(10, 6))
    
    plt.plot(gibbs_scores, label='Gibbs Sampling', color='blue', linewidth=2)
    # Greedy typically flattens out very fast
    plt.plot(greedy_scores, label='Greedy Search', color='red', linestyle='dashed', linewidth=2)
    
    plt.title("Algorithm Convergence Comparison", fontsize=14)
    plt.xlabel("Iterations")
    plt.ylabel("Motif Score (Lower is Better)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

# --- Example Usage ---
if __name__ == "__main__":
    # Mock data representing motifs found by your algorithm
    found_motifs = [
        "ATGCGT", "ATGCGA", "ATGCGT", "ATGCGG", 
        "ATGCGT", "CTGCGT", "ATGAGT", "ATGCGT"
    ]
    
    print("Generating Sequence Logo...")
    plot_sequence_logo(found_motifs, title="Gibbs Sampler: Recovered Motif Logo")
    
    print("Generating PSSM Heatmap...")
    plot_pssm_heatmap(found_motifs, title="Gibbs Sampler: PSSM Heatmap")