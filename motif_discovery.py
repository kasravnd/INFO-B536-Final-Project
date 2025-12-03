import random
import sys
import logging
import argparse
import pandas as pd
from tqdm import tqdm

# Import visualization if available
try:
    from visualize_results import plot_sequence_logo, plot_pssm_heatmap
    VIZ_AVAILABLE = True
except ImportError:
    VIZ_AVAILABLE = False

# ==========================================
# PART 0: Logging & Configuration
# ==========================================
def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S'
    )

# ==========================================
# PART 1: Helper Functions
# ==========================================

def get_profile(motifs):
    k = len(motifs[0])
    profile = {'A': [1.0]*k, 'C': [1.0]*k, 'G': [1.0]*k, 'T': [1.0]*k}
    for motif in motifs:
        for i, char in enumerate(motif):
            profile[char][i] += 1
    for i in range(k):
        col_sum = sum(profile[nuc][i] for nuc in 'ACGT')
        for nuc in 'ACGT':
            profile[nuc][i] /= col_sum
    return profile

def profile_most_probable_kmer(text, k, profile):
    """
    REQUIRED FOR GREEDY ALGORITHM: Finds the most probable k-mer.
    """
    max_prob = -1
    best_kmer = text[:k]
    
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        prob = 1.0
        for j, char in enumerate(kmer):
            prob *= profile[char][j]
        
        if prob > max_prob:
            max_prob = prob
            best_kmer = kmer
    return best_kmer

def profile_randomly_generated_kmer(text, k, profile):
    """
    REQUIRED FOR GIBBS SAMPLER: Stochastically selects a k-mer.
    """
    probabilities = []
    kmers = []
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        prob = 1.0
        for j, char in enumerate(kmer):
            prob *= profile[char][j]
        probabilities.append(prob)
        kmers.append(kmer)
    
    total_prob = sum(probabilities)
    if total_prob == 0:
        return random.choice(kmers)
    
    normalized_probs = [p / total_prob for p in probabilities]
    return random.choices(kmers, weights=normalized_probs, k=1)[0]

def score_motifs(motifs):
    score = 0
    k = len(motifs[0])
    for i in range(k):
        column = [motif[i] for motif in motifs]
        max_count = 0
        for nuc in 'ACGT':
            count = column.count(nuc)
            if count > max_count:
                max_count = count
        score += (len(motifs) - max_count)
    return score

def load_fasta(filename):
    sequences = []
    current_seq = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    if current_seq: sequences.append("".join(current_seq))
                    current_seq = []
                else:
                    current_seq.append(line.upper())
            if current_seq: sequences.append("".join(current_seq))
        logging.info(f"Loaded {len(sequences)} sequences from {filename}")
        return sequences
    except FileNotFoundError:
        logging.error(f"File '{filename}' not found.")
        return []

def save_results_report(filename, csv_filename, best_motifs, best_score, args, t):
    with open(filename, "w") as f:
        f.write("========================================\n")
        f.write(" GIBBS SAMPLING MOTIF DISCOVERY REPORT\n")
        f.write("========================================\n\n")
        f.write(f"Input File:      {args.file}\n")
        f.write(f"Motif Length (k): {args.k}\n")
        f.write(f"Sequences (t):    {t}\n")
        f.write(f"Iterations:       {args.iters}\n")
        f.write(f"Random Restarts:  {args.restarts}\n\n")
        f.write("----------------------------------------\n")
        f.write(f"BEST SCORE:       {best_score}\n")
        f.write(f"CONSENSUS MOTIF:  {best_motifs[0]}\n")
        f.write("----------------------------------------\n")
        f.write("\nTop 10 Motifs Found:\n")
        for i, m in enumerate(best_motifs[:10]):
            f.write(f"{i+1}: {m}\n")
    logging.info(f"Results summary saved to {filename}")

    df = pd.DataFrame(best_motifs, columns=["Motif_Sequence"])
    df.to_csv(csv_filename, index_label="Sequence_ID")
    logging.info(f"Raw motif data saved to {csv_filename}")

# ==========================================
# PART 2: Algorithms
# ==========================================

def randomized_greedy_search(dna, k, t):
    """
    BASELINE ALGORITHM: Deterministic greedy search.
    """
    n = len(dna[0])
    motifs = []
    for seq in dna:
        r = random.randint(0, n - k)
        motifs.append(seq[r:r+k])
        
    best_motifs = motifs[:]
    
    while True:
        profile = get_profile(motifs)
        motifs = []
        for seq in dna:
            motifs.append(profile_most_probable_kmer(seq, k, profile))
            
        if score_motifs(motifs) < score_motifs(best_motifs):
            best_motifs = motifs[:]
        else:
            return best_motifs

def gibbs_sampler_chain(dna, k, t, N):
    """
    Single chain of the Gibbs Sampler.
    """
    n_len = len(dna[0])
    motifs = []
    for seq in dna:
        r = random.randint(0, n_len - k)
        motifs.append(seq[r:r+k])
        
    best_motifs = motifs[:]
    best_score = score_motifs(motifs)
    
    for _ in range(N):
        i = random.randint(0, t - 1)
        current_motifs_except_i = [motifs[j] for j in range(t) if j != i]
        profile = get_profile(current_motifs_except_i)
        motifs[i] = profile_randomly_generated_kmer(dna[i], k, profile)
        
        current_score = score_motifs(motifs)
        if current_score < best_score:
            best_motifs = motifs[:]
            best_score = current_score
            
    return best_motifs, best_score

def run_gibbs_with_restarts(dna, k, t, N, restarts):
    global_best_motifs = []
    global_best_score = float('inf')
    
    logging.info(f"Starting Gibbs Sampling: k={k}, steps={N}, restarts={restarts}")
    
    for i in tqdm(range(restarts), desc="Gibbs Restarts", unit="chain"):
        motifs, score = gibbs_sampler_chain(dna, k, t, N)
        if score < global_best_score:
            global_best_score = score
            global_best_motifs = motifs
            logging.debug(f"New best score found: {score}")

    return global_best_motifs, global_best_score

# ==========================================
# PART 3: Main Execution
# ==========================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, default="MA0018.4.sites", help="Input FASTA file")
    parser.add_argument("--k", type=int, default=8, help="Length of the motif")
    parser.add_argument("--iters", type=int, default=1000, help="Inner Gibbs iterations")
    parser.add_argument("--restarts", type=int, default=20, help="Random restarts")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging")
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    
    dna = load_fasta(args.file)
    if not dna: sys.exit(1)
        
    t = len(dna)
    
    best_motifs, best_score = run_gibbs_with_restarts(dna, args.k, t, args.iters, args.restarts)
    
    save_results_report("results_summary.txt", "found_motifs.csv", best_motifs, best_score, args, t)
    
    if VIZ_AVAILABLE:
        logging.info("Generating Visualizations...")
        plot_sequence_logo(best_motifs, title=f"Sequence Logo (k={args.k})", save_path="final_sequence_logo.png")
        plot_pssm_heatmap(best_motifs, title=f"PSSM Heatmap (k={args.k})", save_path="final_pssm_heatmap.png")
    else:
        logging.warning("Visualization module missing. No images saved.")

if __name__ == "__main__":
    main()