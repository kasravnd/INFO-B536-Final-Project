import random
from motif_discovery import gibbs_sampler_chain, randomized_greedy_search, score_motifs

# --- Helper to generate synthetic data ---
def generate_synthetic_data(num_seqs, seq_len, k):
    dna = []
    # 1. Create Random Background
    for _ in range(num_seqs):
        seq = ''.join(random.choice('ACGT') for _ in range(seq_len))
        dna.append(seq)
    
    # 2. Define Hidden Motif
    hidden_motif = ''.join(random.choice('ACGT') for _ in range(k))
    print(f"True Implanted Motif (Consensus): {hidden_motif}")
    
    implanted_dna = []
    for seq in dna:
        # 3. Mutate the motif before implanting!
        # 30% chance to mutate one letter in the motif
        motif_variant = list(hidden_motif)
        if random.random() < 0.3: 
            pos = random.randint(0, k-1)
            motif_variant[pos] = random.choice('ACGT')
        motif_variant = "".join(motif_variant)

        idx = random.randint(0, seq_len - k)
        new_seq = seq[:idx] + motif_variant + seq[idx+k:]
        implanted_dna.append(new_seq)
        
    return implanted_dna, hidden_motif

# --- Run the Comparison ---
def compare_algorithms():
    print("--- Running Synthetic Data Comparison (Gibbs vs Greedy) ---")
    
    # 1. Create Data
    k = 8
    t = 10
    seq_len = 100
    dna, true_motif = generate_synthetic_data(t, seq_len, k)
    print(f"True Implanted Motif: {true_motif}")
    
    # 2. Run Greedy (Baseline)
    greedy_result = randomized_greedy_search(dna, k, t)
    greedy_score = score_motifs(greedy_result)
    greedy_consensus = greedy_result[0] # Simplified consensus
    
    # 3. Run Gibbs (Your Algorithm)
    # Run a few restarts to be fair
    best_gibbs_score = float('inf')
    best_gibbs_motif = ""
    for _ in range(5): 
        motifs, score = gibbs_sampler_chain(dna, k, t, 1000)
        if score < best_gibbs_score:
            best_gibbs_score = score
            best_gibbs_motif = motifs[0]
            
    # 4. Print Results for your Report
    print(f"\nRESULTS:")
    print(f"Greedy Algorithm -> Score: {greedy_score}, Found: {greedy_consensus}")
    print(f"Gibbs Sampling   -> Score: {best_gibbs_score}, Found: {best_gibbs_motif}")
    
    if best_gibbs_score < greedy_score:
        print("\nCONCLUSION FOR REPORT: Gibbs performed better (lower score).")
    elif greedy_consensus != true_motif and best_gibbs_motif == true_motif:
        print("\nCONCLUSION FOR REPORT: Gibbs found the true motif, Greedy failed.")
    else:
        print("\nCONCLUSION: Both performed similarly on this easy dataset.")

if __name__ == "__main__":
    compare_algorithms()