import random
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from motif_discovery import gibbs_sampler_chain, randomized_greedy_search, score_motifs

# ==========================================
# 1. Setup Data Generator (with Mutation)
# ==========================================
def generate_synthetic_data(num_seqs, seq_len, k):
    dna = []
    # Background
    for _ in range(num_seqs):
        seq = ''.join(random.choice('ACGT') for _ in range(seq_len))
        dna.append(seq)
    
    # Hidden Motif
    hidden_motif = ''.join(random.choice('ACGT') for _ in range(k))
    
    implanted_dna = []
    for seq in dna:
        # 30% Mutation Rate (Simulation of Biological Noise)
        motif_variant = list(hidden_motif)
        if random.random() < 0.3: 
            pos = random.randint(0, k-1)
            motif_variant[pos] = random.choice('ACGT')
        motif_variant = "".join(motif_variant)

        idx = random.randint(0, seq_len - k)
        new_seq = seq[:idx] + motif_variant + seq[idx+k:]
        implanted_dna.append(new_seq)
        
    return implanted_dna, hidden_motif

# ==========================================
# 2. Run Statistical Simulation
# ==========================================
def run_simulation_and_plot():
    TRIALS = 50  # 50 trials is enough for a good plot (and faster)
    k = 8
    t = 10
    seq_len = 100
    
    # Data storage for Box Plot
    results_data = [] 
    
    # Data storage for Success Rate
    success_counts = {"Gibbs": 0, "Greedy": 0}

    print(f"--- Running {TRIALS} Statistical Trials for Visualization ---")
    
    for i in tqdm(range(TRIALS)):
        dna, true_motif = generate_synthetic_data(t, seq_len, k)
        
        # --- A. Run Greedy ---
        greedy_res = randomized_greedy_search(dna, k, t)
        greedy_score = score_motifs(greedy_res)
        
        # Store for Box Plot
        results_data.append({"Algorithm": "Greedy", "Score": greedy_score})
        
        # Check Success
        if greedy_res[0] == true_motif:
            success_counts["Greedy"] += 1
            
        # --- B. Run Gibbs (Best of 5 Restarts) ---
        best_gibbs_score = float('inf')
        best_gibbs_motif = ""
        for _ in range(5):
            motifs, score = gibbs_sampler_chain(dna, k, t, 1000)
            if score < best_gibbs_score:
                best_gibbs_score = score
                best_gibbs_motif = motifs[0]
        
        # Store for Box Plot
        results_data.append({"Algorithm": "Gibbs", "Score": best_gibbs_score})
        
        # Check Success
        if best_gibbs_motif == true_motif:
            success_counts["Gibbs"] += 1

    # Convert to DataFrame for Seaborn
    df_scores = pd.DataFrame(results_data)

    # ==========================================
    # 3. Generate Plot 1: Score Distribution (Box Plot)
    # ==========================================
    plt.figure(figsize=(8, 6))
    sns.boxplot(x="Algorithm", y="Score", data=df_scores, palette="Set2")
    plt.title(f"Score Distribution Comparison (Over {TRIALS} Trials)", fontsize=14)
    plt.ylabel("Motif Score (Lower is Better)", fontsize=12)
    plt.xlabel("Algorithm", fontsize=12)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Save
    plt.savefig("comparison_boxplot.png", dpi=300)
    print("\nSaved Box Plot to 'comparison_boxplot.png'")

    # ==========================================
    # 4. Generate Plot 2: Success Rate (Bar Chart)
    # ==========================================
    algorithms = list(success_counts.keys())
    success_rates = [(success_counts[algo] / TRIALS) * 100 for algo in algorithms]

    plt.figure(figsize=(8, 6))
    colors = ['#FF9999', '#66B2FF'] # Red-ish for Greedy, Blue-ish for Gibbs
    plt.bar(algorithms, success_rates, color=colors, alpha=0.9, width=0.6)
    
    plt.title(f"Algorithm Success Rate (Exact Motif Recovery)", fontsize=14)
    plt.ylabel("Success Rate (%)", fontsize=12)
    plt.ylim(0, 100)
    
    # Add percentage labels on top of bars
    for i, v in enumerate(success_rates):
        plt.text(i, v + 2, f"{v:.1f}%", ha='center', fontsize=12, fontweight='bold')

    # Save
    plt.savefig("comparison_success_rate.png", dpi=300)
    print("Saved Success Rate Chart to 'comparison_success_rate.png'")

if __name__ == "__main__":
    run_simulation_and_plot()