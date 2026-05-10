import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import os

np.random.seed(123)

D0, K0 = 10, 5

def rspdmatrix(D, sparsity):
    """Generate a random sparse positive definite matrix"""
    A = np.random.uniform(0, 1, size=(D, D))
    mask = np.random.binomial(1, 1 - sparsity, size=(D, D))
    A = A * mask
    A = (A + A.T) / 2
    A = A + (abs(np.min(np.linalg.eigvals(A))) + 1) * np.eye(D)
    return A

# Generate covariance matrices with increasing sparsity
Sigma0 = [rspdmatrix(D0, 0.19 * k) for k in range(1, K0 + 1)]

# Get precision matrices (inverses)
PMat0 = [np.linalg.inv(S) for S in Sigma0]

# Convert to adjacency matrices (non-zero off-diagonal elements)
AMat0 = []
for P in PMat0:
    A = (P != 0).astype(int)
    np.fill_diagonal(A, 0)
    AMat0.append(A)

# Create output directory
os.makedirs('figures', exist_ok=True)

# Create low sparsity figure
fig, axes = plt.subplots(1, 2, figsize=(14, 7))
fig.patch.set_facecolor('white')
fig.suptitle('Low Sparsity Precision Networks\n(More Highly Connected)', fontsize=16, fontweight='bold', y=0.98)

for idx, i in enumerate([1, 2]):
    ax = axes[idx]
    ax.set_facecolor('white')
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.15, 1.15)
    ax.set_aspect('equal')
    ax.axis('off')
    
    G = nx.Graph(AMat0[i])
    pos = nx.circular_layout(G)
    
    # Draw edges manually for perfect alignment
    for u, v in G.edges():
        x = [pos[u][0], pos[v][0]]
        y = [pos[u][1], pos[v][1]]
        ax.plot(x, y, color='#2E86AB', linewidth=2.5, alpha=0.7, zorder=1)
    
    # Draw nodes
    node_x = [pos[node][0] for node in G.nodes()]
    node_y = [pos[node][1] for node in G.nodes()]
    ax.scatter(node_x, node_y, c='#A23B72', s=1200, edgecolors='black', linewidths=2.5, zorder=10)
    
    # Draw labels
    for node in G.nodes():
        ax.text(pos[node][0], pos[node][1], str(node), ha='center', va='center',
               fontsize=11, fontweight='bold', color='white', zorder=11)
    
    sparsity_pct = 0.19 * (i + 1) * 100
    ax.set_title(f'Sparsity = {sparsity_pct:.1f}% | {G.number_of_edges()} edges', 
                fontsize=12, fontweight='bold', pad=10)

plt.subplots_adjust(top=0.92, bottom=0.05, left=0.05, right=0.95, hspace=0.3, wspace=0.3)
plt.savefig('figures/precision_networks_low_sparsity.png', dpi=150, bbox_inches='tight', 
           facecolor='white', edgecolor='none')
print('Saved: figures/precision_networks_low_sparsity.png')

# Create high sparsity figure
fig, axes = plt.subplots(1, 2, figsize=(14, 7))
fig.patch.set_facecolor('white')
fig.suptitle('High Sparsity Precision Networks\n(Fewer Connections)', fontsize=16, fontweight='bold', y=0.98)

for idx, i in enumerate([3, 4]):
    ax = axes[idx]
    ax.set_facecolor('white')
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.15, 1.15)
    ax.set_aspect('equal')
    ax.axis('off')
    
    G = nx.Graph(AMat0[i])
    pos = nx.circular_layout(G)
    
    # Draw edges manually for perfect alignment
    for u, v in G.edges():
        x = [pos[u][0], pos[v][0]]
        y = [pos[u][1], pos[v][1]]
        ax.plot(x, y, color='#F18F01', linewidth=2.5, alpha=0.7, zorder=1)
    
    # Draw nodes
    node_x = [pos[node][0] for node in G.nodes()]
    node_y = [pos[node][1] for node in G.nodes()]
    ax.scatter(node_x, node_y, c='#C73E1D', s=1200, edgecolors='black', linewidths=2.5, zorder=10)
    
    # Draw labels
    for node in G.nodes():
        ax.text(pos[node][0], pos[node][1], str(node), ha='center', va='center',
               fontsize=11, fontweight='bold', color='white', zorder=11)
    
    sparsity_pct = 0.19 * (i + 1) * 100
    ax.set_title(f'Sparsity = {sparsity_pct:.1f}% | {G.number_of_edges()} edges', 
                fontsize=12, fontweight='bold', pad=10)

plt.subplots_adjust(top=0.92, bottom=0.05, left=0.05, right=0.95, hspace=0.3, wspace=0.3)
plt.savefig('figures/precision_networks_high_sparsity.png', dpi=150, bbox_inches='tight', 
           facecolor='white', edgecolor='none')
print('Saved: figures/precision_networks_high_sparsity.png')

print('Images created successfully!')
