import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import cdist
from collections import defaultdict
from matplotlib.patches import Patch
import os
import matplotlib.colors as mcolors


file_dir = r'./test_dimer/dt1000_fmt_raw_value' # Directory containing the feature file
file_name = 'features_ca.csv'
file_path = os.path.join(file_dir, file_name)
output_dir = r'./test_dimer'                    # Directory to save plots

# KMeans and PCA-related parameters
scaling_option = 'normalize'  # Scaling options: 'none', 'standard', 'minmax', 'normalize'
num_clusters = 2 
num_bins = 35
default_pc_x = 1
default_pc_y = 2
pc_for_1d_analysis = 1
num_contributions_to_print = 5
num_residue_contributions_to_print = 5
num_pcs_to_plot = 3
num_residue_contributions_to_plot = 6

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

base_filename = os.path.splitext(os.path.basename(file_name))[0]
var_plot_name = os.path.join(output_dir, f'{base_filename}_variance_explained.png')
res_contribution_plot_name = os.path.join(output_dir, f'{base_filename}_residue_contributions.png')
kmeans_plot_name = os.path.join(output_dir, f'{base_filename}_kmeans.png')
pca_plot_name = os.path.join(output_dir, f'{base_filename}_pca_plot.png')


def get_evenly_spaced_colors(colormap, num_colors):

    if num_colors <= 1:
        return [colormap(0.5)] 
    return [colormap(i / (num_colors - 1)) for i in range(num_colors)]


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n))
    )
    return new_cmap

def plot_variance_explained(explained_variance, var_plot_name='variance_explained.png'):
    plt.figure(figsize=(15, 6))
    x_labels = [f'PC{i}' for i in range(1, len(explained_variance) + 1)]
    
    sns.barplot(x=x_labels, y=explained_variance * 100, palette='coolwarm_r')
    
    plt.xticks(rotation=90)
    
    plt.title('Variance Explained by Principal Components')
    plt.xlabel('Principal Component')
    plt.ylabel('Variance Explained (%)')
    plt.tight_layout()
    plt.savefig(var_plot_name, dpi=600, transparent=True)
    plt.show()


def plot_residue_contributions(aggregated_residues, num_pcs_to_plot, num_residue_contributions_to_plot, \
                               res_contribution_plot_name='residue_contributions.png'):
    fig, ax = plt.subplots(figsize=(12, 7)) 

    bar_width = 0.6         
    spacing = 1.5           
    all_positions = []      
    all_labels = []         

    colormap = plt.get_cmap('viridis') 

    pc_colors = truncate_colormap(colormap, 0.2, 0.45, num_pcs_to_plot)
    pc_colors = truncate_colormap(colormap, 0.2, 0.9, num_pcs_to_plot)

    handles = []

    for i, (pc, residue_contributions) in enumerate(list(aggregated_residues.items())[:num_pcs_to_plot]):
        sorted_residues = sorted(residue_contributions.items(), key=lambda x: x[1], reverse=True)
        top_residues = sorted_residues[:num_residue_contributions_to_plot]
        residue_labels = [residue for residue, _ in top_residues]
        contributions = [contribution for _, contribution in top_residues]

        n = len(top_residues)  

        if n > 1:
            step = 0.8 / (n - 1)  
        else:
            step = 0.0            

        alphas = [0.8 - step * j for j in range(n)]
        alphas = [max(a, 0.1) for a in alphas]  

        pos = np.arange(n) + i * (n + spacing)
        all_positions.extend(pos)
        all_labels.extend(residue_labels)

        color = pc_colors(i / num_pcs_to_plot)  # Assign evenly spaced color to each PC

        for j, (contribution, label) in enumerate(zip(contributions, residue_labels)):
            ax.bar(
                pos[j],
                contribution,
                bar_width,
                color=color,          
                alpha=alphas[j]       
            )

        legend_patch = Patch(facecolor=color, alpha=0.7, label=f'PC{i+1}')
        handles.append(legend_patch)

    ax.legend(handles=handles, loc='best')
    ax.set_xticks(all_positions)
    ax.set_xticklabels(all_labels, rotation=45, ha='right')
    ax.set_ylabel('Contribution Percentage (%)')
    ax.set_title('Residue Contribution per Principal Component')
    plt.tight_layout()
    plt.savefig(res_contribution_plot_name, dpi=600, transparent=True)
    plt.show()


df = pd.read_csv(file_path)
df['index_col'] = df.index

residues_df = df.drop(columns=['index_col', 'filename']).copy()  
residue_names = residues_df.columns.tolist() 

if scaling_option == 'none':
    # No scaling
    scaled_df = residues_df.copy()
elif scaling_option == 'standard':
    # Standard Scaling (mean=0, variance=1)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(residues_df)
    scaled_df = pd.DataFrame(scaled_data, columns=residues_df.columns)
elif scaling_option == 'minmax':
    # Min-Max Scaling (range 0 to 1)
    from sklearn.preprocessing import MinMaxScaler
    scaler = MinMaxScaler()
    scaled_data = scaler.fit_transform(residues_df)
    scaled_df = pd.DataFrame(scaled_data, columns=residues_df.columns)
elif scaling_option == 'normalize':
    # Normalize each row (vector magnitude = 1)
    from sklearn.preprocessing import Normalizer
    scaler = Normalizer()
    scaled_data = scaler.fit_transform(residues_df)
    scaled_df = pd.DataFrame(scaled_data, columns=residues_df.columns)
else:
    raise ValueError(f"Invalid scaling option: {scaling_option}")

pca = PCA()
pca_result = pca.fit_transform(scaled_df)
num_components = pca_result.shape[1]
pca_df = pd.DataFrame(pca_result, columns=[f'PC{i}' for i in range(1, num_components + 1)])


print(f"Number of principal components calculated: {num_components}")

explained_variance = pca.explained_variance_ratio_
plot_variance_explained(pca.explained_variance_ratio_, var_plot_name)

for i, var in enumerate(explained_variance[:num_contributions_to_print], start=1):
    print(f"PC{i} explains {var * 100:.2f}% of the variance")

aggregated_residues = {}
for i, component in enumerate(pca.components_[:num_pcs_to_plot], start=1):
    aggregated_residue_contributions = defaultdict(float)
    for residue, contribution in zip(residue_names, component):
        residue_prefix = "_".join(residue.split('_')[:2]) if '_' in residue else residue
        aggregated_residue_contributions[residue_prefix] += abs(contribution)
    
    total_absolute_contribution = sum(aggregated_residue_contributions.values())
    if total_absolute_contribution > 0:
        normalized_residue_contributions = {
            residue: (contrib / total_absolute_contribution) * 100
            for residue, contrib in aggregated_residue_contributions.items()
        }
    else:
        normalized_residue_contributions = {
            residue: 0.0 for residue in aggregated_residue_contributions
        }
    aggregated_residues[f'PC{i}'] = normalized_residue_contributions

    sorted_residues = sorted(normalized_residue_contributions.items(), key=lambda x: x[1], reverse=True)
    for residue, normalized_contribution in sorted_residues[:num_residue_contributions_to_print]:
        print(f"{residue}: {normalized_contribution:.2f}%")

plot_residue_contributions(aggregated_residues, num_pcs_to_plot, num_residue_contributions_to_plot, res_contribution_plot_name)


df[f'PC{pc_for_1d_analysis}'] = pca_df[f'PC{pc_for_1d_analysis}']
kmeans_data = pca_df[[f'PC{pc_for_1d_analysis}']]  
print(kmeans_data)

# # debug part
# kmeans_data = df.drop(columns=['index_col', 'filename', 'PC1']).copy()
# scaler = StandardScaler()
# kmeans_data_scaled = scaler.fit_transform(kmeans_data)
# print(kmeans_data)

kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(kmeans_data)
df['cluster'] = kmeans.labels_ 


for i in range(num_clusters):
    cluster_avg = df[df['cluster'] == i][f'PC{pc_for_1d_analysis}'].mean()
    print(f"Cluster {i+1} average PC{pc_for_1d_analysis}: {cluster_avg:.4f}")

data_range = df[f'PC{pc_for_1d_analysis}'].max() - df[f'PC{pc_for_1d_analysis}'].min()
bin_width = data_range / num_bins


colormap = plt.get_cmap('coolwarm')  # 'plasma', 'viridis'
cluster_colors = get_evenly_spaced_colors(colormap, num_clusters)

plt.figure(figsize=(8, 6))
bins = np.arange(df[f'PC{pc_for_1d_analysis}'].min(), df[f'PC{pc_for_1d_analysis}'].max() + bin_width, bin_width)

for i in range(num_clusters):
    cluster_data = df[df['cluster'] == i][f'PC{pc_for_1d_analysis}']
    color = cluster_colors[i] 
    sns.histplot(cluster_data, bins=bins, kde=True, label=f'Cluster {i+1}', color=color, alpha=0.5)

plt.xlabel(f'PC{pc_for_1d_analysis}')
plt.ylabel('Frequency')
plt.legend()
plt.savefig(kmeans_plot_name, dpi=600, transparent=True)
plt.show()


def plot_pcs(pca_df, pc_x=default_pc_x, pc_y=default_pc_y, pca_plot_name=pca_plot_name):
    plt.figure(figsize=(8, 6))

    color_values = np.linspace(0, 1, len(pca_df)) 
    scatter = plt.scatter(pca_df[f'PC{pc_x}'], pca_df[f'PC{pc_y}'], c=color_values, cmap='rainbow_r', alpha=0.5)

    cbar = plt.colorbar(scatter, label='Time')
    cbar.set_ticks(np.linspace(0, 1, 6))  # 5 ticks from 0 to 1
    cbar.set_ticklabels(['0 ns', '100 ns', '200 ns', '300 ns', '400 ns', '500 ns'])
    plt.xlabel(f'PC{pc_x}')
    plt.ylabel(f'PC{pc_y}')
    
    plt.savefig(pca_plot_name, dpi=600, transparent=True)
    plt.show()

plot_pcs(pca_df, pc_x=default_pc_x, pc_y=default_pc_y, pca_plot_name=pca_plot_name)

closest_files = {}
for i in range(num_clusters):
    cluster_center = kmeans.cluster_centers_[i]
    cluster_indices = np.where(df['cluster'] == i)[0]
    distances = cdist(kmeans_data.iloc[cluster_indices], [cluster_center])  
    closest_indices = cluster_indices[np.argsort(distances.ravel())[:5]]
    closest_files[i] = df.iloc[closest_indices][['index_col', 'filename']]  

for cluster, files in closest_files.items():
    print(f"\nCluster {cluster + 1} - Closest Original Indices and Filenames to Cluster Center:")
    print(files)
