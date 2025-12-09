import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from io import StringIO

data_str_test = """Northing,Easting,Depth
9951,15731,8141
9952,15734,8143
10545,15254,8412
10025,15889,8185
9861,14614,8358"""

df = pd.read_csv("data/event_data.csv")
df.drop('Moment Magnitude', axis=1, inplace=True)

print("Original Data:")
print(df)

# Standardize the Data
# This centers the data at 0 and scales it to unit variance
scaler = StandardScaler()
data_scaled = scaler.fit_transform(df)

pca = PCA(n_components=3)
pca.fit(data_scaled)

components = pca.components_           # The eigenvectors
explained_variance = pca.explained_variance_ # The eigenvalues

print("\nEigenvectors (Principal Components):")
print(components)
print("\nExplained Variance Ratio:")
print(pca.explained_variance_ratio_)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot the standardized data points
ax.scatter(data_scaled[:,0], data_scaled[:,1], data_scaled[:,2], 
           c='b', marker='o', alpha=0.6, label='Data Points')

# Plot the eigenvectors
origin = np.zeros(3)
colors = ['r', 'g', 'orange']
labels = ['PC1', 'PC2', 'PC3']

scale_factor = 3

for i in range(3):
    vector = components[i] * np.sqrt(explained_variance[i]) * scale_factor
    print("PCA 1", vector)
    
    ax.quiver(origin[0], origin[1], origin[2], 
              vector[0], vector[1], vector[2], 
              color=colors[i], label=f'{labels[i]} (Var: {pca.explained_variance_ratio_[i]:.2f})', 
              linewidth=2.5)

ax.set_xlabel('Northing (Standardized)')
ax.set_ylabel('Easting (Standardized)')
ax.set_zlabel('Depth (Standardized)')
ax.set_title('3D PCA: Data and Eigenvectors')
ax.legend()
plt.show()