import numpy as np
import pandas as pd
import json
from rdkit import Chem
from rdkit.Chem import PandasTools
from sklearn.manifold import TSNE
import matplotlib
import matplotlib.pyplot as plt

# Load data from JSON file
with open("resources/CBLB_inhibitors_vsF_added_embedding_mean.json", "r") as file:
    json_data = json.load(file)

df = pd.DataFrame(json_data)
embeddings = np.array(df["embedding_mean"].tolist())  # embeddings -> embeddings_mean へマニュアルで変更
IC50_range_nM = df["IC50_range_nM"].tolist()

# Define Active/Inactive
strong_actives = ['<100', '101-300']
weak_actives = ['301-1000', '101-1000']
df['active'] = df['IC50_range_nM'].apply(lambda x: 2 if x in strong_actives else 1 if x in weak_actives else 0)

# Convert to tsne space
tsne = TSNE(n_components=2, perplexity=30, early_exaggeration=30, random_state=42, init='random', learning_rate=200)
vis_dims = tsne.fit_transform(embeddings)

# Plot
x = [x for x,y in vis_dims]
y = [y for x,y in vis_dims]

color_indices = df['active']
colors = ["darkgreen", "gold", "red"]
colormap = matplotlib.colors.ListedColormap(colors)

plt.scatter(x, y, c=color_indices, cmap=colormap, alpha=0.3)

plt.scatter([],[], color=colormap(0), label='inactive', alpha=0.5)
plt.scatter([],[], color=colormap(1), label='week-active', alpha=0.5)
plt.scatter([],[], color=colormap(2), label='active', alpha=0.5)

plt.legend()
plt.title("Embedded SMILES with text-embedding-ada-002 using t-SNE")
plt.savefig("resources/CBLB_inhibitors_vsF_added_embedding_mean.png")  # マニュアルで追加
