import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv('dS_values_distances.txt', sep='\t', header=None, names=['Gene1', 'Gene2', 'Distance', 'Scaffold1', 'Scaffold2'])

# Add a new column 'SameChromosome' indicating whether genes are on the same chromosome
data['SameChromosome'] = data['Distance'].notna()

# Count the occurrences of genes on the same chromosome
same_chromosome_counts = data['SameChromosome'].value_counts()

# Bar plot
same_chromosome_counts.plot(kind='bar', color=['red', 'green'])
plt.xticks([0, 1], ['Different Chromosomes', 'Same Chromosome'], rotation=0)
plt.xlabel('Chromosome Relationship')
plt.ylabel('Count')
plt.title('Count of Genes on Different and Same Chromosomes')
plt.show()

# Filter out rows where 'Distance' is NaN
computed_distances = data[data['Distance'].notna()]

# Boxplot of distances for each scaffold
plt.figure(figsize=(12, 6))
ax = sns.boxplot(x='Scaffold1', y='Distance', data=computed_distances)
ax.set_yscale('log')
plt.xlabel('Scaffold')
plt.ylabel('Distance')
plt.title('Distribution of Distances for Each Scaffold')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()