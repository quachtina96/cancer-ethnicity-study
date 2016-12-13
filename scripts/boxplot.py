import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import seaborn as sns
import csv
from glm import load_data

filename = '../data/BRCA/BRCA'
cancer = 'BRCA'

data, patients, genes, master_df, clinical, indices_to_delete = load_data(filename, cancer)
expression_data = [item for a, item in enumerate(list(data[18286])) if a not in indices_to_delete]
df = master_df.copy()
df['expression'] = pd.Series(expression_data, index=master_df.index)
print df
ax = sns.boxplot(x="race", y="expression", data=df)
sns.plt.show()



