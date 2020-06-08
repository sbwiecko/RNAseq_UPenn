# %%
import pandas as pd

# %%
# read in the experimental plan
targets = pd.read_csv("./covid_metadata.txt", sep='\t')
print(targets.columns)

# %%
meta = targets[['sample', 'cell_line', 'treatment', 'timepoint_postTreament']]

meta.groupby(['cell_line', 'treatment', 'timepoint_postTreament']).describe()

# %%
