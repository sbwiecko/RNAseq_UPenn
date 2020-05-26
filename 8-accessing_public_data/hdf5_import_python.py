#%%
import h5py

#%%
f = h5py.File("mouse_matrix_v8.h5", 'r')

# %%
f

# %%
list(f)

# %%
list(f['meta'])

# %%
data = f['meta']['Sample_geo_accession']
print(data[:10])

# %%
f.close()

# %%
