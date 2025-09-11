import pandas as pd
import math as m

df = pd.read_csv('../systematic_study/comp_Alu.csv')

n_variations = len(df.index) / 3 # nominal is considered 1st variation

t_bin_indices = [0, 1, 2]

nom_vals = {}

for i in t_bin_indices: # collect nominal valuse
    nom_vals[i] = df.iloc[i, 4]


rms_sum = [0.0, 0.0, 0.0]


for i in t_bin_indices:
    for index, row in df.iterrows():
        if index % 3 == i and index < 3:
            continue
    
        if index % 3 == i:
            rms_sum[i] += (nom_vals[i] - row['LTp']) ** 2


dLTp = [
    m.sqrt(rms_sum[0] / n_variations),
    m.sqrt(rms_sum[1] / n_variations),
    m.sqrt(rms_sum[2] / n_variations)
    ]
print("dLTp_cuts:", dLTp)