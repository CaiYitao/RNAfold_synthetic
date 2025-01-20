import varnaapi as api
import pandas as pd
import RNA

rnas = pd.read_csv("/home/mescalin/yitao/Documents/Code/RNAfold/rna_dataset.csv")
print(rnas)

rna = rnas.to_numpy()[0]

# rnafig = api.Structure(rna[0], rna[1])
# rnafig.show()
RNA.plot(rna[0], rna[1])
