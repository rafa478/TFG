import pandas as pd
import matplotlib.pyplot as plt
import sys

files = []
for file in sys.argv[1:]:
  if 'media' not in file:
    files.append(file)
print(files)


ninguna = []
muertas = []
metastasis = []
metastasis_y_proliferacion = []


for file in files:
    df = pd.read_csv(file)
    
    df_final = df.iloc[-1]
    
    ninguna.append(df_final['Ninguna'])
    muertas.append(df_final['Muertas'])
    metastasis.append(df_final['Metastasis'])
    metastasis_y_proliferacion.append(df_final['MetastasisYProliferacion'])

plt.figure(figsize=(8, 6))
plt.boxplot([ninguna, muertas, metastasis], labels=['Ninguna', 'Muertas', 'Metastasis'])
plt.ylabel('Número de células')

plt.ylim(min(ninguna + muertas + metastasis) - 20, max(ninguna + muertas + metastasis) + 20)

plt.show()

plt.figure(figsize=(8, 6))
plt.boxplot([metastasis_y_proliferacion], labels=['Metastasis y Proliferacion'])
plt.ylabel('Número de células')


plt.ylim(min(metastasis_y_proliferacion) - 200, max(metastasis_y_proliferacion) + 200)

plt.show()
