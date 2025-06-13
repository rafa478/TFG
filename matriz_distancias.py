import pandas as pd
import sys

data = pd.read_csv(sys.argv[1])

WT = data.iloc[-1, 1:].values

distancias = []

for i in range(len(data) - 1):
    distancias.append(data.iloc[i, 1:].values - WT)

distancias_df = pd.DataFrame(distancias, columns=data.columns[1:], index=data['Gen'][:-1])

result_df = data.copy()
result_df = result_df.iloc[:-1]
result_df.iloc[:, 1:] = distancias_df.round(2)

result_df_sorted = result_df.iloc[result_df['Metástasis y Proliferación'].argsort()]

result_df_sorted.to_csv('matriz_dist.csv', index=False)
