import csv
import matplotlib.pyplot as plt
import sys

archivo_csv = sys.argv[1]

tiempos = []
proliferacion = []
metastasis = []
metyprolif = []
ninguna = []
muertes = []

with open(archivo_csv, 'r') as file_csv:
    lector_csv = csv.reader(file_csv)
    next(lector_csv)
    for row in lector_csv:
        if len(row) > 4:
            tiempos.append(row[0])
            proliferacion.append(float(row[1]))
            metastasis.append(float(row[2]))
            metyprolif.append(float(row[3]))
            ninguna.append(float(row[4]))
            muertes.append(float(row[5]))

plt.figure(figsize=(10, 6))

plt.plot(tiempos, proliferacion, label='Proliferación', marker='o')
plt.plot(tiempos, metastasis, label='Metástasis', marker='o')
plt.plot(tiempos, metyprolif, label='Metástasis y Proliferación', marker='o')
plt.plot(tiempos, ninguna, label='Ninguna', marker='o')
plt.plot(tiempos, muertes, label='Muertes', marker='o')

plt.xlabel('Tiempo')
plt.ylabel('Cantidad')
plt.title('Evolución de Proliferación, Metástasis, LasDos, Muertes y Ninguna')

plt.legend()
plt.grid(True)
plt.xticks(rotation=45)
plt.tight_layout()

plt.show()
