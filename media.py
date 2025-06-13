import csv
import numpy as np

archivos_G3 = [f'G3_{i}.csv' for i in range(1, 11)]
archivos_G4 = [f'G4_{i}.csv' for i in range(1, 11)]

def leer_y_guardar_media(archivos, nombre_salida):
    datos = []

    for archivo in archivos:
        with open(archivo, 'r') as f:
            lector = csv.reader(f)
            next(lector)
            filas = [list(map(float, row)) for row in lector if len(row) > 5]
            datos.append(filas)

    datos = np.array(datos)
    medias = np.mean(datos, axis=0)

    encabezado = ['Tiempo', 'Proliferacion', 'Metastasis', 'MetastasisYProliferacion', 'Ninguna', 'Muertas']
    with open(nombre_salida, 'w', newline='') as salida:
        escritor = csv.writer(salida)
        escritor.writerow(encabezado)
        for fila in medias:
            escritor.writerow(fila)

leer_y_guardar_media(archivos_G3, 'G3_media.csv')
leer_y_guardar_media(archivos_G4, 'G4_media.csv')
