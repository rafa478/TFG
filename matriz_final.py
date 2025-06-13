import pandas as pd
import os

dir = 'RESULTADOS/Single_KO'

archivos_g3 = [f for f in os.listdir(dir) if f.startswith('G3')]
archivos_g4 = [f for f in os.listdir(dir) if f.startswith('G4')]

def procesar_archivos(archivos, directorio):
    datos = []
    for archivo in archivos:
        df = pd.read_csv(os.path.join(directorio, archivo))
        ultima_fila = df.iloc[-1]
        gen = archivo.replace('G3_KO_', '').replace('G4_KO_', '').replace('.csv', '')
        fila = [gen, ultima_fila['Proliferacion'], ultima_fila['Metastasis'], 
                ultima_fila['MetastasisYProliferacion'], ultima_fila['Ninguna'], ultima_fila['Muertas']]
        datos.append(fila)
    
    return datos

datos_g3 = procesar_archivos(archivos_g3, dir)
datos_g4 = procesar_archivos(archivos_g4, dir)

archivo_wt_G3 = 'RESULTADOS/WT/G3_media.csv'
archivo_wt_G4 = 'RESULTADOS/WT/G4_media.csv'

df_wt = pd.read_csv(archivo_wt_G3)
ultima_fila_wt = df_wt.iloc[-1]
wt_fila_G3 = ['WT', ultima_fila_wt['Proliferacion'], ultima_fila_wt['Metastasis'], 
           ultima_fila_wt['MetastasisYProliferacion'], ultima_fila_wt['Ninguna'], ultima_fila_wt['Muertas']]

df_wt = pd.read_csv(archivo_wt_G4)
ultima_fila_wt = df_wt.iloc[-1]
wt_fila_G4 = ['WT', ultima_fila_wt['Proliferacion'], ultima_fila_wt['Metastasis'], 
           ultima_fila_wt['MetastasisYProliferacion'], ultima_fila_wt['Ninguna'], ultima_fila_wt['Muertas']]

datos_g3.append(wt_fila_G3)
datos_g4.append(wt_fila_G4)

columnas = ['Gen', 'Proliferación', 'Metástasis', 'Metástasis y Proliferación', 'Ninguna', 'Muertas']
df_g3 = pd.DataFrame(datos_g3, columns=columnas)
df_g4 = pd.DataFrame(datos_g4, columns=columnas)

df_g3.to_csv('matriz_G3.csv', index=False)
df_g4.to_csv('matriz_G4.csv', index=False)
