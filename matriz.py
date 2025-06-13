import os
import csv
import scipy.io
import h5py

directorio = 'output'

archivos_csv = [file for file in os.listdir(directorio) if file.endswith('_boolean_intracellular.csv')]
archivos_csv.sort()
archivos_csv.remove('initial_boolean_intracellular.csv')
archivos_csv.remove('final_boolean_intracellular.csv')

archivos_mat = [file for file in os.listdir(directorio) if file.endswith('_cells.mat')]
archivos_mat.sort()
archivos_mat.remove('initial_cells.mat')
archivos_mat.remove('final_cells.mat')
archivos_mat = ['initial_cells.mat'] + archivos_mat + ['final_cells.mat']

resultados = {}

def convertir_a_hdf5(archivo_mat):
    try:
        with h5py.File(archivo_mat, 'r'):
            return archivo_mat
    except:
        mat_data = scipy.io.loadmat(archivo_mat)
        nuevo_archivo_mat = archivo_mat.replace('.mat', '_hdf5.mat')
        
        with h5py.File(nuevo_archivo_mat, 'w') as f:
            for key in mat_data:
                f.create_dataset(key, data=mat_data[key])
                
        return nuevo_archivo_mat

contador_muertes = 0

for archivo_mat, archivo_csv in zip(archivos_mat, archivos_csv):
    archivo_mat = os.path.join(directorio, archivo_mat)
    archivo_csv = os.path.join(directorio, archivo_csv)
    
    archivo_mat_hdf5 = convertir_a_hdf5(archivo_mat)
    
    proliferacion = 0
    metastasis = 0
    metyprolif = 0
    ninguna = 0
    muertes = 0

    with open(archivo_csv, 'r') as file_csv:
        lector_csv = csv.reader(file_csv)       
        for line in lector_csv:
            line = [campo.strip().lower() for campo in line]
            if line[0] == 'id':
                continue
            line = line[1].split(' -- ')
            if 'metastasis' and 'proliferation' in line:
                metyprolif += 1
            elif 'metastasis' in line:
                metastasis += 1
            elif 'proliferation' in line:
                proliferacion += 1
            else:
                ninguna += 1

    resultados[archivo_csv] = [proliferacion, metastasis, metyprolif, ninguna, contador_muertes]

tiempo = 0
with open('RESULTADOS/matriz.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Tiempo', 'Proliferacion', 'Metastasis', 'MetastasisYProliferacion', 'Ninguna', 'Muertas'])
    for archivo, conteos in resultados.items():
        if archivo == 'output/initial_boolean_intracellular.csv' or archivo == 'output/final_boolean_intracellular.csv':
            writer.writerow([archivo, conteos[0], conteos[1], conteos[2], conteos[3], int(conteos[4])])
        else:
            writer.writerow([str(tiempo), conteos[0], conteos[1], conteos[2], conteos[3], int(conteos[4])])
            tiempo += 1
