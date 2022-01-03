import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
import math
import statsmodels.api as sm
import subprocess
import pathlib
import csv
import seaborn as sns

def getValuesFromCoords(path, lonLatList):
    #Coordenadas das estações são traduzidas em pixel da imagem usando 
    # a ferramenta gdallocationinfo com opção -wgs84
    valueList = []
    for lonLat in lonLatList:
        value = str(subprocess.Popen(f"gdallocationinfo -wgs84 -valonly {path} {lonLat[0]} {lonLat[1]}",
                stdout=subprocess.PIPE).stdout.read())
        valueList.append(float(''.join(filter(str.isdigit, value))))
    return valueList

def getToaCalArgs(path):
    dataDict = {}
    with open(path) as f:
        for line in f.readlines():
            if "SUN_ELEVATION" in line:
                dataDict.update({
                    "SUN_ELEVATION": float(str(line).split("= ")[1])
                })
            elif "REFLECTANCE_MULT_BAND_1" in line:
                dataDict.update({
                    "REFLECTANCE_MULT_BAND_1": float(str(line).split("= ")[1])
                })
            elif "REFLECTANCE_MULT_BAND_2" in line:
                dataDict.update({
                    "REFLECTANCE_MULT_BAND_2": float(str(line).split("= ")[1])
                })
            elif "REFLECTANCE_MULT_BAND_3" in line:
                dataDict.update({
                    "REFLECTANCE_MULT_BAND_3": float(str(line).split("= ")[1])
                })
            elif "REFLECTANCE_MULT_BAND_4" in line:
                dataDict.update({
                    "REFLECTANCE_MULT_BAND_4": float(str(line).split("= ")[1])
                })
            elif "REFLECTANCE_ADD_BAND_1" in line:
                dataDict.update({
                    "REFLECTANCE_ADD_BAND_1": float(str(line).split("= ")[1])
                })
            elif "REFLECTANCE_ADD_BAND_2" in line:
                dataDict.update({
                    "REFLECTANCE_ADD_BAND_2": float(str(line).split("= ")[1])
                })
            elif "REFLECTANCE_ADD_BAND_3" in line:
                dataDict.update({
                    "REFLECTANCE_ADD_BAND_3": float(str(line).split("= ")[1])
                })
            elif "REFLECTANCE_ADD_BAND_4" in line:
                dataDict.update({
                    "REFLECTANCE_ADD_BAND_4": float(str(line).split("= ")[1])
                })
    return dataDict


imagePath = str(pathlib.Path(__file__).parent.resolve()) + "\\"
level1identifier = "LC08_L1TP_"
level2identifier = "LC08_L2SP_"

#Dados de média móvel de PM2.5 das estações no mesmo dia e hora da obtenção das imagens, são
# carregados de um arquivo csv, para fácil utilização. Além disso, o arquivo de csv tambem
# inclui a latitude e longitude em wgs84 e o levelString, equivalente ao nome comum entre os arquivos
coords = []
pm = []
levelString = []


with open(imagePath + '\\data.csv', newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    pmAux = []
    for row in reader:
        coords.append([float(row['Longitude'].replace(',', '.')), float(row['Latitude'].replace(',', '.'))])
        row.pop('Estação')
        row.pop('Longitude')
        row.pop('Latitude')
        levelString = list(row.keys())
        pmAux.append(list(row.values()))
    for i in range(len(pmAux[0])):
        for pmValue in pmAux:
            pm.append(int(pmValue[i]))



#variaveis para armazenar refletância da ultima imagem calculada
band1ar = np.zeros([]);
band3ar = np.zeros([]);
band2ar = np.zeros([]);
band4ar = np.zeros([]);

refl = []

for i in range(len(levelString)):
    #Obtendo caminho da imagem de nível 1
    level1folder = imagePath + level1identifier + levelString[i] + "\\"
    #Valores da imagem de nível 1, em números digitais
    band1dn = getValuesFromCoords(level1folder + level1identifier + levelString[i] + "_B1.TIF", coords)
    band2dn = getValuesFromCoords(level1folder + level1identifier + levelString[i] + "_B2.TIF", coords)
    band3dn = getValuesFromCoords(level1folder + level1identifier + levelString[i] + "_B3.TIF", coords)
    band4dn = getValuesFromCoords(level1folder + level1identifier + levelString[i] + "_B4.TIF", coords)

    band1DNarray = np.array(band1dn)
    band2DNarray = np.array(band2dn)
    band3DNarray = np.array(band3dn)
    band4DNarray = np.array(band4dn)
    
    #Obtendo dados necessários para calculo de TOA Reflectance nos pontos
    toaCalArgs = getToaCalArgs(level1folder + level1identifier + levelString[i] + "_MTL.txt");

    #Calculo de TOA Reflectance nos pontos
    sunCos = math.cos(toaCalArgs['SUN_ELEVATION']*math.pi/180)
    band1toa = (band1DNarray*toaCalArgs['REFLECTANCE_MULT_BAND_1'] + toaCalArgs['REFLECTANCE_ADD_BAND_1'])/sunCos
    band2toa = (band2DNarray*toaCalArgs['REFLECTANCE_MULT_BAND_2'] + toaCalArgs['REFLECTANCE_ADD_BAND_2'])/sunCos
    band3toa = (band3DNarray*toaCalArgs['REFLECTANCE_MULT_BAND_3'] + toaCalArgs['REFLECTANCE_ADD_BAND_3'])/sunCos
    band4toa = (band4DNarray*toaCalArgs['REFLECTANCE_MULT_BAND_4'] + toaCalArgs['REFLECTANCE_ADD_BAND_4'])/sunCos


    #Obtendo caminho da imagem de nível 2
    level2folder = imagePath + level2identifier + levelString[i] + "\\"
    
    #Valores da imagem de nível 2, em refletância de superfície
    band1sr = getValuesFromCoords(level2folder + level2identifier + levelString[i] + "_SR_B1.TIF", coords)
    band2sr = getValuesFromCoords(level2folder + level2identifier + levelString[i] + "_SR_B2.TIF", coords)
    band3sr = getValuesFromCoords(level2folder + level2identifier + levelString[i] + "_SR_B3.TIF", coords)
    band4sr = getValuesFromCoords(level2folder + level2identifier + levelString[i] + "_SR_B4.TIF", coords)

    band1SRarray = np.array(band1sr)
    band2SRarray = np.array(band2sr)
    band3SRarray = np.array(band3sr)
    band4SRarray = np.array(band4sr)

    #Calculo de Atmospheric Reflectance nos pontos
    band1ar = (band1toa - band1SRarray)
    band2ar = (band2toa - band2SRarray)
    band3ar = (band3toa - band3SRarray)
    band4ar = (band4toa - band4SRarray)

    #Criação de lista de Atmospheric Reflectance nos pontos
    for i in range(band1ar.size):
        refl.append([band1ar[i], band2ar[i], band3ar[i], band4ar[i]])
        

#Criação e impressão do modelo. Dividimos os valores de Atmospheric Reflectance por 1000 para o modelo obter
# menores coeficientes, não afetando o calculo pois dividimos novamente depois, necessário para
# o statsmodels não acusar erro.
modelo = sm.OLS(pm, sm.add_constant(np.array(refl)/1000))
res = modelo.fit()
print(res.summary())

#Calculando a previsão de PM2.5 referente a ultima imagem carregada (ultima coluna do .csv)
predicted = res.params[0]*np.ones(np.shape(band1ar))
predicted += (res.params[1]/1000)*band1ar
predicted += (res.params[2]/1000)*band2ar
predicted += (res.params[3]/1000)*band3ar
predicted += (res.params[4]/1000)*band4ar
predictedArray = np.array(predicted)

#Obtendo dados de PM2.5 medidos em solo referentes a ultima imagem carregada (ultima coluna do .csv)
measured = np.zeros(np.shape(coords))
measuredArray = np.array(pm[len(coords)*(len(levelString)-1) : len(coords)*len(levelString)])

#Calculando RMSE
print("RMSE: " + str(sm.tools.eval_measures.rmse(measuredArray, predictedArray)) + "µg/m³")

#Criação de gráfico da regressão entre dados calculados e medidos
fig, axe = plt.subplots()
axe.set_ylabel("Dados em solo (µg/m³)")
axe.set_xlabel("Dados calculados (µg/m³)")
axe.set_title("Regressão Linear")
sns.regplot(x=predictedArray, y=measuredArray)
plt.show()