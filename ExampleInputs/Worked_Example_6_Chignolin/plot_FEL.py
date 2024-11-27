import os
from os import path as p

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns


def inputs():
    reportersDir = p.join(os.getcwd(),"outputs", "Chignolin", "06_Metadynamics", "00_reporters_and_plots")
    freeEnergyCsv = p.join(reportersDir, "freeEnergy.csv")
    outDir = os.getcwd()
    distanceRange = (2.7, 10)


    return freeEnergyCsv, outDir, distanceRange

def main():
    freeEnergyCsv, outDir, distanceRange = inputs()


    freeEnergyDf = pd.read_csv(freeEnergyCsv, skiprows=1)

    colNames = np.linspace(2.7, 10, len(freeEnergyDf.columns))
    indexNames = np.linspace(2.7, 10, num=len(freeEnergyDf.index))
    
    freeEnergyDf.columns = colNames
    freeEnergyDf.index = indexNames

    globalMinimumEnergy = freeEnergyDf.min().min()
    freeEnergyDf = freeEnergyDf - globalMinimumEnergy

    freeEnergyDf = freeEnergyDf.iloc[::-1]
    # Assuming freeEnergyDf is your DataFrame
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(freeEnergyDf, cmap='plasma', annot=False)

    # Set x and y tick labels to 2 decimal places
    ax.set_xticklabels([f'{float(label):.2f}' for label in ax.get_xticks()], rotation=45)
    ax.set_yticklabels([f'{float(label):.2f}' for label in ax.get_yticks()], rotation=0)

    plt.title('Heatmap of Free Energy Data')
    plt.xlabel('Collective Variable 1')
    plt.ylabel('Collective Variable 2')
    plt.savefig(p.join(outDir, "freeEnergy.png"))

if __name__ == "__main__":
    main()