import os
import sys
from os import path as p
import numpy as np
import re
import math
import random
import subprocess
from subprocess import run



def DropletFormation(PDBfile, Mode, charge):
    with open(PDBfile, 'r') as pdb:
        lines= pdb.readlines()
        AtomNumber=0
        for row in lines:
            if row.find("ATOM") != -1:
                AtomNumber +=1

        Coords= np.zeros((AtomNumber,3))
        
        for row in range(AtomNumber):
            if (lines[row]).find("ATOM") != -1:
                AtomVals= (lines[row]).split()
                Coords[row,:]= [AtomVals[6], AtomVals[7], AtomVals[8]]
    DropletCenter, DropletRadius= DefineDroplet(Coords)

    RLim=round((8*math.pi)*math.sqrt(0.072*8.8541878188*(10**-12)*((DropletRadius*(10**-9))**3))/(1.602*10**(-19)))
    print(f"\n-->    Created a droplet centered around {round(DropletCenter[0],2)}, {round(DropletCenter[1],2)}, {round(DropletCenter[2],2)}, of radius {DropletRadius} nm")
    if Mode == "Positive":
        Hydronium_Number= RLim
        Hydroxide_Number= RLim- charge - round(RLim*0.1)
    elif Mode == "Negative":
        Hydronium_Number= RLim+ charge + round(RLim*0.1)
        Hydroxide_Number= RLim

    ESI_Path= (os.path.dirname(os.path.abspath(sys.argv[0])))+ "/ESIWing/SpecialistEquipment"
    

    print(f"\n-->    Adding {Hydronium_Number} Sodium ions and {Hydroxide_Number} Chloride ions to the droplet")

    TleapInputs: str = p.join( "Droplet.in")
    with open(TleapInputs, "w") as f:
        f.write(f"source leaprc.protein.ff19SB \n")
        f.write(f"source leaprc.water.tip3p \n")
        f.write(f"loadamberparams frcmod.ions1lm_126_tip3p\n")
        f.write(f"loadamberparams frcmod.ions234lm_126_tip3p\n")
        #f.write(f"loadamberparams {ESI_Path}/Hydroxide.frcmod \n")
        #f.write(f"loadamberprep {ESI_Path}/Hydroxide.prep \n")
        #f.write(f"loadamberparams {ESI_Path}/Hydronium.frcmod \n")
        #f.write(f"loadamberprep {ESI_Path}/Hydronium.prep \n")
        f.write(f"mol = loadpdb {PDBfile} \n")
        f.write(f"solvateShell mol TIP3PBOX 55 2.0 \n")
        f.write("solvateCap mol TIP3PBOX {0.00, 0.00, 0.00} " + str(DropletRadius*10)+" 1.0 \n")
        #f.write(f"addions mol +O3 {Hydronium_Number}\n")
        #f.write(f"addions mol OH- {Hydroxide_Number}\n")
        f.write(f"addions mol Na+ {Hydronium_Number}\n")
        f.write(f"addions mol Cl- {Hydroxide_Number}\n")
        f.write(f"savepdb mol {PDBfile} \n")
        f.write(f"quit \n")
    
    print("\n-->    preparing the droplet")
    TleapCommand = "tleap -f Droplet.in"
    try: 
       # Execute the command and capture its output
        result: subprocess.CompletedProcess[str] = subprocess.run(
            TleapCommand.split(),
            capture_output=True,
            check=True,
            text=True, 
            env = os.environ
            )
    except Exception as errorMessage:
            print(errorMessage)
        
    return
                


def DefineDroplet(Coordinates):
    Center= np.mean(Coordinates, axis=0)
    Vectors= np.zeros((np.size(Coordinates,0),1))
    for I in range(np.size(Coordinates,0)):
        Vectors[I]= math.sqrt((Coordinates[I,0]-Center[0])*(Coordinates[I,0]-Center[0])+(Coordinates[I,1]-Center[1])*(Coordinates[I,1]-Center[1])+(Coordinates[I,2]-Center[2])*(Coordinates[I,2]-Center[2]))
    Radius= round((np.max(Vectors)+55)/10)
    return Center, Radius

