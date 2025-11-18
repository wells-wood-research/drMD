from subprocess import call
from pdbUtils import pdbUtils


"""
1. rename desired residues to AIB
2. delete sidechains
"""

pdb4amber = "/home/esp/anaconda3/envs/drMD/bin/pdb4amber"

inputPdb = "NonaPeptide_3AIB_input.pdb"

pdbDf = pdbUtils.pdb2df(inputPdb)
## rename resids
aibLocations = [3,6,9]
pdbDf.loc[pdbDf["RES_ID"].isin(aibLocations), "RES_NAME"] = "AIB"

# Remove sidechain atoms for AIB residues, keeping only backbone atoms
aibSideChainIndexes = pdbDf[(pdbDf["RES_NAME"] == "AIB") & (~pdbDf["ATOM_NAME"].isin(["C", "O", "CA", "HA", "N", "H"]))].index
pdbDf = pdbDf[~pdbDf.index.isin(aibSideChainIndexes)]
ch3Indexes = pdbDf[(pdbDf["RES_NAME"] == "NME") & (pdbDf["ATOM_NAME"] == "CH3")].index
pdbDf = pdbDf[~pdbDf.index.isin(ch3Indexes)]


# Save the modified DataFrame to a PDB file
renamedPdb = "NonaPeptide_3AIB_renamed.pdb"
pdbUtils.df2pdb(pdbDf, renamedPdb)

## add missing and protonate with tleap
ncaaMol2 = "/home/esp/scriptDevelopment/drMD/ExampleInputs/Worked_Example_7_Non_Canonical_Amino_Acids/AIB.mol2"
ncaaFrcmod = "/home/esp/scriptDevelopment/drMD/ExampleInputs/Worked_Example_7_Non_Canonical_Amino_Acids/AIB.frcmod"
ncaaLib = "/home/esp/scriptDevelopment/drMD/ExampleInputs/Worked_Example_7_Non_Canonical_Amino_Acids/AIB.lib"

tleapIn = "tleap.in"
with open(tleapIn, "w") as f:
    f.write("source leaprc.protein.ff19SB\n")
    f.write("source leaprc.gaff2\n")
    f.write("source leaprc.water.tip3p\n")
    f.write("loadamberparams frcmod.ions1lm_126_tip3p\n")
    f.write("loadamberparams frcmod.ions234lm_126_tip3p\n")
    f.write(f"ncaa = loadmol2 {ncaaMol2}\n")
    f.write(f"loadamberparams {ncaaFrcmod}\n")
    f.write(f"loadoff {ncaaLib}\n")
    f.write(f"mol = loadpdb {renamedPdb}\n")
    f.write("savepdb mol NonaPeptide_3AIB_protonated.pdb\n")
    f.write("quit\n")
tleapOut = "tleap.out"
call(["tleap", "-f", tleapIn], stdout=open(tleapOut, "w"))

pdbDf = pdbUtils.pdb2df("NonaPeptide_3AIB_protonated.pdb")
pdbUtils.df2pdb(pdbDf, "NonaPeptide_3AIB.pdb")


## rename MOL to AIB
