import os
import sys
from os import path as p
import numpy as np
import re
import math
import random
import subprocess
from itertools import combinations
from subprocess import run
from pdbUtils import pdbUtils
from UtilitiesCloset import drListInitiator
import pandas as pd

def DynamicCharger (PDBfile: str, TargetCharge: int ):
    
    print("\n-->    Adding charges to protein to reach target charge of: ", TargetCharge)
    def BuildResidueOrder(pdbDf: pd.DataFrame):
        residueOrder = []
        seenResidues = set()
        for _, atom in pdbDf.iterrows():
            residueKey = (atom["CHAIN_ID"], atom["RES_ID"])
            if residueKey in seenResidues:
                continue
            seenResidues.add(residueKey)
            residueOrder.append((atom["CHAIN_ID"], int(atom["RES_ID"]), atom["RES_NAME"]))

        renumberedDf = pdbDf.copy()
        for newResidueId, (chainId, originalResidueId, _) in enumerate(residueOrder, start=1):
            renumberedDf.loc[
                (renumberedDf["CHAIN_ID"] == chainId) & (renumberedDf["RES_ID"] == originalResidueId),
                "RES_ID"
            ] = newResidueId
        return renumberedDf, residueOrder

    def BuildChainTermini(seqChains: list):
        firstResidueByChain = {}
        lastResidueByChain = {}
        for residueIndex, chainId in enumerate(seqChains):
            firstResidueByChain.setdefault(chainId, residueIndex)
            lastResidueByChain[chainId] = residueIndex
        return set(firstResidueByChain.values()), set(lastResidueByChain.values())
    def GetProteinCharge(Seq:list, SeqChains: list):
        ResNumber= len(Seq)
        ResCharge=np.zeros((ResNumber))
        PosResCodes= "LYS", "ARG", "HIP"
        NegResCodes= "ASP", "GLU", "CYM"
        ## gets chareg for each residue in the sequence
        for Res in range(ResNumber):
            if Seq[Res] in PosResCodes:
                ResCharge[Res] = 1.0
            elif Seq[Res] in NegResCodes:
                ResCharge[Res] = -1.0
            else:
                ResCharge[Res] = 0.0 


        ## Checks capping groups on each chain rather than only the full complex
        ChainStarts, ChainEnds = BuildChainTermini(SeqChains)
        for ChainStart in ChainStarts:
            if Seq[ChainStart] != "ACE":
                ResCharge[ChainStart] = ResCharge[ChainStart] +1
        for ChainEnd in ChainEnds:
            if Seq[ChainEnd] != "NME":
                ResCharge[ChainEnd] = ResCharge[ChainEnd] -1   
        TotalCharge=int(sum(ResCharge))
        return TotalCharge, ResCharge
    
    def GetTitrateSites(pdbDf: pd.DataFrame):
        ## gets the alpha carbon for each titrateable site
        TitrateSites= pd.DataFrame(pdbDf.iloc[0:0])
        TitratableResidues= "LYS", "LYN", "HIS", "HIE", "HID", "HIP", "ASP", "GLU", "ASH", "GLH", "ARG", "CYS", "CYM"
        for AtmNumb, Atom in pdbDf.iterrows():
            if  Atom["ATOM_NAME"] == "CA":
                if Atom["RES_NAME"] in TitratableResidues:
                    new_row = pd.DataFrame([Atom.copy()])
                    TitrateSites= TitrateSites._append(Atom.copy(), ignore_index=True)
            elif Atom["ATOM_NAME"] == "OXT" or Atom["RES_NAME"] == "NME":
                new_row = pd.DataFrame([Atom.copy()])
                TitrateSites= TitrateSites._append(Atom.copy(), ignore_index=True)
        return TitrateSites
    
    def AddCharge(TitrateSites: pd.DataFrame, CurrentCharge: int, TargetCharge: int, ResCharge, Seq: list):

        ## finds the number of additional charges needed
        ChargeToAdd= int(TargetCharge-CurrentCharge)
        ProtResCodes= "LYS", "HIP", "ASH", "GLH", "CYS"
        DeprotResCodes= "ASP", "GLU", "LYN", "HID", "HIS", "HIE", "CYM"

        ## Adds random charges to the protein
        for AddedCharge in range(abs(ChargeToAdd)):
            SitesToChange=[]
            for SiteIndex, Site in TitrateSites.iterrows():
                if ChargeToAdd >= 1:
                    if Site.loc["RES_NAME"] in DeprotResCodes:
                        SitesToChange.append(Site.copy())
                elif ChargeToAdd < 0:
                    if Site.loc["RES_NAME"] in ProtResCodes:
                        SitesToChange.append(Site.copy())

            if not SitesToChange:
                break

            if ChargeToAdd >= 1:
                ## adds capping group if positive
                #if AddedCharge >= 1 and ResCharge[-1] == -1:
                #    TargetSite= TitrateSites.iloc[-1].copy()
                #    TitrateSites.loc[len(TitrateSites)-1, "ATOM_NAME"] = "N"
                #    TitrateSites.loc[len(TitrateSites)-1, "RES_NAME"] = "NME"
                #    TitrateSites.loc[len(TitrateSites)-1, "RES_ID"] = len(Seq) + 1
                #    ResCharge[-1]= 0
                #    Seq.append("NME")
                #    TitrateSites.iloc[[-1]]= TargetSite.copy()
                ## picks a random site to protonate
                #else:
                    TargetSite=random.choice(SitesToChange).copy()
                    if TargetSite.loc["RES_NAME"] == "LYN":
                        TitrateSites.loc[TargetSite.name, "RES_NAME"] = "LYS"
                        ResCharge[TargetSite["RES_ID"]-1] = 1.0
                        Seq[TargetSite["RES_ID"]-1] = "LYS"
                    elif TargetSite.loc["RES_NAME"] == "HID" or TargetSite.loc["RES_NAME"] == "HIS" or TargetSite.loc["RES_NAME"] == "HIE":
                        TitrateSites.loc[TargetSite.name, "RES_NAME"] = "HIP"
                        ResCharge[TargetSite["RES_ID"]-1] = 1.0
                        Seq[TargetSite["RES_ID"]-1] = "HIP"
                    elif TargetSite.loc["RES_NAME"] == "ASP":
                        TitrateSites.loc[TargetSite.name, "RES_NAME"] = "ASH"
                        ResCharge[TargetSite["RES_ID"]-1] = 0
                        Seq[TargetSite["RES_ID"]-1] = "ASH"
                    elif TargetSite.loc["RES_NAME"] == "GLU":
                        TitrateSites.loc[TargetSite.name, "RES_NAME"] = "GLH"
                        ResCharge[TargetSite["RES_ID"]-1] = 0
                        Seq[TargetSite["RES_ID"]-1] = "GLH"
                    elif TargetSite.loc["RES_NAME"] == "CYM":
                        TitrateSites.loc[TargetSite.name, "RES_NAME"] = "CYS"
                        ResCharge[TargetSite["RES_ID"]-1] = 0
                        Seq[TargetSite["RES_ID"]-1] = "CYS"
            ## same but for negative charges (negative mode)        
            elif ChargeToAdd <= -1:
                #if AddedCharge == 0 and ResCharge[-1] == 1:
                #    TargetSite= TitrateSites.iloc[0].copy()
                #    TargetSite["ATOM_NAME"] = "C"
                #    TargetSite["RES_NAME"] = "ACE"
                #    TargetSite["RES_ID"] = 0
                #    ResCharge[0]= 0
                #    Seq.append("ACE")
                #    TitrateSites.iloc[0]= TargetSite
                #else:
                    TargetSite=random.choice(SitesToChange).copy()
                    if TargetSite.loc["RES_NAME"] == "LYS":
                        TitrateSites.loc[TargetSite.name, "RES_NAME"] = "LYN"
                        ResCharge[TargetSite["RES_ID"]-1] = 0
                        Seq[TargetSite["RES_ID"]-1] = "LYN"
                    elif TargetSite.loc["RES_NAME"] == ("HIP"):
                        TitrateSites.loc[TargetSite.name, "RES_NAME"] = "HIS"
                        ResCharge[TargetSite["RES_ID"]-1] = 0
                        Seq[TargetSite["RES_ID"]-1] = "HIS"
                    elif TargetSite.loc["RES_NAME"] == "ASH":
                        TitrateSites.loc[TargetSite.name, "RES_NAME"] = "ASP"
                        ResCharge[TargetSite["RES_ID"]-1] = -1
                        Seq[TargetSite["RES_ID"]-1] = "ASP"
                    elif TargetSite.loc["RES_NAME"] == "GLH":
                        TitrateSites.loc[TargetSite.name, "RES_NAME"] = "GLU"
                        ResCharge[TargetSite["RES_ID"]-1] = -1
                        Seq[TargetSite["RES_ID"]-1] = "GLU"
                    elif TargetSite.loc["RES_NAME"] == "CYS":
                        TitrateSites.loc[TargetSite.name, "RES_NAME"] = "CYM"
                        ResCharge[TargetSite["RES_ID"]-1] = -1
                        Seq[TargetSite["RES_ID"]-1] = "CYM"
        return TitrateSites, Seq

    def UpdateResidues(pdbDf, Seq):
                ## Handels capping groupds
        
        
        
        ## gets amino acid residue names from the list initiator
        amino_acid_residues = drListInitiator.get_amino_acid_residue_names()   
        for AtmNumb, Atom in pdbDf.iterrows():
            
            ## Changes RES_ID to the same as the sequence and removes any hydrogens
            if Atom["RES_NAME"] in amino_acid_residues and Atom["RES_ID"] <= len(Seq):
                if Atom["RES_NAME"] != Seq[Atom["RES_ID"]-1]:
                    pdbDf.loc[AtmNumb, "RES_NAME"] = Seq[Atom["RES_ID"]-1]
                if Atom["ATOM_NAME"].find("H") != -1:
                    pdbDf= pdbDf[pdbDf.index!= AtmNumb]

        pdbDf= pdbDf.sort_index().reset_index(drop=True)
        outDf= pdbDf.copy()
        uniqueResidues = []
        seen_residues = set() 
        for index, row in pdbDf.iterrows():
            res_tuple = (row["CHAIN_ID"], row["RES_ID"])
            if res_tuple not in seen_residues:
                uniqueResidues.append(res_tuple)
                seen_residues.add(res_tuple)

        new_res_id_counter = 0
        for chain_id, original_res_id in uniqueResidues:
            new_res_id_counter += 1
            outDf.loc[(pdbDf["CHAIN_ID"] == chain_id) & (pdbDf["RES_ID"] == original_res_id), "RES_ID"] = new_res_id_counter
        pdbDf= outDf.copy()
        return pdbDf, Seq
            
    def CoulombCalc(ResCharge, TitrateSites):       
        CoulombMatrix= np.zeros((len(TitrateSites), len(TitrateSites)))
        for SiteIndexI, SiteI in TitrateSites.iterrows():
            Ires= int(SiteI["RES_ID"])-1
            ChargeI= ResCharge[Ires]
            CoordsI= float(SiteI["X"]), float(SiteI["Y"]), float(SiteI["Z"])
            for SiteIndexJ, SiteJ in TitrateSites.iterrows():
                if SiteIndexI > SiteIndexJ and SiteIndexI != SiteIndexJ:
                    Jres= int(SiteJ["RES_ID"])-1
                    ChargeJ= ResCharge[Jres]
                    CoordsJ= float(SiteJ["X"]), float(SiteJ["Y"]), float(SiteJ["Z"])
                    VectorIJ= [CoordsI[0]-CoordsJ[0], CoordsI[1]-CoordsJ[1], CoordsI[2]-CoordsJ[2]]
                    DistIJ= math.sqrt(VectorIJ[0]*VectorIJ[0]+VectorIJ[1]*VectorIJ[1]+VectorIJ[2]*VectorIJ[2])
                    CoulombMatrix[SiteIndexI, SiteIndexJ]= (ChargeI*ChargeJ)*138.9/(DistIJ/10)
        SiteEnergies= np.zeros((len(TitrateSites)))
        for SiteIndex, Site in TitrateSites.iterrows():
            for I in range(len(TitrateSites)):
                if I == SiteIndex:
                    for J in range(len(TitrateSites)):
                        SiteEnergies[SiteIndex] = SiteEnergies[SiteIndex]+ CoulombMatrix[I,J]
                else:
                    for J in range(len(TitrateSites)):
                        if J== SiteIndex:
                            SiteEnergies[SiteIndex] = SiteEnergies[SiteIndex]+ CoulombMatrix[I,J]  
        
        TotalCoulombInt= np.sum(CoulombMatrix)
        return TotalCoulombInt, SiteEnergies
    
    def PACalc(TitrateSites, SeqChains: list):
        PAAffinitites= 886.6, 918.0, 1002.0, 1452.7, 1453.5, 952.7, 1450.0 #N-term, LYS, ARG, ASP, GLU, HIS, C-term
        PAEnergy= 0
        SitePA= np.zeros(len(TitrateSites))
        ChainStarts, ChainEnds = BuildChainTermini(SeqChains)
        I=0
        for SiteIndex, Site in TitrateSites.iterrows():
            SiteResidueIndex = int(Site["RES_ID"]) - 1
            if SiteResidueIndex in ChainStarts and Site["RES_NAME"] != "ACE":
                PAEnergy -= PAAffinitites[0]
                SitePA[I]= -PAAffinitites[0]
            elif Site["RES_NAME"] == "LYS":
               PAEnergy -= PAAffinitites[1]
               SitePA[I]= -PAAffinitites[1]
            elif Site["RES_NAME"] == "ARG":
               PAEnergy -= PAAffinitites[2]
               SitePA[I]= -PAAffinitites[2]
            elif Site["RES_NAME"] == "ASH":
               PAEnergy -= PAAffinitites[3]
               SitePA[I]= -PAAffinitites[3]
            elif Site["RES_NAME"] == "GLH":
               PAEnergy -= PAAffinitites[4]
               SitePA[I]= -PAAffinitites[4]
            elif Site["RES_NAME"] == "HIP":
                PAEnergy -= PAAffinitites[5]
                SitePA[I]= -PAAffinitites[5]
            elif SiteResidueIndex in ChainEnds and Site["RES_NAME"] == "NME":
                PAEnergy -= PAAffinitites[6]
                SitePA[I]= -PAAffinitites[6]
            I += 1
        return PAEnergy, SitePA
                

    def MoveCharge(TitrateSites, SiteCoulomb, OldResCharge, InitialEnergy, SeqChains: list):
        IndexArr=np.where(SiteCoulomb == SiteCoulomb.max())   
        Index= IndexArr[0]
        OldSite= TitrateSites.iloc[(Index)].copy()
        ProtResCodes= "LYS", "HIP", "ASH", "GLH", "CYS"
        DeprotResCodes= "ASP", "GLU", "LYN", "HID", "HIS", "HIE", "CYM"
        PosResCodes= "LYS", "HIP"
        NegResCodes= "ASP", "GLU", "CYM"
        protState= "NEUTRAL"
        ## find if the old site is positive, negative or neutral and what the new charge will be
        if str(OldSite.iloc[0]["RES_NAME"]) in ProtResCodes:
            protState= "PROT"
        elif str(OldSite.iloc[0]["RES_NAME"]) in DeprotResCodes:
            protState= "DEPROT"
        if str(OldSite.iloc[0]["RES_NAME"]) in PosResCodes:
            Charge = 1
            DCharge= 0
        elif str(OldSite.iloc[0]["RES_NAME"]) in NegResCodes:
            Charge = -1
            DCharge= 0
        else:
            Charge = 0
            if str(OldSite.iloc[0]["RES_NAME"]) in ("ASH", "GLH"):
                DCharge= -1
            elif str(OldSite.iloc[0]["RES_NAME"]) in ("LYN", "HIS", "HID", "HIE"):
                DCharge= 1
            else:
                DCharge= 0
        
        ##find the new residue name
        if str(OldSite.iloc[0]["RES_NAME"]) == "LYS":
            ChangedRes= "LYN"
        elif str(OldSite.iloc[0]["RES_NAME"]) == "LYN":
            ChangedRes= "LYS"
        elif str(OldSite.iloc[0]["RES_NAME"]) == "HIS" or str(OldSite.iloc[0]["RES_NAME"]) == "HID" or str(OldSite.iloc[0]["RES_NAME"]) == "HIE":
            ChangedRes= "HIP"
        elif str(OldSite.iloc[0]["RES_NAME"]) == "HIP":
            ChangedRes= "HIS"
        elif str(OldSite.iloc[0]["RES_NAME"]) == "GLH":
            ChangedRes= "GLU"
        elif str(OldSite.iloc[0]["RES_NAME"]) == "GLU":
            ChangedRes= "GLH"
        elif str(OldSite.iloc[0]["RES_NAME"]) == "ASP":
            ChangedRes= "ASH"
        elif str(OldSite.iloc[0]["RES_NAME"]) == "ASH":
            ChangedRes= "ASP"
        elif str(OldSite.iloc[0]["RES_NAME"]) == "CYM":
            ChangedRes= "CYS"
        elif str(OldSite.iloc[0]["RES_NAME"]) == "CYS":
            ChangedRes= "CYM"

        

        
        ## find all sites within 15 angstroms that can exchange protons with the old site
        CandidateSites = pd.DataFrame(columns=TitrateSites.columns)
        OldCoords= [OldSite.iloc[0]["X"], OldSite.iloc[0]["Y"], OldSite.iloc[0]["Z"]]
        for SiteIndex, Site in TitrateSites.iterrows():
            if Site.loc["RES_ID"] != OldSite.iloc[0]["RES_ID"]:
                SiteCoords= [Site.loc["X"], Site.loc["Y"], Site.loc["Z"]]
                SiteDistance= math.sqrt((float(SiteCoords[0])-float(OldCoords[0]))**2+(float(SiteCoords[1])-float(OldCoords[1]))**2+ ((float(SiteCoords[2]))-float(OldCoords[2]))**2)
                if SiteDistance <= 15.0:
                    if protState == "PROT":
                        if Site.loc["RES_NAME"] in DeprotResCodes:
                            CandidateSites= CandidateSites._append(Site, ignore_index=True)
                    elif protState == "DEPROT":
                        if Site.loc["RES_NAME"] in ProtResCodes:
                            CandidateSites= CandidateSites._append(Site, ignore_index=True)
        if CandidateSites.empty:
            return TitrateSites, OldResCharge
        ## finds the most optimal site to move to

        ## Initialises the initial energy and sites
        OptEnergy= InitialEnergy
        OptSites= TitrateSites.copy()
        OptResCharge= np.copy(OldResCharge)
        InitialRes= np.copy(OldResCharge) 
        InitialRes[int(OldSite.iloc[0]["RES_ID"])-1] = DCharge 
        InitialSites= TitrateSites.copy()
        InitialSites.loc[int(Index[0]), "RES_NAME"] = ChangedRes
        ## Finds most optimal site to exchange with
        TestRCharge= np.copy(InitialRes)
        TestSites= InitialSites.copy()
        for CandidaSiteIndex, Sites in CandidateSites.iterrows():
            ## finds the index of the residue in the charge array and the index of the site in the titrateable sites dataframe                 
            ResIndex= int(Sites.loc["RES_ID"])-1
            for I, Tsite in TitrateSites.iterrows():
                if Tsite.loc["RES_ID"]== Sites.loc["RES_ID"]:
                    SiteIndex= I
            ## Changes the residue name and updates the charge array for testing
            if Sites.loc["RES_NAME"] == "LYS":
                TestSites.loc[SiteIndex, "RES_NAME"] = "LYN"
                TestRCharge[ResIndex] -= 1.0 
            elif Sites.loc["RES_NAME"] == "LYN":
                TestSites.loc[SiteIndex, "RES_NAME"] = "LYS"
                TestRCharge[ResIndex] += 1.0
            elif Sites.loc["RES_NAME"] == "HIP":
                TestSites.loc[SiteIndex, "RES_NAME"] = "HIS"
                TestRCharge[ResIndex]  -= 1.0
            elif Sites.loc["RES_NAME"] == "HIS" or Sites.loc["RES_NAME"] =="HID" or Sites.loc["RES_NAME"] =="HIE":
                TestSites.loc[SiteIndex, "RES_NAME"] = "HIP"
                TestRCharge[ResIndex]  += 1.0
            elif Sites.loc["RES_NAME"] == "GLU":
                TestSites.loc[SiteIndex, "RES_NAME"] = "GLH"
                TestRCharge[ResIndex]  += 1.0
            elif Sites.loc["RES_NAME"] == "GLH":
                TestSites.loc[SiteIndex, "RES_NAME"] = "GLU"
                TestRCharge[ResIndex]  -= 1.0
            elif Sites.loc["RES_NAME"] == "ASP":
                TestSites.loc[SiteIndex, "RES_NAME"] = "ASH"
                TestRCharge[ResIndex]  += 1.0
            elif Sites.loc["RES_NAME"] == "ASH":
                TestSites.loc[SiteIndex, "RES_NAME"] = "ASP"
                TestRCharge[ResIndex]  -= 1.0
            elif Sites.loc["RES_NAME"] == "CYM":
                TestSites.loc[SiteIndex, "RES_NAME"] = "CYS"
                TestRCharge[ResIndex] += 1.0
            elif Sites.loc["RES_NAME"] == "CYS":
                TestSites.loc[SiteIndex, "RES_NAME"] = "CYM"
                TestRCharge[ResIndex] -= 1.0
            TestCoulombEnergy, TestSiteCoulomb= CoulombCalc(TestRCharge, TestSites)
            TestPAEnergy, TestsitePAEnergy= PACalc(TestSites, SeqChains)
            TotalEnergy= TestCoulombEnergy+ TestPAEnergy
            if float(TotalEnergy) < float(OptEnergy):
                OptEnergy= TotalEnergy
                OptSites= TestSites.copy()
                OptResCharge= np.copy(TestRCharge) 
            np.copyto(TestRCharge, InitialRes)
            TestSites= InitialSites.copy()
        NewTitrateSites= OptSites.copy()
        NewResSites= np.copy(OptResCharge)
        return NewTitrateSites, NewResSites 




    ## Initilaise the pdb DataFrame
    pdbDf = pdbUtils.pdb2df(PDBfile)
    pdbDf, ResidueOrder = BuildResidueOrder(pdbDf)
    NewpdbDf= pdbDf.copy()

    ## Get protein Sequence
    amino_acid_residues = drListInitiator.get_amino_acid_residue_names()
    Seq=[]
    SeqChains=[]
    for ChainId, _, ResidueName in ResidueOrder:
        Seq.append(ResidueName)
        SeqChains.append(ChainId)
    ## Get protein Charge and titrateable sites
    TotalCharge, ResCharge= GetProteinCharge(Seq, SeqChains)
    if TargetCharge == 0:
        TargetCharge= TotalCharge
    TitrateSites= GetTitrateSites(pdbDf)

    ## if the charge of the pdb is not the same as the target charge
    if TargetCharge != TotalCharge:
        TitrateSites, Seq =AddCharge(TitrateSites, TotalCharge, TargetCharge, ResCharge, Seq)
        NewpdbDf, Seq= UpdateResidues(NewpdbDf, Seq)
        TotalCharge, ResCharge= GetProteinCharge(Seq, SeqChains)
    ## Get the Coulombic and proton affinity energy for the current sequence
    CoulombEnergy, SiteCoulomb= CoulombCalc(ResCharge, TitrateSites)
    PAEnergy, sitePAEnergy= PACalc(TitrateSites, SeqChains)
    TotalEnergy= CoulombEnergy+ PAEnergy

    ## optimise proton location to reduce total energy
    OptEnergy= TotalEnergy
    TestTotalEnergy= TotalEnergy
    TestResCharge= np.copy(ResCharge)
    TestSiteCoulomb= np.copy(SiteCoulomb)
    TestTitrateSites= TitrateSites.copy()
    
    for I in range(10):
        TestTitrateSites, TestResCharge= MoveCharge(TestTitrateSites, TestSiteCoulomb, TestResCharge, TestTotalEnergy, SeqChains)
        TestCoulombEnergy, TestSiteCoulomb= CoulombCalc(TestResCharge, TestTitrateSites)
        TestPAEnergy, TestsitePAEnergy= PACalc(TestTitrateSites, SeqChains)
        TestTotalEnergy= TestCoulombEnergy+ TestPAEnergy
        if TestTotalEnergy < OptEnergy:
            TitrateSites= TestTitrateSites.copy()
            OptEnergy= TestTotalEnergy
            SiteCoulomb= np.copy(TestSiteCoulomb)
            ResCharge= np.copy(TestResCharge)
        TestTitrateSites= TitrateSites.copy()
        TestResCharge= np.copy(ResCharge)
        TestTotalEnergy= OptEnergy
    

      

        
    ## updates the residues in the pdbDf
    for IndexSite, Site in TitrateSites.iterrows():
        ResNumber= int(Site.loc["RES_ID"])-1
        Seq[ResNumber] = Site.loc["RES_NAME"]

    ## Runs pdb2pqr to fix and clean up the charged pdb file
    NewpdbDf, Seq= UpdateResidues(NewpdbDf, Seq)
    TotalCharge, ResCharge= GetProteinCharge(Seq, SeqChains)
    outPdbFile = PDBfile.removesuffix('.pdb') + '_Charged.pdb'
    pdbUtils.df2pdb(NewpdbDf, outPdbFile)
    protPqr= PDBfile.removesuffix(".pdb")+f"_{TargetCharge}.pdb"
    pdb2pqrCommand: str = ["pdb2pqr",
                         "--ffout", "AMBER",
                          "--keep-chain",
                             outPdbFile, protPqr]
    try: 
        # Execute the command and capture its output
            result: subprocess.CompletedProcess[str] = subprocess.run(
            pdb2pqrCommand,
            capture_output=True,
            check=True,
            text=True, 
            env = os.environ
            )
    except Exception as errorMessage:
            print(errorMessage)
    os.remove(outPdbFile)
    os.remove(protPqr.removesuffix(".pdb")+".log")
    chargedPdbfile= protPqr
    return chargedPdbfile

def get_total_charge(protPqr: str) -> int:

    with open(protPqr, 'r') as f:
        lines = f.readlines()
        charge= 0
        for line in lines:
            if line.startswith("ATOM"):
                parts = line.split()
                if len(parts) >= 2:
                    charge += float(parts[-2])
    return int(round(charge))


def SimpleCharger( PDBfile: str, TargetCharge: int):
    print(f"\n-->    Using the Simple charge method to get to {TargetCharge} charge")

    ##removes hydrogens
    pdbDf = pdbUtils.pdb2df(PDBfile)
    for AtmNumb, Atom in pdbDf.iterrows():
        if Atom["ATOM_NAME"].find("H") != -1:
            pdbDf= pdbDf[pdbDf.index!= AtmNumb]
    strippedPdb= PDBfile.removesuffix("_RENAMED.pdb")+"stripped.pdb"
    pdbUtils.df2pdb(pdbDf, strippedPdb)
    chargedPqr=PDBfile.removesuffix("_RENAMED.pdb")+f"_initial.pqr"

    def run_pdb2pqr(command: list[str], outputPqr: str):
        try:
            subprocess.run(
                command,
                capture_output=True,
                check=True,
                text=True,
                env=os.environ,
            )
        except Exception as errorMessage:
            print(errorMessage)
            return None, None

        outputDf = pdbUtils.pdb2df(outputPqr)
        totalCharge = get_total_charge(outputPqr)
        return outputDf, totalCharge

    def cleanup_outputs(outputPqr: str):
        if os.path.exists(outputPqr):
            os.remove(outputPqr)
        outputLog = outputPqr.removesuffix(".pqr") + ".log"
        if os.path.exists(outputLog):
            os.remove(outputLog)

    ## initial pdb2pqr
    pdb2pqrCommand: str = ["pdb2pqr",
                         "--ffout", "AMBER",
                          "--keep-chain",
                             strippedPdb, chargedPqr]

    NewpdbDf, TotalCharge = run_pdb2pqr(pdb2pqrCommand, chargedPqr)
    if NewpdbDf is None:
        raise RuntimeError(f"pdb2pqr failed while charging {PDBfile}")
            
    if TargetCharge == 0:
        TargetCharge= TotalCharge
    if TotalCharge == TargetCharge:
        chargedPDB=PDBfile.removesuffix(".pdb")+f"_{TargetCharge}.pdb"
        pdbUtils.df2pdb(NewpdbDf, chargedPDB)
        cleanup_outputs(chargedPqr)
        return chargedPDB



    pH: float = 7.4
    UBpH= 14
    LBpH= 0
    
    while TotalCharge != TargetCharge:

        protPqr= PDBfile.removesuffix(".pdb")+f"_{pH}.pqr"
        pdb2pqrCommand: str = ["pdb2pqr",
                         "--ffout", "AMBER",
                         "--titration-state-method", "propka",
                          "--keep-chain",
                           "--with-ph", str(float(pH)),
                             strippedPdb, protPqr]

        NewpdbDf, TotalCharge = run_pdb2pqr(pdb2pqrCommand, protPqr)
        if NewpdbDf is None:
            cleanup_outputs(protPqr)
            raise RuntimeError(f"pdb2pqr failed while probing pH {pH} for {PDBfile}")
        if TotalCharge < TargetCharge:
            UBpH= pH
            pH= UBpH- 0.5*(UBpH-LBpH)
        elif TotalCharge > TargetCharge:
            LBpH= pH
            pH= LBpH+ 0.5*(UBpH-LBpH)
        
        elif TotalCharge == TargetCharge:
            chargedPDB=PDBfile.removesuffix(".pdb")+f"_{TargetCharge}.pdb"
            pdbUtils.df2pdb(NewpdbDf, chargedPDB)
            cleanup_outputs(protPqr)
            return chargedPDB

        if UBpH-LBpH <= 0.01:
            for probePH in np.linspace(LBpH, UBpH, num=11):
                probePqr= PDBfile.removesuffix(".pdb")+f"_{probePH}.pqr"
                probeCommand: str = ["pdb2pqr",
                                    "--ffout", "AMBER",
                                    "--titration-state-method", "propka",
                                    "--keep-chain",
                                    "--with-ph", str(float(probePH)),
                                    strippedPdb, probePqr]
                ProbeDf, ProbeCharge = run_pdb2pqr(probeCommand, probePqr)
                if ProbeDf is None:
                    cleanup_outputs(probePqr)
                    continue
                if ProbeCharge == TargetCharge:
                    chargedPDB=PDBfile.removesuffix(".pdb")+f"_{TargetCharge}.pdb"
                    pdbUtils.df2pdb(ProbeDf, chargedPDB)
                    cleanup_outputs(probePqr)
                    cleanup_outputs(protPqr)
                    return chargedPDB
                cleanup_outputs(probePqr)

            cleanup_outputs(protPqr)
            raise RuntimeError(
                f"Unable to reach target charge {TargetCharge} for {PDBfile}; "
                f"last observed charge was {TotalCharge} at pH {pH}"
            )

        cleanup_outputs(protPqr)
    raise RuntimeError(f"Unable to reach target charge {TargetCharge} for {PDBfile}")