

def Forcefield(Protein: str, Water: str ) -> dict:
    """
        specify which force fields the simulation should use based on user input

        Args:
        Protein:
        Water:

    Returns:
        Dict {
            Protein force field
            Water force field
            Ion1 force field
            Ion234 force field
        }
    
    If users want to add additional force fields, use the following template

    if Protein == "Protein forcefield name"
        ForcefieldDict= {"protein":"protein leaprc file"}
            if Water == "Water forcefield name"
            ForcefieldDict.update({"Water": "water leaprc file})
            ForcefieldDict.update({"ion1": "ion1 frcmod file"})
            ForcefieldDict.update({"ion2": "ion234 frcmod"})
    """
    if Protein == "a99SB-disp":
        ForcefieldDict= {"protein":"leaprc.protein.ff99SBdisp"}
        if Water == "Tip4p":
            ForcefieldDict.update({"water": "leaprc.water.tip4pd-a99SBdisp"})
            ForcefieldDict.update({"ion1": "frcmod.ions1lm_126_tip4pew"})
            ForcefieldDict.update({"ion2": "frcmod.ions234lm_126_tip4pew"})
            ForcefieldDict.update({"box": "TIP4PBOX"})
        elif Water == "Tip3p":
            ForcefieldDict.update({"water": "leaprc.water.tip3p"})
            ForcefieldDict.update({"ion1": "frcmod.ions1lm_126_tip3p"})
            ForcefieldDict.update({"ion2": "frcmod.ions234lm_126_tip3p"})
            ForcefieldDict.update({"box": "TIP3PBOX"})
        elif Water == "OPC":
            ForcefieldDict.update({"water": "leaprc.water.opc"})
            ForcefieldDict.update({"ion1": "frcmod.ions1lm_126_spce"})
            ForcefieldDict.update({"ion2": "frcmod.ions234lm_126_spce"})
            ForcefieldDict.update({"box": "OPCBOX"})
    
    if Protein == "a99SB-ILDN":
        ForcefieldDict= {"protein": "oldff/leaprc.ff99SBildn"}
        if Water == "Tip4p":
            ForcefieldDict.update({"water": "leaprc.water.tip4pew"})
            ForcefieldDict.update({"ion1": "frcmod.ions1lm_126_tip4pew"})
            ForcefieldDict.update({"ion2": "frcmod.ions234lm_126_tip4pew"})
            ForcefieldDict.update({"box": "TIP4PEWBOX"})
        elif Water == "Tip3p":
            ForcefieldDict.update({"water": "leaprc.water.tip3p"})
            ForcefieldDict.update({"ion1": "frcmod.ions1lm_126_tip3p"})
            ForcefieldDict.update({"ion2": "frcmod.ions234lm_126_tip3p"})
            ForcefieldDict.update({"box": "TIP3PBOX"})
        elif Water == "OPC":
            ForcefieldDict.update({"water": "leaprc.water.opc"})
            ForcefieldDict.update({"ion1": "frcmod.ions1lm_126_spce"})
            ForcefieldDict.update({"ion2": "frcmod.ions234lm_126_spce"})
            ForcefieldDict.update({"box": "OPCBOX"})
    
    if Protein == "a19SB":
        ForcefieldDict= {"protein":"leaprc.protein.ff19SB"}
        if Water == "Tip4p":
            ForcefieldDict.update({"water": "leaprc.water.tip4pew"})
            ForcefieldDict.update({"ion1": "frcmod.ions1lm_126_tip4pew"})
            ForcefieldDict.update({"ion2": "frcmod.ions234lm_126_tip4pew"})
            ForcefieldDict.update({"box": "TIP4PEWBOX"})
        elif Water == "Tip3p":
            ForcefieldDict.update({"water": "leaprc.water.tip3p"})
            ForcefieldDict.update({"ion1": "frcmod.ions1lm_126_tip3p"})
            ForcefieldDict.update({"ion2": "frcmod.ions234lm_126_tip3p"})
            ForcefieldDict.update({"box": "TIP3PBOX"})
        elif Water == "OPC":
            ForcefieldDict.update({"water": "leaprc.water.opc"})
            ForcefieldDict.update({"ion1": "frcmod.ions1lm_126_spce"})
            ForcefieldDict.update({"ion2": "frcmod.ions234lm_126_spce"})
            ForcefieldDict.update({"box": "OPCBOX"})

    return ForcefieldDict
    
