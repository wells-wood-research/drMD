## BASIC PYTHON LIBRARIES
import xml.etree.ElementTree as ET
from os import path as p
import math

## OPENMM LIBRARIES
import openmm as openmm
import openmm.unit as unit
from openmm import app

## drMD LIBRARIES
from ExaminationRoom import drLogger
from UtilitiesCloset import drSelector

## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils

## CLEAN CODE
from typing import  Dict, List
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath

###########################################################################################
def restraints_handler(
        system: openmm.System,
        prmtop: app.AmberPrmtopFile,
        inpcrd: app.AmberInpcrdFile,
        sim: Dict,
        saveFile: FilePath,
        pdbFile: FilePath) -> openmm.System:
    """
    Handle restraints in the system.

    Args:
        system (openmm.System): The OpenMM system object.
        prmtop (AmberPrmtopFile): The path to the prmtop file.
        inpcrd (AmberInpcrdFile): The path to the inpcrd file.
        sim (dict): The simulation dictionary containing restraint information.
        saveFile (FilePath): The path to save the system.
        pdbFile (FilePath): The path to the pdb file.

    Returns:
        openmm.System: The system with restraints applied.
    """
    ## if we are loading from a XML file, we need to clear all restraints
    if saveFile:
        if p.splitext(saveFile)[1] == ".xml":
            clear_all_restraints(saveFile)
    ## check if there are any restraints specified in simulation config
    if "restraintInfo" in sim:
        ## load restraintInfo from simluation config
        restraintInfo: List[Dict] = sim["restraintInfo"]
        ## create a counter for naming restraint parameters
        kNumber: int = 0
        ## loop through restraints 
        ## create position, distance, angle, and torsion restraints
        for restraint in restraintInfo:
            selection: List = restraint["selection"]
            parameters: Dict = restraint["parameters"]
            ## add a position restraint
            if restraint["restraintType"] == "position":
                system: openmm.System = create_position_restraint(system, inpcrd, selection, parameters, kNumber, pdbFile)
            ## add a distance restraint
            elif restraint["restraintType"] == "distance":
                system: openmm.System = create_distance_restraint(system, selection, parameters, kNumber, pdbFile)
            ## add an angle restraint
            elif restraint["restraintType"] == "angle":
                system: openmm.System = create_angle_restraint(system, selection, parameters, kNumber, pdbFile)
            ## add a torsion restraint
            elif restraint["restraintType"] == "torsion":
                system: openmm.System = create_torsion_restraint(system, selection, parameters, kNumber, pdbFile)
            ## increment kNumber
            kNumber += 1

    else:
        ## if we have a checkpoint file, we are continuing a simulation
        ## just return the checkpoint file as-is
        if saveFile:
            if p.splitext(saveFile)[1] == ".chk":
                return system

    inspect_restraints(system)
    return system




def inspect_restraints(system):
    forces = system.getForces()

    # Print detailed information about each force
    for i, force in enumerate(forces):
        drLogger.log_info(f"Force {i}: {force.__class__.__name__}")
        
        if isinstance(force, openmm.CustomBondForce):
            for j in range(force.getNumBonds()):
                particleA, particleB, bondParameters = force.getBondParameters(j)

                drLogger.log_info(f"Bond {j}: particles ({particleA}, {particleB}), length {bondParameters[1] * 10} Å, force constant {bondParameters[0]}")
        
        elif isinstance(force, openmm.CustomAngleForce):
            for j in range(force.getNumAngles()):
                particleA, particleB, particleC, angleParameters = force.getAngleParameters(j)
                drLogger.log_info(f"Angle {j}: particles ({particleA}, {particleB}, {particleC}), angle {round(angleParameters[1] * 180 /math.pi, 2)} degrees, force constant {angleParameters[0]}")
        
        elif isinstance(force, openmm.CustomTorsionForce):
            for j in range(force.getNumTorsions()):
                particleA, particleB, particleC, particleD, torsionParameters = force.getTorsionParameters(j)   
                drLogger.log_info(f"Torsion {j}: particles ({particleA}, {particleB}, {particleC}, {particleD}), angle {round(torsionParameters[1] * 180 /math.pi, 2)} degrees, force constant {torsionParameters[0]}")      
###########################################################################################
def create_position_restraint(
    system: openmm.System,
    inpcrd: app.AmberInpcrdFile,
    selection: str,
    parameters: Dict,
    kNumber: int,
    pdbFile: FilePath) -> openmm.System:
    """
    Creates a position restraint for a given system.

    Parameters:
        system (openmm.System): The system to add the position restraint to.
        inpcrd (Inpcrd): The Inpcrd object containing the positions of the atoms.
        selection (str): The selection string specifying the atoms to be restrained.
        kNumber (int): The number used to identify the force constant parameter.
        pdbFile (str): The path to the PDB file.

    Returns:
        openmm.System: The system with the position restraint added.
    """
    ## create a position restraint object, add parameters
    positionRestraint = openmm.CustomExternalForce(f"k{str(kNumber)}*periodicdistance(x, y, z, x0, y0, z0)^2")

    kForceConstant: int =  parameters["k"]
    positionRestraint.addGlobalParameter(f"k{str(kNumber)}", kForceConstant* unit.kilojoules_per_mole / unit.nanometer**2)
    positionRestraint.addPerParticleParameter("x0")
    positionRestraint.addPerParticleParameter("y0")
    positionRestraint.addPerParticleParameter("z0")

    ## use selectio to get atom indexes, add them to the restraint object
    restraintAtomIndexes: List[int] = drSelector.get_atom_indexes(selection, pdbFile)
    for restraintAtomIndex in restraintAtomIndexes:
        positionRestraint.addParticle(restraintAtomIndex,
                                       inpcrd.getPositions()[restraintAtomIndex])
        
    ## add position restrain to system
    system.addForce(positionRestraint)
    return system

###########################################################################################
def create_distance_restraint(system: openmm.System, selection: list, parameters: Dict, kNumber: int, pdbFile: FilePath) -> openmm.System:
    """
    Creates a distance restraint between two atoms for a given system.

    Parameters:
        system (openmm.System): The system to add the position restraint to.
        inpcrd (Inpcrd): The Inpcrd object containing the positions of the atoms.
        selection (list): The selection string specifying the atoms to be restrained.
        kNumber (int): The number used to identify the force constant parameter.
        pdbFile (str): The path to the PDB file.

    Returns:
        openmm.System: The system with the distance restraint added.
    """

    # Create the distance restraint object
    distanceRestraint: openmm.CustomBondForce = openmm.CustomBondForce(f"0.5 * k{str(kNumber)} * (r - r0)^2")
    ## add per bond parameters, k for force constant, r0 for desired distance
    distanceRestraint.addPerBondParameter(f"k{str(kNumber)}")
    distanceRestraint.addPerBondParameter("r0")

    # Get the indices of the two atoms to be restrained
    restraintAtomIndexes: List[int] = drSelector.get_atom_indexes(selection, pdbFile)
    if len(restraintAtomIndexes) != 2:
        raise ValueError("Expected exactly two atom indices for a distance restraint.")

    # Get target distance from the parameters dictionary, and convert from angstroms to nanometers
    kForceConstant: float = parameters["k"]
    targetDistance_nm: float = parameters["r0"] * unit.angstroms 

    # Add the atom pair and the calculated target distance in nanometers to the bond restraint
    distanceRestraint.addBond(restraintAtomIndexes[0], restraintAtomIndexes[1], [kForceConstant, targetDistance_nm])

    ## add force to system
    system.addForce(distanceRestraint)

    return system
###########################################################################################
def create_angle_restraint(system: openmm.System, selection: list, parameters: Dict, kNumber: int, pdbFile: FilePath) -> openmm.System:
    """
    Creates an angle restraint between three atoms for a given system.

    Parameters:
        system (openmm.System): The system to add the angle restraint to.
        selection (str): The selection string specifying the atoms to be restrained.
        parameters (dict): The parameters dictionary containing the force constant and target angle.
        kNumber (int): The number used to identify the force constant parameter.
        pdbFile (str): The path to the PDB file.

    Returns:
        openmm.System: The system with the angle restraint added.
    """

    # Create the angle restraint object
    angleRestraint: openmm.CustomAngleForce = openmm.CustomAngleForce(f"0.5 * k{str(kNumber)} * (theta - theta0)^2")

    ## add per angle parameters, k for force constant, theta0 for desired angle
    angleRestraint.addPerAngleParameter(f"k{str(kNumber)}")
    angleRestraint.addPerAngleParameter("theta0")

    ## add force to system
    system.addForce(angleRestraint)

    # Get the indices of the three atoms to be restrained
    restraintAngleAtoms: list = drSelector.get_atom_indexes(selection, pdbFile)
    if len(restraintAngleAtoms) != 3:
        raise ValueError("Expected exactly three atom indices for an angle restraint.")

    # Get target angle from the parameters dictionary and convert from degrees to radians
    kForceConstant: float = parameters["k"] * unit.kilojoules_per_mole / unit.radians**2
    targetAngle_degrees: float = parameters["theta0"] * unit.degrees


    # Add the atom triplet and the calculated target angle in radians to the angle restraint
    angleRestraint.addAngle(restraintAngleAtoms[0],
                             restraintAngleAtoms[1],
                               restraintAngleAtoms[2],
                                 [kForceConstant, targetAngle_degrees])

    return system
###########################################################################################
def create_torsion_restraint(system: openmm.System, selection: list, parameters: dict, kNumber: int, pdbFile: FilePath) -> openmm.System:
    """
    Creates a torsion restraint between four atoms for a given system.

    Parameters:
        system (openmm.System): The system to add the torsion restraint to.
        selection (str): The selection string specifying the atoms to be restrained.
        parameters (dict): The parameters dictionary containing the force constant and target torsion angle.
        kNumber (int): The number used to identify the force constant parameter.
        pdbFile (str): The path to the PDB file.

    Returns:
        openmm.System: The system with the torsion restraint added.
    """
    ## create the torsion restraint object  
    torsionRestraint: openmm.CustomTorsionForce = openmm.CustomTorsionForce(f"0.5*k{str(kNumber)}*(1-cos(theta-theta0))")
    
    ## Add parameters: k for force constant, phi0 for desired torsion angle
    torsionRestraint.addPerTorsionParameter(f"k{str(kNumber)}")
    torsionRestraint.addPerTorsionParameter("theta0")
    
    ## Add force to system
    system.addForce(torsionRestraint)
    
    ## Get indices of the four atoms to be restrained
    restraintTorsionAtoms: list = drSelector.get_atom_indexes(selection, pdbFile)
    if len(restraintTorsionAtoms) != 4:
        raise ValueError("Expected exactly four atom indices for a torsion restraint.")
    
    ## Get target torsion angle from parameters (in radians)
    kForceConstant: float = parameters.get("k", 1000) * unit.kilojoules_per_mole / unit.radians**2  # Default k value if not specified
    targetTorsion_degrees: float = parameters["phi0"] * unit.degrees
    
    ## Add the atom quartet and settings to the torsion restraint
    torsionRestraint.addTorsion(restraintTorsionAtoms[0], restraintTorsionAtoms[1],
                                restraintTorsionAtoms[2], restraintTorsionAtoms[3],
                                [kForceConstant, targetTorsion_degrees])
    
    return system


###########################################################################################
def clear_all_restraints(saveXml: FilePath) -> None:
    """
    Remove all custom force constants from the given XML file.

    Parameters:
        saveXml (str): The path to the XML file.

    Returns:
        None
    """

    # Parse the XML file
    tree: ET.ElementTree = ET.parse(saveXml)
    root: ET.Element = tree.getroot()
    ## get the parameters section of the XML file
    parametersElement: ET.Element = root.find("Parameters")

    # Safely remove any custom force constants
    ## NB. our custom restrants are the only ones that use "k"
    paramsToPop: list = []
    if parametersElement is not None:
        for param in parametersElement.attrib:
            if param.startswith("k"):
                paramsToPop.append(param)
    for param in paramsToPop:
        parametersElement.attrib.pop(param)

    # Write the modified XML file
    tree.write(saveXml)
    
###########################################################################################