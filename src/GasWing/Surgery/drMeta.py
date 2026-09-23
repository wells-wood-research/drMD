## BASIC PYTHON LIBRARIES
import os
from os import path as p
import numpy as np
import pandas as pd

## OPENMM LIBRARIES
import openmm.app as app
from openmm.app import metadynamics
import openmm as openmm
import  openmm.unit  as unit

## drMD LIBRARIES
from Surgery import drSim, drRestraints, drFirstAid
from ExaminationRoom import drLogger, drCheckup
from UtilitiesCloset import drSelector, drFixer, drMethodsWriter

########################################################################################################
########################################################################################################
@drLogger.monitor_progress_decorator()
@drFirstAid.firstAid_handler()
@drCheckup.check_up_handler()
def run_metadynamics(prmtop: app.Topology,
                      inpcrd: any,
                        sim: dict,
                          saveFile: str,
                            outDir: str,
                              platform: openmm.Platform,
                                refPdb: str,
                                    config:dict) -> str:

    """
    Run a simulation at constant pressure (NpT) step with biases.

    Args:
        prmtop (str): The path to the topology file.
        inpcrd (str): The path to the coordinates file.
        sim (dict): The simulation parameters.
        saveFile (str): The path to the checkpoint or XML file.
        simDir (str): The path to the simulation directory.
        platform (openmm.Platform): The simulation platform.
        pdbFile (str): The path to the PDB file.

    Returns:
        str: The path to the XML file containing the final state of the simulation.

    This function runs a metadynamics simulation at constant pressure (NpT) step with biases.
    It initializes the system, handles any constraints, adds a constant pressure force,
    sets up the integrator and system, loads the state from a checkpoint or XML file,
    sets up reporters, runs the simulation, saves the final geometry as a PDB file,
    resets the chain and residue Ids, runs drCheckup, performs trajectory clustering
    if specified in the simulation parameters, and saves the simulation state as an
    XML file.
    """
    stepName = sim["stepName"]
    drLogger.log_info(f"Running MetaDynamics Step: {stepName}", True)
    drMethodsWriter.add_simulation_step_to_log(stepName)
    drMethodsWriter.add_parameter_to_simulation_log(stepName, "simulationType", str(sim["simulationType"]))

    ## make a simulation directory
    simDir: str = p.join(outDir, stepName)
    os.makedirs(simDir, exist_ok=True)

    sim = drSim.process_sim_data(sim)
    # Define the nonbonded method and cutoff.
    nonbondedMethod: openmm.NonbondedForce = app.NoCutoff
    drMethodsWriter.add_parameter_to_simulation_log(stepName, "nonbondedMethod", "NoCutoff")

    # Define the restraints.
    hBondconstraints: openmm.Force = app.HBonds
    drMethodsWriter.add_parameter_to_simulation_log(stepName, "constraints", str(hBondconstraints))
    # Create the system.
    system: openmm.System = prmtop.createSystem(nonbondedMethod=nonbondedMethod,
                                                constraints=hBondconstraints)

    # Deal with restraints (clear all lurking restraints and constants)
    system: openmm.System = drRestraints.restraints_handler(system, prmtop, inpcrd, sim, saveFile, refPdb)
    system: openmm.System = scale_dielectric_constant(system)

    drMethodsWriter.add_parameter_to_simulation_log(stepName, "temperature", str(sim["temperature"]))
    # Read metaDynamicsInfo from sim config
    metaDynamicsInfo: dict = sim["metaDynamicsInfo"]


    # Read metaDynamicsInfo from sim config
    metaDynamicsInfo: dict = sim["metaDynamicsInfo"]
    # Read biases from sim config and create bias variables
    biases: list = metaDynamicsInfo["biases"]
    biasVariables: list = []

    for bias in biases:
        # Get atom indexes and coordinates for the biases
        atomIndexes: list = drSelector.get_atom_indexes(bias["selection"], refPdb)
        atomCoords: list = get_atom_coords_for_metadynamics(prmtop, inpcrd)
        # Create bias variable based on the type of bias
        if bias["biasVar"].upper() == "RMSD":
            biasVariable, addedBias= gen_rmsd_bias_variable(bias, atomCoords, atomIndexes)
            biasVariables.append(biasVariable)
        elif bias["biasVar"].upper() == "TORSION":
            biasVariable, addedBias= gen_dihedral_bias_variable(bias, atomCoords, atomIndexes)
            biasVariables.append(biasVariable)
        elif bias["biasVar"].upper() == "DISTANCE":
            biasVariable, addedBias= gen_distance_bias_variable(bias, atomCoords, atomIndexes)
            biasVariables.append(biasVariable)
        elif bias["biasVar"].upper() == "ANGLE":
            biasVariable, addedBias= gen_angle_bias_variable(bias, atomCoords, atomIndexes)
            biasVariables.append(biasVariable)
        elif bias["biasVar"].upper() == "RG":
            biasVariable, addedBias= gen_gyration_bias_variable(bias, atomCoords, atomIndexes)
            biasVariables.append(biasVariable)

    drMethodsWriter.add_parameter_to_simulation_log(stepName, "biases", addedBiases)    

        
    meta: metadynamics.Metadynamics = metadynamics.Metadynamics(
        system=system,
        variables=biasVariables,
        temperature=sim["temperature"],
        biasFactor=metaDynamicsInfo["biasFactor"],
        height=metaDynamicsInfo["height"],
        frequency=50,
        saveFrequency=50,
        biasDir=simDir,
    )

    drMethodsWriter.add_parameter_to_simulation_log(stepName, "biasFactor", str(metaDynamicsInfo["biasFactor"]))
    drMethodsWriter.add_parameter_to_simulation_log(stepName, "height", str(metaDynamicsInfo["height"]))
    drMethodsWriter.add_parameter_to_simulation_log(stepName, "frequency", str(50))
    drMethodsWriter.add_parameter_to_simulation_log(stepName, "saveFrequency", str(50))

    # Set up integrator
    integrator: openmm.LangevinMiddleIntegrator = openmm.LangevinMiddleIntegrator(
        sim["temperature"], 1/unit.picosecond, sim["timestep"]
    )
    drMethodsWriter.add_parameter_to_simulation_log(stepName, "integrator", "LangevinMiddleIntegrator")
    drMethodsWriter.add_parameter_to_simulation_log(stepName, "temperature", str(sim["temperature"]))
    drMethodsWriter.add_parameter_to_simulation_log(stepName, "timeStep", str(sim["timestep"]))
    drMethodsWriter.add_parameter_to_simulation_log(stepName, "friction", "1 ps^-1")

    # Create new simulation
    simulation: app.Simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)
    simulation: app.Simulation = drSim.load_simulation_state(simulation, saveFile)

    # Set up reporters
    totalSteps: int = simulation.currentStep + sim["nSteps"]
    reportInterval: int = sim["logInterval"]


    simulation: app.Simulation = drSim.init_reporters(
        simDir=simDir,
        nSteps=totalSteps,
        reportInterval=reportInterval,
        simulation=simulation,
        dcdAtomSelections=config["miscInfo"]["trajectorySelections"],
        refPdb=refPdb,
    )


    
    drMethodsWriter.add_parameter_to_simulation_log(stepName, "totalSteps", totalSteps)
    drMethodsWriter.add_parameter_to_simulation_log(stepName, "reportInterval", reportInterval)
    drMethodsWriter.add_parameter_to_simulation_log(stepName, "simulationTime", totalSteps * sim['timestep'])

    # Run metadynamics simulation
    meta.step(simulation, sim["nSteps"])
 
    # find name to call outFiles
    protName: str = p.basename(p.dirname(simDir))
    # save result as pdb - reset chain and residue Ids
    endPointPdb: str = p.join(simDir, f"{protName}.pdb")

    drSim.write_pdb(endPointPdb, simulation)
    # with open(endPointPdb, 'w') as output:
    #     app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)
    ## reset the chain and residue ids to original
    drFixer.reset_chains_residues(refPdb, endPointPdb, config)

    ## create a PDB file with the same atoms as the trajectory
    trajectoryPdb = p.join(simDir, "trajectory.pdb")
    drSelector.slice_pdb_file(config["miscInfo"]["trajectorySelections"], endPointPdb, trajectoryPdb)
    # save simulation as XML
    saveXml: str = p.join(simDir, f"{stepName}.xml")
    simulation.saveState(saveXml)
    drFixer.reset_chains_residues(refPdb, endPointPdb, config)

    ## get free energy and write to csv
    metadynamicsFreeEnergy = meta.getFreeEnergy()
    freeEnergyCsv = p.join(simDir, "freeEnergy.csv")
    np.savetxt(freeEnergyCsv, metadynamicsFreeEnergy, delimiter=",", header="Bias Variable Free Energy")

    # Return checkpoint file for continuing simulation
    return saveXml
########################################################################################################
def scale_dielectric_constant(system: openmm.System) -> openmm.System:
    """
    Scale the dielectric constant of the system.

    Parameters
    ----------
    system : openmm.System
        The system to scale.
    scaleFactor : float
        The factor by which to scale the dielectric constant.

    Returns
    -------
    system : openmm.System
        The system with the scaled dielectric constant.
    """
    scaleFactor: float = 0.8
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            originalDielectric: float = force.getReactionFieldDielectric()
            newDielectric: float = originalDielectric * scaleFactor
            force.setReactionFieldDielectric(newDielectric)
            drMethodsWriter.update_methods_log(f"Scaled dielectric constant from {originalDielectric} to {newDielectric}")
            break  # Assuming only one NonbondedForce exists
    return system
########################################################################################################
def get_atom_coords_for_metadynamics(prmtop: app.Topology, inpcrd: any) -> list:
    """
    Get the coordinates of all atoms in the system.

    Parameters
    ----------
    prmtop : file
        The topology file.
    inpcrd : file
        The coordinates file.

    Returns
    -------
    atomCoords : list
        The coordinates of all atoms in the system.
    """
    # Get the topology and positions from the input files
    topology: app.Topology = prmtop.topology
    positions: np.ndarray = inpcrd.positions

    # Initialize an empty list to store atom coordinates
    atomCoords: list = []

    # Loop over all atoms in the topology
    for i, atom in enumerate(topology.atoms()):
        # Append the coordinates of the current atom to the atomCoords list
        atomCoords.append(positions[i])

    # Return the list of atom coordinates
    return atomCoords
########################################################################################################
def gen_angle_bias_variable(bias: dict, atomCoords: np.ndarray, atomIndexes: list) -> metadynamics.BiasVariable:
    """
    Generate an angle bias variable.

    Parameters
    ----------
    bias : dict
        The bias dictionary containing the angle parameters.
    atomCoords : np.ndarray
        The coordinates of all atoms in the system.
    atomIndexes : list
        The indexes of the atoms involved in the angle calculation.

    Returns
    -------
    angleBiasVariable : openmm.CustomTorsionForce
        The angle bias variable.
    """

    # Create an angle bias force
    angleForce: openmm.CustomAngleForce = openmm.CustomAngleForce("theta")


    # Add the atom indexes for the angle
    angleForce.addAngle(atomIndexes[0],
                          atomIndexes[1],
                          atomIndexes[2])
    
    angleBiasVariable: metadynamics.BiasVariable = metadynamics.BiasVariable(force = angleForce,
                                                    minValue = bias["minValue"] * unit.degrees,
                                                    maxValue = bias["maxValue"] * unit.degrees,
                                                    biasWidth = bias["biasWidth"] * unit.degrees,
                                                     periodic = False)
    bias.update({"periodic": False})
    bias.update({"energy": angleForce.getEnergyExpression()})
    return angleBiasVariable, bias

########################################################################################################
def gen_dihedral_bias_variable(bias: dict, atomCoords: np.ndarray, atomIndexes: list) -> metadynamics.BiasVariable:
    """
    Generate a dihedral bias variable.

    Parameters
    ----------
    bias : dict
        The bias dictionary containing the dihedral angle parameters.
    atomCoords : list
        The coordinates of all atoms in the system.
    atomIndexes : list
        The indexes of the atoms involved in the dihedral angle.

    Returns
    -------
    dihedralBiasVariable : openmm.CustomTorsionForce
        The dihedral bias variable.
    """
    ## express dihedral as a harmonic potential
    dihedralEnergyExpression: str = "theta"
    
    ## generate a dihedral bias force:
    # Create a custom torsion force object
    dihedralForce: openmm.CustomTorsionForce = openmm.CustomTorsionForce(dihedralEnergyExpression)

    drMethodsWriter.update_methods_log(f"              force: {dihedralForce}")
 
    # Add the atom indexes for the dihedral angle
    dihedralForce.addTorsion(atomIndexes[0],
                              atomIndexes[1],
                                atomIndexes[2],
                                  atomIndexes[3])

    # Create a dihedral bias variable
    dihedralBiasVariable: metadynamics.BiasVariable = metadynamics.BiasVariable(force = dihedralForce,
                                                    minValue = bias["minValue"] * unit.degrees,
                                                    maxValue = bias["maxValue"] * unit.degrees,
                                                    biasWidth = bias["biasWidth"] * unit.degrees,
                                                    periodic = True)
    bias.update({"periodic": True})
    bias.update({"energy": dihedralForce.getEnergyExpression()})
    return dihedralBiasVariable, bias

########################################################################################################
def gen_distance_bias_variable(bias: dict, atomCoords: np.ndarray, atomIndexes: list) -> metadynamics.BiasVariable:
    """
    Generate a distance bias variable.

    Parameters
    ----------
    bias : dict
        The bias dictionary containing the distance parameters.
    atomCoords : np.ndarray
        The coordinates of all atoms in the system.
    atomIndexes : list
        The indexes of the atoms involved in the distance calculation.

    Returns
    -------
    distanceBiasVariable : openmm.CustomTorsionForce
        The distance bias variable.
    """

    # Create a distance bias force
    distanceForce: openmm.CustomBondForce = openmm.CustomBondForce("r")
    drMethodsWriter.update_methods_log(f"              force: {distanceForce}")

    distanceForce.addBond(atomIndexes[0],
                          atomIndexes[1])
    
    distanceBiasVariable: metadynamics.BiasVariable = metadynamics.BiasVariable(force = distanceForce,
                                                    minValue = bias["minValue"] * unit.angstrom,
                                                    maxValue = bias["maxValue"] * unit.angstrom,
                                                    biasWidth = bias["biasWidth"] * unit.angstrom,
                                                    periodic = False)
    bias.update({"periodic": False})
    bias.update({"energy": distanceForce.getEnergyExpression()})
    return distanceBiasVariable, bias
########################################################################################################
def gen_gyration_bias_variable(bias: dict, atomCoords: np.ndarray, atomIndexes: list) -> metadynamics.BiasVariable:
    """
    Generate a distance bias variable.

    Parameters
    ----------
    bias : dict
        The bias dictionary containing the distance parameters.
    atomCoords : np.ndarray
        The coordinates of all atoms in the system.
    atomIndexes : list
        The indexes of the atoms involved in the distance calculation.

    Returns
    -------
    distanceBiasVariable : openmm.CustomTorsionForce
        The distance bias variable.
    """

     # Generate a Rg bias force
    rgForce: openmm.RGForce = openmm.RGForce(atomIndexes)
    drMethodsWriter.update_methods_log(f"              force: {rgForce}")

    # Create a Rg bias variable
    gyrationBiasVariable: metadynamics.BiasVariable = metadynamics.BiasVariable(
        force=rgForce,
        minValue = bias["minValue"] * unit.angstrom,
        maxValue = bias["maxValue"] * unit.angstrom,
        biasWidth = bias["biasWidth"] * unit.angstrom,
        periodic=False)
    bias.update({"periodic": False})
    bias.update({"energy": rgForce.getEnergyExpression()})
    return gyrationBiasVariable, bias



########################################################################################################
def gen_rmsd_bias_variable(bias: dict, atomCoords: np.ndarray, atomIndexes: list) -> metadynamics.BiasVariable:
    """
    Generate a RMSD (Root Mean Square Deviation) bias variable.

    Parameters
    ----------
    bias : dict
        The bias dictionary containing the RMSD parameters.
    atomCoords : np.ndarray
        The coordinates of all atoms in the system.
    atomIndexes : list
        The indexes of the atoms involved in the RMSD calculation.

    Returns
    -------
    rmsdBiasVariable : openmm.RMSDForce
        The RMSD bias variable.
    """

    # Generate a RMSD bias force
    rmsdForce: openmm.RMSDForce = openmm.RMSDForce(atomCoords, atomIndexes)
    drMethodsWriter.update_methods_log(f"              force: {rmsdForce}")

    # Create a RMSD bias variable
    rmsdBiasVariable: metadynamics.BiasVariable = metadynamics.BiasVariable(
        force=rmsdForce,
        minValue = bias["minValue"] * unit.angstrom,
        maxValue = bias["maxValue"] * unit.angstrom,
        biasWidth = bias["biasWidth"] * unit.angstrom,
        periodic=False)
    bias.update({"periodic": False})
    bias.update({"energy": rmsdForce.getEnergyExpression()})
    return rmsdBiasVariable, bias



########################################################################################################
