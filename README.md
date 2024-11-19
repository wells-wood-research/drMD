# :medical_symbol: drMD :medical_symbol:
Automated workflow for running molecular dynamics simulations with Amber and Openmm
# :medical_symbol: README Contents :medical_symbol:

1. [GitHub Installation (recommended)](#github-installation)
2. [Pip Installation (for advanced users)](#pip-installation)
3. [Usage](#usage)
4. [Config Syntax](#config-syntax)
    - [pathInfo](#pathinfo)
      - [inputDir](#inputdir)
      - [outputDir](#outputdir)
    - [hardwareInfo](#hardwareinfo)
      - [platform](#platform)
      - [parallelCPU](#parallelcpu)
      - [subprocessCpus](#subprocesscpus)
    - [miscInfo](#miscinfo)
      - [pH](#pH)
      - [firstAidMaxRetries](#firstaidmaxretries)
      - [boxGeometry](#boxgeometry)
      - [writeMyMethodsSection](#writemymethodssection)
    - [ligandInfo](#ligandinfo)
      - [ligandName](#ligandname)
      - [protons](#protons)
      - [charge](#charge)
      - [toppar](#toppar)
      - [mol2](#mol2)
    - [simulationInfo](#simulationinfo)
      - [stepName](#stepname)
      - [simulationType](#simulationtype)
      - [temperature](#temperature)
      - [temperatureRange](#temperaturerange)
      - [maxIterations](#maxiterations)
      - [duration](#duration)
      - [timestep](#timestep)
      - [logInterval](#loginterval)
    - [postSimulationInfo](#postsimulationinfo)
      - [endPointInfo](#endpointinfo)
        - [stepNames](#stepnamesendpoint)
        - [removeAtoms](#removeatomsendpoint)
      - [clusterInfo](#clusterinfo)
        - [stepNames](#stepnamescluster)
        - [removeAtoms](#removeatomscluster)
        - [nClusters](#nclusters)
        - [clusterBy](#clusterby)
    - [drMD Selection Syntax](#drmdselectionsyntax)
      - [keyword](#keyword)
      - [customSelection](#customselection)
         - [CHAIN_ID](#chainid)
         - [RES_NAME](#resname)
         - [RES_ID](#resid)
         - [ATOM_NAME](#atomname)
    - [Adding Restraints in drMD](#addingrestraints)
      - [restraintInfo](#restraintinfo)
        - [restraintType](#restrainttype)
        - [parameters](#parameters)
          - [k](#k)
          - [r0](#r0)
          - [theta0](#theta0)
          - [phi0](#phi0)
        - [selection](#selectionrestraints)

    - [Running Metadynamics with drMD](#runningmetadynamics)
      - [metaDynamicsInfo](#metadynamicsinfo)
        - [height](#height)
        - [biasFactor](#biasfactor)
        - [frequency](#frequency)
        - [biases](#biases)
          - [biasVar](#biasvar)
          - [minValue](#minvalue)
          - [maxValue](#maxvalue)
          - [biasWidth](#biaswidth)
          - [selection](#selectiometadynamics)


# GitHub Installation 
We reccomned that you use the following steps to install drMD:
1. Clone this repository
```bash
git clone https://github.com/ESPhoenix/drMD.git
```
2. Create and activate conda environment
```bash
conda create -n drMD python=3.10
```
```bash
conda activate drMD
```
3. Install AmberTools (needs to be before OpenMM) with conda
```bash
conda install -c conda-forge ambertools=23
``` 
4. Install OpenMM with conda
```bash
conda install -c omnia openmm
``` 
5. Install OpenBabel with conda
```bash
conda install -c conda-forge openbabel
```
6. Install other python libraries with pip
```bash
pip install -r requirements.txt
```

# Pip Installation
If you want to intergate drMD into a python-based pipeline, you can install drMD with pip and use it as a python module:

1. Create and activate conda environment
```bash
conda create -n drMD python=3.10
```
```bash
conda activate drMD
```
2. Install drMD with pip
```bash
pip install drMD
```
3. Install AmberTools (needs to be before OpenMM) with conda
```bash
conda install -c conda-forge ambertools=23
``` 
4. Install OpenMM with conda
```bash
conda install -c omnia openmm
``` 
5. Install OpenBabel with conda
```bash
conda install -c conda-forge openbabel
```


# Usage

Now that you have sucessfully set up the dependancies for drMD, you are nearly ready to run some biomolecular simulations!

If you have used the GitHub installation method, you can run drMD using the following command:

```bash
python /path/to/drMD.py --config config.yaml
```

If you have used the Pip installation method, you can import drMD as a python module, and as following:

```python
import drMD

myBatchConfig = "/path/to/config.yaml"

drMD.main(myBatchConfig)
```

This config file contains all of the user inputs drMD needs to run a series of biomolecular simulations.
The following section will detail the correct formatting of this config.yaml file

# Config syntax
The config.yaml file is in the YAML format *(https://en.wikipedia.org/wiki/YAML)* 
Inputs are grouped by theme and are stored as nested dictionaries and lists.
The next few sections will detail the correct formatting of the config.yaml file
<a id="pathinfo"></a>
## :brain: pathInfo
The **pathInfo** entry in the config file is a dictionary containing two parameters:
<a id="inputdir"></a>
### :anatomical_heart: inputDir
*(DirectoryPath)* This is the absoloute path towards a directory containing PDB files that will be used as starting points for your simulations.
            
  > :medical_symbol:
  > **To Perform Replicate** simulations, simply create copies of your starting PDB files in the inputDir, with each copy
  > named with a unique number. For example, your inputDir could contain my_protein_1.pdb, my_protein_2.pdb, etc.

<a id="outputdir"></a>
### :anatomical_heart:  outputDir  
*(DirectoryPath)*  This is the absoloute path towards a directory that you want your drMD outputs to be written to.

  > :medical_symbol:
  > The outputDir will be created if it does not already exist at the point of running drMD

  > :medical_symbol:
  > Within outputDir, a directory will be created for each PDB file contained in inputDir, in this document, these subdirectories will be refered to
  > as **runDirs**

Example pathInfo:
```yaml
pathInfo:
  inputDir: "/home/esp/scriptDevelopment/drMD/01_inputs"
  outputDir: "/home/esp/scriptDevelopment/drMD/02_outputs"
```
<a id="hardwareinfo"></a>
## :brain: hardwareInfo
This config entry tells drMD about your computer hardware and how you want to use it to run your simulations
The **hardwareInfo** entry in the config file is a dictionary containing three parameters:

<a id="inputdir"></a>
### :anatomical_heart:  platform
*(str)* This is the platform that will be used to run simulations in OpenMM. Accepted arguments for **platform** are *"CUDA"*, *"OpenCL"*, and *"CPU"*

  > :medical_symbol:
  > If you have access to GPU acceleration using CUDA, we recommend this option. If you cant use CUDA but have access to OpenCL, this is a close second.
  > If you don't have a GPU you can use the CPU option, this will be a lot slower.
  > Energy minimisation calculations do not benefit from GPU acceleration, so you should use the CPU option for these

<a id="parallelcpu"></a>
### :anatomical_heart:  parallelCPU
*(int)* This is the number  of simulations that will be run in paralell

<a id="subprocesscpus"></a>
### :anatomical_heart:  subprocessCpus
 *(int)* This is the number of cpu cores that will be allocated to each simulation.  

  > :medical_symbol:
  > The total CPU usage will be parallelCPU * subprocessCpus, so make sure you have enough CPUs when you set these parameters

Example hardwareInfo:
```yaml
hardwareInfo:
  parallelCPU: 16
  platform: "CUDA"
  subprocessCpus: 2
```
This will use CUDA git achive GPU acceleration and run 16 simulations in paralell using 2 cores each for a total useage of 32 cores.
---
<a id="miscinfo"></a>
## :brain: miscInfo
This section allows you to set some general options for your simulations:

<a id="ph"></a>
### :anatomical_heart:  pH
 *(int or float)* This is the pH of your simulation, this will affect the protonation states of your protein and any ligands present in your simulation

<a id="firstaidmaxretries"></a>
### :anatomical_heart:  firstAidMaxRetries
*(int)* This is the maximum number of times that drMD will attempt to recover from an error in a simulation

> :medical_symbol: This option can be very helpful for rescuing crashed simulations. However 
> don't rely on it too much. If your simulation keeps crashing you may want to reduce the 
> temperature or timestep parameters instead to make it more stable

<a id="boxgeometry"></a>
### :anatomical_heart:  boxGeometry 
*(str)*  This is the shape of the solvation box that will be used in your simulations. Accepted arguments for **boxGeometry** are *"cubic" or "octahedral"

<a id="writemymethodsection"></a>
### :anatomical_heart:  writeMyMethodsSection
*(bool)* If set to TRUE, drMD will automatically write a methods section for you to use in your publications or thesis.

> :medical_symbol: drMD methods sections contain all of the information one might need to replicate your simulations.
> The formatting of these methods section may be too robotic and repetative for you, feel free to reformat them as you see fit. 

Example miscInfo:
```yaml
miscInfo:
  pH: 7.4
  firstAidMaxRetries: 10
  boxGeometry: "cubic"
  writeMyMethodsSection: True
```
Simulations will be run with a pH of 7.4 in a cubic solvation box. The maximum number of first-aid retries will be 10. A methods section will automatically be generated. 


<a id="ligandinfo"></a>
## :brain: ligandInfo
The **ligandInfo** entry in the config file is optional and may be used if your PDB files have organic ligand or cofactors.
These small molecules will not have parameters in the AMBER forcefield, drMD will run an automated protocol to generate these parameters for you.
To do this, you will need to tell drMD some things about each ligand you whish tp simulate.

> :medical_symbol:
> The **ligandInfo** entry is *optional*. drMD will automatically detect ligand in your PDB files. It will also detect
> parameter files in your input directory. If you have frcmod and mol2 files for your ligand already made, they must be located in your **inputDir**

**ligandInfo** is a list of dictionaries that contain the following parameters:

<a id="ligandname"></a>
### :anatomical_heart:  ligandName
*(str)*  This is the three letter name of the ligand, this will be used to match the residue names in your PDB files

<a id="protons"></a>
### :anatomical_heart:  protons
  *(bool)*  This is a to tell drMD whether you have protons on your ligand. 
              If set to FALSE, drMD will run an automated protonation protocol to add protons to your ligand

  > :medical_symbol:
  > The automatic protonation protocol only works reliably for simple organic ligands.

  > :medical_symbol:
  > For more complex ligand, we recommended that you manually add protons in your input PDB file prior to running drMD

<a id="charge"></a>
### :anatomical_heart:  charge
*(int)*  This is the formal charge of the ligand 

<a id="toppar"></a>
### :anatomical_heart:  toppar
*(bool)*  This is to tell drMD whether you have an frcmod file for your ligand already made.
                If you already have one, it must be located in the 01_ligand_parameters directory within your outputDir

<a id="mol2"></a>
### :anatomical_heart:  mol2
*(bool)*   This is to tell drMD whether you have a mol2 file for your ligand already made.
                If you already have one, it must be located in the 01_ligand_parameters directory within your outputDir

Example ligandInfo:
```yaml
ligandInfo:
  - ligandName: "FMN"
    protons: True
    charge: -1
    toppar: False
    mol2: False
  - ligandName: "TPA"
    protons: True
    charge: -2
    toppar: False
    mol2: False
```
This **ligandInfo** tells drMD to expect two ligands: FMN and TPA. FMN has a formal charge of -1 and TPA has a formal charge of -2. Both ligands already have protons, so drMD will not add any. For both ligands the toppar and mol2 parameters are set to False, drMD will automatically generate these files for you

---

<a id="simulationinfo"></a>
## :brain: simulationInfo
This is the real meat and potatoes of the drMD config file. 

The **simulationInfo** entry in the config file is a list of dictionaries containing information about each simulation.

Each simulation detailed in **simulationInfo** will be run in sequence, with the output of the previous simulation being the starting point for the next simulation.
Each simulation dictionary contains the following parameters:

<a id="stepnamesiminfo"></a>
### :anatomical_heart:  stepName
*(str)* This is the name of the step that will be used to create a subdirectory in the runDir, we reccomend prefixing these names with numbers to make them order nicely

<a id="simulationtype"></a>
### :anatomical_heart: simulationType
*(str)* This is the type of simulation that will be run. Accepted arguments are:

    - **"EM"**:         This will run a steepest-decent Energy Minimisation step. 
    > :medical_symbol:
    > We reccomended that you run one of these steps before any other simulation steps
    - **"NVT"**:        This will run an NVT (constant volume) molecular dynamics simulation
    - **"NPT"**:        This will run an NPT (constant pressure) molecular dynamics simulation
    > :medical_symbol:
    > For the majority of protein simulations, the NPT ensemble is used for production MD simulations, while the NVT ensemble is only used in equilibration steps
    - **"META"**:       This will run a Metadynamics simulation 

### Selecting simulation temperature 
For most simulations, a constant temperature is used. In this case the following parameter is required:

<a id="temperature"></a>
#### :anatomical_heart: temperature
*(int)* This is the temperature of the simulation in Kelvin 

If you wish to change the temperature throughout the simulation, the following parameter is required:


<a id="temperaturerange"></a>
#### :anatomical_heart: temperatureRange
*(list of int)* This is a list of integers (again, in Kelvin) that will be used to change the temperature throughout the simulation. 

### Energy Minimisation Pararameters
For Energy Minimisation steps, the following additional parameters are required:

<a id="maxiterations"></a>
#### :anatomical_heart: maxIterations
*(int)* This is the maximum number of iterations that will be run in the Energy Minimisation step.
If this parameter is set to -1, the step will run until the energy converges.

Example Energy Minimisation syntax:
```yaml
simulationInfo:
  - stepName: "01_energy_minimisation"
    type: "EM"
    temp: 300
    maxIterations: -1
```
This will run a energy minimisation until the energy converges

### Generic Simulation Parameters
For "normal" MD simulations using NVT or NpT ensembles, as well as for Metadynamics simulations, the following additional parameters are required:

<a id="duration"></a>
#### :anatomical_heart: duration
: *(str)* This is the duration of the simulation step, as a string "int unit" eg. "1000 ps"

<a id="timestep"></a>
#### :anatomical_heart: timestep
 *(str)* This is the timestep of the simulation, as a string "int unit" eg. "2 fs"

<a id="loginterval"></a>
#### :anatomical_heart: logInterval
*(str)* This is the frequency that the simulation will write to file using built-in OpemMM reporters. As a string "int unit" eg. "100 ps"

Example NVT simulation syntax:
```yaml
simulationInfo:
  - stepName: "02_NVT_pre-equilibraition"
    type: "NVT"
    duration: "100 ps"
    timestep: "2 fs"
    temp: 300
    logInterval: "10 ps"
```
This will run a 100 ps NVT molecular dynamics simulation with a timestep of 2 fs, a temp of 300 and a logInterval of 10 ps

---

### Post-simulation processing
After all of your simulations have been run, drMD contains some simple utilities for organising your output files and deleting any unwanted files.

If you want to do any post-processing, you will need to provide the following parameter in your config file:

<a id="postsimulationinfo"></a>
### :brain: postSimulationInfo
 *(dict)* This is a dictionary containing the parameters for the post-simulation processing

If you wish to collect PDB files that represent the last frame of each simulation, you may include the following parameter in **postSimulationInfo**:

<a id="endpointinfo"></a>
### :brain: endpointInfo
*(dict)* This is a dictionary containing the parameters the following parameters:

<a id="stepnamesendpoint"></a>
#### :anatomical_heart: stepNames
*(list)* This is a list of strings containing the names of the steps in the simulation, these should match the stepNames that you have used in your simulationInfo dictionary (described above). Endpoint PDB files will be gathered for these steps

<a id="removeatomsendpoint"></a>
#### :anatomical_heart: removeAtoms
*(list)* This is a list of dictionaries containing the selections of atoms to be removed from the PDB files. For a full description of how to do this, see [drMD Selection syntax](#drmd-selection-syntax)

Molecular Dyamics simulations can generate very large output files that can become rather unweildy and difficult to anaylse. One way to quickly see the most important parts of your simulation is to perform clustering on your simulation trajectories. To do this with drMD, include the following parameter in your config file:

<a id="clusterinfo"></a>
### :brain: clusterInfo
*(dict)* This is a dictionary containing the parameters following parameters:

<a id="stepnamescluster"></a>
#### :anatomical_heart: stepNames
*(list)* This is a list of strings containing the names of the steps in the simulation, these should match the stepNames that you have used in your simulationInfo dictionary (described above). Clustering will be performed on trajectories of these steps

<a id="nclusters"></a>
#### :anatomical_heart: nClusters
*(int)* This is the number of clusters PDB files that will be generated

<a id="clusterby"></a>
#### :anatomical_heart: clusterBy
*(list)* This is a list of selections of atoms to cluster on. If you want to explore the motions of one particular group of atoms in your system (e.g. the backbone of a protein), you can include a selection of these atoms in this parameter. For a full description of how to do this, see [drMD Selection syntax](#drmd-selection-syntax)

<a id="removeatomscluster"></a>
#### :anatomical_heart: removeAtoms
*(list)* This is a list of dictionaries containing the selections of atoms to be removed from the output cluster PDB files.  For a full description of how to do this, see [drMD Selection syntax](#drmd-selection-syntax)


<a id="collatevitalsreports"></a>
#### :anatomical_heart: collateVitalsReports
*(bool)* If True, will collate vitals reports from the trajectories generated by the MD simulations into the **00_vitals_reports** directory in your specified output directory.

---

<a id="drmdselectionsyntax"></a>
## drMD Selection syntax
When creating restraints, metadynamics bias variables or running post-simulation clustering, you will need to specify the selection of atoms that the restraints will be applied to. To do this, you will need to supply a "selection" dictionary. This dictionary must contain the following parameter:

<a id="keyword"></a>
### :anatomical_heart: keyword
*(str)* This is the keyword that will be used to specify the selection. Accepted arguments are:
  - **"protein"** : This will select all protein atoms in your system
  - **"water"** : This will select all water molecules in your system
  - **"ions"**: This will select all ions in your system
  - **"ligand"** : This will select all non-protein, non-water and non-ion atoms in your system
  - **"custom"** : This will select all atoms that match the customSelection (see below)

Example use of keywords in the selection dictionary:

```yaml
selection:
  keyword: "water"
```

This will select all water molecules in your system


If you have used the **custom** keyword, you will need to use an additional parameter the selection dictionary:
<a id="customselection"></a>
### :anatomical_heart: customSelection
*(list)* This is a list of dictionaries containing details of the atoms that will be selected.
Each dictionary in the list must contain the following parameters:

<a id="chainid"></a>
#### :anatomical_heart: CHAIN_ID
 *(str, list of str, or "all")* This is the chain ID of the atom to be selected

 <a id="resname"></a>
#### :anatomical_heart: RES_NAME
 *(str, llist of str, or "all")*  This is the three-letter residue name of the atom to be selected

<a id="resid"></a>
#### :anatomical_heart: RES_ID
 *(int, list of int, or "all")* This  is the residue ID of the atom to be selected

<a id="atomname"></a>
#### :anatomical_heart: ATOM_NAME
*(str, list of str, or "all")*  This is the atom name of the atom to be selected


For the above parameters:
- if a single string (or int for **RES_ID**) is provided, the selection will match that value for the given parameter
- if a list is provided, the selection will match any value in the list for the given parameter
- if the wildcard string "all" is provided, the selection will match all values for the given parameter

:biohazard: WARNING: PLEASE DON'T CALL ANY OF YOUR RESIDUE NAMES IN YOUR PDB FILE "ALL" AS THIS WILL CAUSE PROBLEMS

> :medical_symbol:
> This selection method is used to match atoms using columns of your PDB files. To work out how to identify what inputs you need to select your atoms of interest, simply open your input PDB file in PyMOL and use labels to identify the relavent **CHAIN_ID**, **RES_NAME**, **RES_ID**, and **ATOM_NAME** values

Example customSelection syntax, (this is probably more complex than you will need, but shows off what you can do):

```yaml
selection:
  keyword: "custom"
  customSelection:
    - {CHAIN_ID: ["A","B"], RES_NAME: "all", RES_ID:  "all", ATOM_NAME: "CA"}
    - {CHAIN_ID: "A", RES_NAME: "SER", RES_ID: "131", ATOM_NAME: ["CB", "HB1", "HB2", "OG", "HG1"]}
    - {CHAIN_ID: "C", RES_NAME: "FMN", RES_ID: "all", ATOM_NAME: "all"}
```
In the above example, a selection containing the following:
- all CA atoms in chains A and B
- the atoms CB, HB1, HB2, OG, and HG1 in residue Ser131 in chain A 
- all atoms in the residue FMN in chain C
---
<a id="addingrestraints"></a>
### Adding Restraints with drMD
If you whish to perform simulations with restraints, create a *restraintsInfo* dictionary in the simulation step:

<a id="restraintinfo"></a>
#### :brain: restraintInfo
*(list of dict)*  This is a list of dictionaries containing information about each restraints. 

Within the restraintInfo list, you must provied at least one dictionary that contains the following parameters:

<a id="restrainttype"></a>
##### :anatomical_heart: restraintType
 *(str)* This is the type of restraints that will be added. Accepted arguments are: "distance", "angle", "dihedral", "position"

<a id="parameters"></a>
##### :anatomical_heart: parameters
*(dict)*  This is a dictionary containing the parameters for the restraints.

Within the parameters dictionary, all restraint types require the following parameter:

<a id="k"></a>
###### :anatomical_heart:  k 
*(int)*  This is the force constant of the restraint (int), given in kJ/mol A^-2


Additional entries in the parameters dictionary depend on the type of restraints:

<a id="r0"></a>
###### :anatomical_heart: r0
 *(int or float)*  *Required for distance restraints*. This is the distance in Angstroms that the restraint should be applied to 

<a id="theta0"></a>
###### :anatomical_heart: theta0
*(int or float)*   *Required for angle restraints*. This is the angle in degrees that the angle should be constrained to 

<a id="phi0"></a>
###### :anatomical_heart: phi0
*(int or float)*  *Required for torsion restraints*. This is the angle in degrees that the dihedral should be constrained to

All restraints require the selection parameter. This tells drMD what atoms to apply the restraint to

<a id="selectionrestraints"></a>
#### :anatomical_heart: selection
    - **selection**:  *(list of dicts)*  This is a dictionary containing information on the selection of atoms that the restraints will be applied to.

> :medical_symbol:
>The selection method is shared between multiple different inputs in the drMD config file. This is described in more detail in the next section


Example restraints syntax:
```yaml
    restraints:
    - restraintType: "position"
      selection:
        keyword: "protein"
        parameters:
          k: 1000
    - type: "distance"
      selection: 
        keyword: "custom"  
        customSelection:
          - {CHAIN_ID: "A", RES_NAME: "ALA", RES_ID: 1, ATOM_NAME: "CA"}
          - {CHAIN_ID: "A", RES_NAME: "ALA", RES_ID: 2, ATOM_NAME: "CA"}
      parameters:
        k: 1000
        r0: 3

```

This example will add the following restraints:
- Position restraints to the protein atoms with a force constant of 1000 kJ/mol 
- A 3 Angstrom distance restraint between the CA atoms of residues 1 and 2 of the protein, with a force constant of 1000 kJ/mol

For a detailed explaination of how to select chains, residues, and atoms for restraints, see the [drMD Selection syntax](#drmd-selection-syntax) section.

---

<a id="runningmetadynamics"></a>
### Running Metadynamics with drMD
To run a metadynamics simulation, first set the **simulationType** to **"META"**.

Once you have selected your **simulationType**, you will need to include an additional **metaDynamicsInfo** dictionary in your simulation dictionary:

<a id="metadynamicsinfo"></a>
### :brain: metaDynamicsInfo
*(dict)* This is a dictionary containing the parameters for the Metadynamics simulation,

Within the **metaDynamicsInfo** dictionary, you must provide the following parameters:

<a id="height"></a>
#### :anatomical_heart: height
*(int)* This is the height  parameter used in the Metadynamics simulation

<a id="biasfactor"></a>
#### :anatomical_heart: biasFactor
*(int)* This is the bias factor  parameter used in the Metadynamics simulation

<a id="frequency"></a>
#### :anatomical_heart: frequency
*(int)* How often (in time steps) gaussians will be added to the bias potential

<a id="biases"></a>
#### :anatomical_heart: biases
*(list of dicts)* This is a list of dictionaries containing information about each biasVariable.

Within each dictionary in **biases** you must provide the following parameters:

<a id="biasvar"></a>
###### :anatomical_heart: biasVar
*(str)* This is the type of biasVariable that will be added.

> Accepted arguments for **biasVar** are **"RMSD"**, **"torsion"** **"distance"** and **"angle"**

<a id="minvalue"></a>
###### :anatomical_heart: minValue
*(float)* This is the minimum value of the biasVariable.

<a id="maxvalue"></a>
###### :anatomical_heart: maxValue
*(float)* This is the maximum value of the biasVariable.

<a id="biaswidth"></a>
###### :anatomical_heart: biasWidth
*(float)* This determines the width of gaussians added to the bias potential

<a id="selectiometadynamics"></a>
##### :anatomical_heart: selection
*(dict)*  This is a dictionary containing information on the selection of atoms that the biasVariable will be applied to. The selection syntax is identical to that used for the restraints. For a full description of how to do this, see [drMD Selection syntax](#drmd-selection-syntax)

> :medical_symbol:
> Depending on the type of bias variable, different numbers of atoms must be selected:
> - Distance bias variables require two atoms to be selected
> - Angle bias variables require three atoms to be selected
> - Torsion bias variables require four atoms to be selected

> :medical_symbol:
> You will also need to specify at least one biasVariable for the simulation to sample.
> You can specify as many biasVariables as you wish, with one dictionary per biasVariable


Example MetaDynamics syntax:
```yaml
    metaDynamicsInfo:
      height: 2
      biasFactor: 5
      frequency: 50
      biases: 
        - biasVar: "RMSD"
          minValue: 0
          maxValue: 10
          biasWidth: 1
          selection: 
            keyword: "backbone"
        - biasVar: "torsion"
          selection: 
            keyword: "custom"
            customSelection:
            - {CHAIN_ID: "A", RES_NAME: "ALA", RES_ID: 1, ATOM_NAME: "CA"}
            - {CHAIN_ID: "A", RES_NAME: "ALA", RES_ID: 2, ATOM_NAME: "CA"}
            - {CHAIN_ID: "A", RES_NAME: "ALA", RES_ID: 3, ATOM_NAME: "CA"}
            - {CHAIN_ID: "A", RES_NAME: "ALA", RES_ID: 4, ATOM_NAME: "CA"}

```

This example will add a RMSD bias to the backbone of the protein
and a torsion bias between the CA atoms of residues 1, 2, 3, and 4 of the protein
---



## Advanced YAML-ing with variables
If you are running multiple simulation steps that share the same parameters, you can use variables in the YAML config file. This will most commonly come up when you are applying position restraints during equilibriation steps. Below is a standard syntax for a pair of equilibiation steps:
```yaml
simulationInfo:
  - stepName: "01_NVT_pre-equilibraition"
    simulationType: "NVT"
    duration: "100 ps"
    timestep: "4 fs"
    heavyProtons: True
    temperature: 300
    logInterval: "10 ps"
    restraintInfo: 
    - restraintType: "position"
      parameters:
        k: 1000
      selection:
        keyword: "protein"
    - restraintType: "position"
      selection: 
        keyword: "ligand"  
      parameters:
        k: 1000

  - stepName: "02_NPT_pre-equilibraition"
    simulationType: "NPT"
    duration: "100 ps"
    timestep: "4 fs"
    heavyProtons: True
    temperature: 300
    logInterval: "10 ps"
    restraintInfo:
      - restraintType: "position"
        parameters:
          k: 1000
        selection:
          keyword: "protein"
      - restraintType: "position"
        selection: 
          keyword: "ligand"  
        parameters:
          k: 1000

``` 
Instead of repeating the restraintsInfo section each time, **equilibriumRestraints** can be defined as a variable, then re-used in each equilibriation step:

```yaml
equilibriationRestraints: &equilibriationRestraints
    - restraintType: "position"
      parameters:
        k: 1000
      selection:
        keyword: "protein"

    - restraintType: "position"
      selection: 
        keyword: "ligand"  
      parameters:
        k: 1000

simulationInfo:
  - stepName: "01_NVT_pre-equilibraition"
    simulationType: "NVT"
    duration: "100 ps"
    timestep: "4 fs"
    heavyProtons: True
    temperature: 300
    logInterval: "10 ps"
    restraintInfo: *equilibriationRestraints

  - stepName: "02_NPT_pre-equilibraition"
    simulationType: "NPT"
    duration: "100 ps"
    timestep: "4 fs"
    heavyProtons: True
    temperature: 300
    logInterval: "10 ps"
    restraintInfo: *equilibriationRestraints

```

# WORKED EXAMPLES
In this section we will go through a series of configuration files 