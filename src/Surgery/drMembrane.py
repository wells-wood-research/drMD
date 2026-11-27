## BASIC PYTHON LIBRARIES
import os
from os import path as p
import sys
import glob
import json
import shutil
import subprocess

## drMD LIBRARIES
from ExaminationRoom import drLogger
from UtilitiesCloset import drSplash, drListInitiator, drDependencies

## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils

## CLEAN CODE
from typing import Dict, List, Optional, Tuple
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath

"""
drMembrane: builds protein-in-bilayer systems for drMD

Workflow:
    1. (optional) orient the protein with MemPrO
    2. pack lipids, water and salt around the protein with PACKMOL-Memgen
       (MemPrO or MEMEMBED orientation unless the protein is pre-oriented)
    3. return the packed PDB and box dimensions to drPrep.make_amber_params,
       which parameterises the system with tleap (Lipid21 + ff14SB/ff19SB + TIP3P)

Protonation is done by drMD (pdb2pqr / propka) before this step, so PACKMOL-Memgen
is always run with --notprotonate --nottrim.
"""
#########################################################################
def prepare_membrane_inputs(batchConfig: Dict) -> Dict:
    """
    Writes cleaned copies of the input PDB files to <outputDir>/00_membrane_inputs and
    points inputDir at them. User files are never modified.
        - removes OPM / PPM membrane dummy atoms (resname DUM) and crystallographic waters
        - keeps only ATOM / HETATM / TER records
        - fills blank occupancy / B-factor columns (OPM files leave these empty)
    """
    inputDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
    cleanDir: DirectoryPath = p.join(batchConfig["pathInfo"]["outputDir"], "00_membrane_inputs")
    os.makedirs(cleanDir, exist_ok=True)
    for pdbName in sorted(os.listdir(inputDir)):
        if not pdbName.endswith(".pdb"):
            continue
        with open(p.join(inputDir, pdbName)) as fin, open(p.join(cleanDir, pdbName), "w") as fout:
            for line in fin:
                if line.startswith("TER"):
                    fout.write("TER\n")
                    continue
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                ## OPM / PPM dummy membrane atoms and crystallographic waters are not simulated
                ## (the system is re-solvated by PACKMOL-Memgen)
                if line[17:20].strip() in ["DUM", "HOH", "WAT"]:
                    continue
                line = line.rstrip("\n").ljust(80)
                occupancy = line[54:60].strip() or "1.00"
                bFactor = line[60:66].strip() or "0.00"
                element = line[76:78].strip()
                if element == "":
                    atomName = line[12:16].strip()
                    element = atomName[1] if atomName[:1].isdigit() else atomName[0]
                fout.write(f"{line[:54]}{float(occupancy):6.2f}{float(bFactor):6.2f}          {element:>2}\n")
            fout.write("END\n")
    ## also copy any parameter files (ligand / ncAA mol2, frcmod, lib) so drMD can still find them
    for paramFile in os.listdir(inputDir):
        if paramFile.endswith((".mol2", ".frcmod", ".lib", ".prepi")):
            shutil.copy(p.join(inputDir, paramFile), cleanDir)
    batchConfig["pathInfo"]["originalInputDir"] = inputDir
    batchConfig["pathInfo"]["inputDir"] = cleanDir
    return batchConfig
#########################################################################
def build_membrane_system(pdbFile: FilePath,
                          outDir: DirectoryPath,
                          config: Dict,
                          keepHeteroAtoms: bool = False) -> Tuple[FilePath, Tuple[float, float, float]]:
    """
    Packs a lipid bilayer, water and ions around a protein using PACKMOL-Memgen

    Args:
        pdbFile (FilePath): protonated protein (+ ligand) PDB file
        outDir (DirectoryPath): directory where prep files are written
        config (Dict): run config containing membraneInfo
        keepHeteroAtoms (bool): keep ligands / ions in the input PDB during orientation

    Returns:
        packedPdb (FilePath): protein + bilayer + water + ions, AMBER (Lipid21) residue names
        boxDims (Tuple[float,float,float]): periodic box lengths in Angstrom
    """
    membraneInfo: Dict = config["membraneInfo"]
    memDir: DirectoryPath = p.join(outDir, "MEMBRANE")
    os.makedirs(memDir, exist_ok=True)

    packedPdb: FilePath = p.join(memDir, "MEMBRANE_packed.pdb")
    boxJson: FilePath = p.join(memDir, "MEMBRANE_box.json")
    ## skip if already built (e.g. on resume)
    if p.isfile(packedPdb) and p.isfile(boxJson):
        drLogger.log_info("Membrane system already built: skipping PACKMOL-Memgen", True)
        with open(boxJson) as f:
            boxDims = tuple(json.load(f)["boxDims"])
        return packedPdb, boxDims

    _make_sure_amberhome_is_set()

    ## copy the input pdb to the membrane directory
    memInputPdb: FilePath = p.join(memDir, "MEMBRANE_input.pdb")
    _write_clean_input_pdb(pdbFile, memInputPdb)

    ## orient the protein in the membrane ourselves (MemPrO / MEMEMBED), so failures are caught
    ## and the fully protonated drMD structure is kept; PACKMOL-Memgen then treats it as pre-oriented
    orientation: str = membraneInfo["orientation"]
    if orientation in ["auto", "mempro", "memembed"]:
        memInputPdb, orientation = orient_protein(memInputPdb, memDir, orientation)
        membraneInfo["orientationUsed"] = orientation

    ## construct the PACKMOL-Memgen command
    memgenLog: FilePath = p.join(memDir, "packmol-memgen.log")
    command: List[str] = ["packmol-memgen",
                          "--pdb", p.basename(memInputPdb),
                          "--lipids", membraneInfo["lipids"],
                          "--ratio", str(membraneInfo["lipidRatio"]),
                          "--dist", str(membraneInfo["padding"]),
                          "--dist_wat", str(membraneInfo["waterThickness"]),
                          "--salt_c", membraneInfo["cation"],
                          "--salt_a", membraneInfo["anion"],
                          "--notprotonate", "--nottrim",
                          "--noprogress", "--overwrite",
                          "--log", p.basename(memgenLog)]
    ## PACKMOL-Memgen only adds neutralising counter-ions here; bulk salt is added by tleap afterwards
    ## (n pairs = concentration x n waters / 55.5), so the requested concentration is reached on top of
    ## neutralisation whatever the lipid / protein charge
    if membraneInfo.get("boxSizeXY"):
        command += ["--distxy_fix", str(membraneInfo["boxSizeXY"])]
    ## orientation is always handled above (or supplied by the user), so PACKMOL-Memgen never re-orients
    command += ["--preoriented"]
    if keepHeteroAtoms:
        command += ["--keepligs"]

    drLogger.log_info(f"Building {membraneInfo['lipids']} bilayer with PACKMOL-Memgen (orientation: {orientation})...", True)
    _run_in_dir(command, memDir, "Membrane building with PACKMOL-Memgen", memgenLog)

    ## find PACKMOL-Memgen output
    memgenOutputs = [f for f in glob.glob(p.join(memDir, "bilayer_*.pdb")) if not f.endswith("_FORCED")]
    if len(memgenOutputs) == 0:
        drSplash.print_prep_failed(FileNotFoundError("PACKMOL-Memgen did not write a bilayer_*.pdb file"),
                                   "Membrane building with PACKMOL-Memgen", memgenLog)
    shutil.copy(memgenOutputs[0], packedPdb)

    ## box dimensions from the packing restraints (equivalent to PACKMOL-Memgen --tight_box)
    boxDims = get_box_dimensions(p.join(memDir, "packmol.inp"), packedPdb)
    nWaters: int = count_waters(packedPdb)
    nSaltPairs: int = int(round(membraneInfo["saltConcentration"] * nWaters / 55.5))
    with open(boxJson, "w") as f:
        json.dump({"boxDims": list(boxDims), "nWaters": nWaters, "nSaltPairs": nSaltPairs,
                   "membraneInfo": membraneInfo, "orientationUsed": orientation,
                   "nLipidResidues": count_lipid_residues(packedPdb)}, f, indent=2)
    drLogger.log_info(f"Membrane system built: box {boxDims[0]:.1f} x {boxDims[1]:.1f} x {boxDims[2]:.1f} A", True)
    return packedPdb, boxDims
#########################################################################
def get_salt_pairs(outDir: DirectoryPath) -> int:
    """Number of cation/anion pairs to add for the requested salt concentration (stored at build time)"""
    with open(p.join(outDir, "MEMBRANE", "MEMBRANE_box.json")) as f:
        return int(json.load(f).get("nSaltPairs", 0))
#########################################################################
def count_waters(pdbFile: FilePath) -> int:
    """Number of water molecules (one oxygen each) in a PDB file"""
    n = 0
    with open(pdbFile) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() in ["WAT", "HOH"] and line[12:16].strip() == "O":
                n += 1
    return n
#########################################################################
def orient_protein(inputPdb: FilePath, memDir: DirectoryPath, orientation: str) -> Tuple[FilePath, str]:
    """
    Orients the protein so the membrane normal is along z and the membrane centre is at z = 0.

    orientation:
        "auto"     -> MemPrO if installed, otherwise MEMEMBED; if the first fails, try the other
        "mempro"   -> MemPrO only
        "memembed" -> MEMEMBED only

    Returns:
        orientedPdb (FilePath): full (protonated) protein in the membrane frame
        methodUsed (str): "mempro" or "memembed"
    """
    memproAvailable: bool = shutil.which("mempro") is not None
    memembedAvailable: bool = _find_memembed() is not None
    if orientation == "auto":
        methods = (["mempro"] if memproAvailable else []) + (["memembed"] if memembedAvailable else [])
    else:
        methods = [orientation]

    errors = []
    for method in methods:
        try:
            if method == "mempro":
                orientedPdb = run_mempro_orientation(inputPdb, memDir)
            else:
                orientedPdb = run_memembed_orientation(inputPdb, memDir)
            drLogger.log_info(f"Protein oriented in the membrane with {method}", True)
            return orientedPdb, method
        except Exception as e:
            drLogger.log_info(f"WARNING: membrane orientation with {method} failed: {e}", True)
            errors.append(f"{method}: {e}")

    advice = ("Could not orient the protein in the membrane.\n"
              + "\n".join(errors) +
              "\nEither:\n"
              "  - install MemPrO (uv sync --extra membrane) and use orientation: mempro (or auto), or\n"
              "  - supply an oriented structure (e.g. from the OPM database) and use orientation: preoriented")
    drSplash.print_prep_failed(RuntimeError(advice), "Membrane orientation", None)
#########################################################################
def run_memembed_orientation(inputPdb: FilePath, memDir: DirectoryPath) -> FilePath:
    """
    Orients the protein with MEMEMBED (bundled with PACKMOL-Memgen / AmberTools).
    MEMEMBED is run on a heavy-atom copy with standard residue names; the full protonated
    protein is then superposed onto the MEMEMBED result.
    """
    memembedExe = _find_memembed()
    if memembedExe is None:
        raise FileNotFoundError("MEMEMBED executable not found")
    heavyPdb: FilePath = p.join(memDir, "MEMEMBED_input_heavy.pdb")
    embedPdb: FilePath = p.join(memDir, "MEMEMBED_oriented_heavy.pdb")
    memembedLog: FilePath = p.join(memDir, "memembed.log")
    _write_heavy_atom_pdb(inputPdb, heavyPdb)
    if p.isfile(embedPdb):
        os.remove(embedPdb)
    with open(memembedLog, "w") as log:
        result = subprocess.run([memembedExe, "-s", "3", "-n", "in", "-o", p.basename(embedPdb), p.basename(heavyPdb)],
                                cwd=memDir, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0 or not p.isfile(embedPdb):
        raise RuntimeError(f"MEMEMBED failed (exit code {result.returncode}), see {memembedLog}")
    orientedPdb: FilePath = p.join(memDir, "MEMBRANE_input_oriented.pdb")
    _superpose_onto(inputPdb, embedPdb, orientedPdb)
    return orientedPdb
#########################################################################
def _find_memembed() -> Optional[str]:
    return drDependencies.find_memembed()
#########################################################################
def _atom_key(line: str) -> Tuple[str, str, str]:
    """chain, residue number + insertion code, atom name"""
    return (line[21], line[22:27].strip(), line[12:16].strip())
#########################################################################
def _write_heavy_atom_pdb(inPdb: FilePath, outPdb: FilePath) -> None:
    """Heavy atoms only, AMBER protonation-state residue names mapped back to standard names"""
    standardNames = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "CYX": "CYS", "CYM": "CYS",
                     "ASH": "ASP", "GLH": "GLU", "LYN": "LYS"}
    aminoAcids = drListInitiator.get_amino_acid_residue_names()
    with open(inPdb) as fin, open(outPdb, "w") as fout:
        for line in fin:
            if line.startswith("TER"):
                fout.write("TER\n")
                continue
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resName = line[17:20].strip()
            atomName = line[12:16].strip()
            if resName not in aminoAcids or resName in ["ACE", "NME", "NHE"]:
                continue
            if atomName.startswith("H") or (atomName[0].isdigit() and atomName[1:2] == "H"):
                continue
            resName = standardNames.get(resName, resName)
            ## standard PDB atom-name columns: names shorter than 4 characters start in column 14 (" CA ")
            ## MEMEMBED reads atom names by fixed column and silently scores nothing (then crashes) otherwise
            pdbAtomName = f" {atomName:<3}" if len(atomName) < 4 else atomName[:4]
            fout.write(f"ATOM  {line[6:11]} {pdbAtomName} {resName:>3}{line[20:54]}  1.00  0.00           {atomName[0]}\n")
        fout.write("END\n")
#########################################################################
def _read_coords(pdbFile: FilePath) -> Dict[Tuple[str, str, str], List[float]]:
    coords = {}
    with open(pdbFile) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                coords.setdefault(_atom_key(line), [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return coords
#########################################################################
def _superpose_onto(mobilePdb: FilePath, referencePdb: FilePath, outPdb: FilePath) -> None:
    """
    Rigid-body (Kabsch) superposition of every atom in mobilePdb onto the matching atoms in referencePdb.
    Used to carry the full protonated protein into the frame of an orientation tool's output.
    """
    import numpy as np
    mobile = _read_coords(mobilePdb)
    reference = _read_coords(referencePdb)
    shared = [k for k in reference if k in mobile]
    if len(shared) < 3:
        raise RuntimeError("could not match atoms between the input and the oriented structure")
    P = np.array([mobile[k] for k in shared])
    Q = np.array([reference[k] for k in shared])
    pCentre, qCentre = P.mean(0), Q.mean(0)
    H = (P - pCentre).T @ (Q - qCentre)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    rotation = Vt.T @ np.diag([1, 1, d]) @ U.T
    fitted = (P - pCentre) @ rotation.T + qCentre
    rmsd = float(np.sqrt(((fitted - Q) ** 2).sum(1).mean()))
    if rmsd > 0.5:
        raise RuntimeError(f"orientation output does not match the input structure (RMSD {rmsd:.2f} A)")
    with open(mobilePdb) as fin, open(outPdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                xyz = (np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]) - pCentre) @ rotation.T + qCentre
                line = f"{line[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}"
            fout.write(line)
    drLogger.log_info(f"Superposed input onto oriented structure: {len(shared)} atoms, RMSD {rmsd:.3f} A")
#########################################################################
def run_mempro_orientation(inputPdb: FilePath, memDir: DirectoryPath) -> FilePath:
    """
    Orients a protein in a membrane with MemPrO (optional dependency, Python >= 3.11)
    Returns a protein-only PDB with the membrane normal along z and the membrane centre at z = 0
    """
    memproExe: Optional[str] = shutil.which("mempro")
    if memproExe is None:
        raise FileNotFoundError("MemPrO is not installed (uv sync --extra membrane)")
    memproDir: DirectoryPath = p.join(memDir, "MemPrO")
    memproLog: FilePath = p.join(memDir, "mempro.log")
    env = dict(os.environ)
    env.setdefault("JAX_PLATFORMS", "cpu")
    command = [memproExe, "-f", p.abspath(inputPdb), "-o", memproDir + os.sep, "-ng", "16", "-ni", "150"]
    drLogger.log_info("Orienting protein with MemPrO...", True)
    with open(memproLog, "w") as log:
        result = subprocess.run(command, cwd=memDir, stdout=log, stderr=subprocess.STDOUT, env=env)
    orientedPdb: FilePath = p.join(memproDir, "Rank_1", "oriented_rank_1.pdb")
    if result.returncode != 0 or not p.isfile(orientedPdb):
        raise RuntimeError(f"MemPrO failed (exit code {result.returncode}), see {memproLog}")
    ## carry the full input structure (all atoms, original records) into MemPrO's membrane frame
    ## (MemPrO's dummy membrane beads, resname DUM, are not matched and so are dropped)
    outPdb: FilePath = p.join(memDir, "MEMBRANE_input_oriented.pdb")
    _superpose_onto(inputPdb, orientedPdb, outPdb)
    return outPdb
#########################################################################
def get_box_dimensions(packmolInp: FilePath, packedPdb: FilePath, tolerance: float = 2.0) -> Tuple[float, float, float]:
    """
    Reads the "inside box" restraints from the PACKMOL input to get the box size.
    PACKMOL keeps molecules `tolerance` apart, so this is added to each length to avoid
    clashes across periodic boundaries. Falls back to coordinate extents.
    """
    mins = [None, None, None]
    maxs = [None, None, None]
    if p.isfile(packmolInp):
        with open(packmolInp) as f:
            for line in f:
                parts = line.split()
                if len(parts) == 8 and parts[0] == "inside" and parts[1] == "box":
                    values = [float(x) for x in parts[2:]]
                    for i in range(3):
                        mins[i] = values[i] if mins[i] is None else min(mins[i], values[i])
                        maxs[i] = values[i + 3] if maxs[i] is None else max(maxs[i], values[i + 3])
    if None in mins:
        pdbDf = pdbUtils.pdb2df(packedPdb)
        mins = [pdbDf[c].min() for c in ["X", "Y", "Z"]]
        maxs = [pdbDf[c].max() for c in ["X", "Y", "Z"]]
    return tuple(float(maxs[i] - mins[i] + tolerance) for i in range(3))
#########################################################################
def count_lipid_residues(pdbFile: FilePath) -> Dict[str, int]:
    """Counts lipid fragment residues (Lipid21 naming) in a PDB file"""
    lipidNames = drListInitiator.get_lipid_residue_names()
    counts: Dict[str, int] = {}
    seen = set()
    with open(pdbFile) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resName = line[17:21].strip()
            if resName not in lipidNames:
                continue
            key = (line[21], line[22:27], resName)
            if key in seen:
                continue
            seen.add(key)
            counts[resName] = counts.get(resName, 0) + 1
    return counts
#########################################################################
def _write_clean_input_pdb(inPdb: FilePath, outPdb: FilePath) -> None:
    """Copies ATOM/HETATM/TER records only (PACKMOL-Memgen is strict about PDB format)"""
    with open(inPdb) as fin, open(outPdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM", "TER")):
                fout.write(line)
        fout.write("END\n")
#########################################################################
def _make_sure_amberhome_is_set() -> None:
    """PACKMOL-Memgen needs AMBERHOME and the tools environment on PATH (see drDependencies)"""
    drDependencies.set_up_tool_paths()
#########################################################################
def _run_in_dir(command: List[str], workDir: DirectoryPath, stepName: str,
                logFile: Optional[FilePath], env: Optional[dict] = None) -> None:
    drLogger.log_info(f"Running: {stepName}")
    drLogger.log_info(" ".join(command))
    try:
        result = subprocess.run(command, cwd=workDir, capture_output=True, text=True,
                                check=True, env=env or os.environ)
    except Exception as errorMessage:
        stdout = getattr(errorMessage, "stdout", "") or ""
        stderr = getattr(errorMessage, "stderr", "") or ""
        with open(p.join(workDir, f"{stepName.split()[0].lower()}_failed.txt"), "w") as f:
            f.write(stdout + "\n" + stderr)
        drSplash.print_prep_failed(errorMessage, stepName, logFile)
    if result.stdout:
        drLogger.log_info(f"Command output:\n{result.stdout[-4000:]}")
    drLogger.log_info(f"{stepName} Completed Successfully")
