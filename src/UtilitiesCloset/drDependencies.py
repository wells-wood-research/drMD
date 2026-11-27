"""
drDependencies
==============
Single place where drMD knows which optional *features* need which dependencies.

drMD's Python dependencies are managed with uv (pyproject.toml + uv.lock):
    uv sync                          -> core drMD
    uv sync --extra <feature>        -> core drMD + an optional feature (e.g. membrane)

Some programs are not on PyPI (AmberTools: tleap, antechamber, parmchk2, packmol-memgen;
OpenBabel: obabel). These live in a small conda environment (drMD_tools.yaml, default name
"drMD-tools"). drMD finds that environment automatically and puts it on PATH for every
subprocess it launches, so it never has to be "activated" alongside the uv environment.

Search order for the tools environment:
    1. $DRMD_TOOLS_PREFIX
    2. $AMBERHOME
    3. tleap already on PATH
    4. a conda / mamba / micromamba environment called "drMD-tools"
       (listed in ~/.conda/environments.txt or found in a standard envs directory)

To add a new feature with its own dependencies:
    1. add an extra to [project.optional-dependencies] in pyproject.toml and run "uv lock"
    2. register the feature in FEATURES below (python modules, executables, and a
       "required" function that says when a config needs it)
"""
## BASIC PYTHON LIBRARIES
import os
import sys
import glob
import shutil
import importlib.util
from os import path as p

## CLEAN CODE
from typing import Callable, Dict, List, Optional

TOOLS_ENV_NAME: str = "drMD-tools"
TOOLS_ENV_FILE: str = "drMD_tools.yaml"

#########################################################################
## FEATURE REGISTRY
#########################################################################
def _always(config: dict) -> bool:
    return True

def _uses_ligands(config: dict) -> bool:
    return bool(config.get("ligandInfo"))

def _uses_membrane(config: dict) -> bool:
    return bool(config.get("membraneInfo"))

def _uses_cuda(config: dict) -> bool:
    return str(config.get("hardwareInfo", {}).get("platform", "")).upper() == "CUDA"

## Each feature:
##   extra       - name of the uv / pip extra that provides its Python packages (None = core)
##   python      - {import name: PyPI name} that must be importable
##   executables - programs that must be found (PATH, uv env or drMD-tools env)
##   required    - function(config) -> True if this config needs the feature
FEATURES: Dict[str, dict] = {
    "core": {
        "description": "standard drMD simulations",
        "extra": None,
        "python": {"openmm": "openmm", "pdbUtils": "pdbutils", "mdtraj": "mdtraj",
                   "MDAnalysis": "mdanalysis", "pandas": "pandas", "yaml": "pyyaml",
                   "argpass": "argpass", "inflect": "inflect", "jinja2": "jinja2",
                   "sklearn": "scikit-learn", "propka": "propka", "pdb2pqr": "pdb2pqr"},
        "executables": ["pdb2pqr", "tleap"],
        "required": _always,
    },
    "ligands": {
        "description": "small-molecule ligands (ligandInfo)",
        "extra": None,
        "python": {},
        "executables": ["obabel", "antechamber", "parmchk2"],
        "required": _uses_ligands,
    },
    "membrane": {
        "description": "membrane-protein simulations (membraneInfo)",
        "extra": "membrane",
        "python": {},                      ## MemPrO is checked below: MEMEMBED is a fallback for orientation "auto"
        "executables": ["packmol-memgen"],
        "required": _uses_membrane,
        "check": "_check_membrane_orientation",
    },
    "cuda": {
        "description": "OpenMM CUDA platform",
        "extra": "cuda12  (or cuda13, to match your NVIDIA driver)",
        "python": {},
        "executables": [],
        "required": _uses_cuda,
        "check": "_check_cuda_platform",
    },
}

#########################################################################
## TOOL DISCOVERY
#########################################################################
def _conda_env_candidates(envName: str) -> List[str]:
    candidates: List[str] = []
    envsTxt = p.join(p.expanduser("~"), ".conda", "environments.txt")
    if p.isfile(envsTxt):
        with open(envsTxt) as f:
            candidates += [line.strip() for line in f if line.strip()]
    roots = [os.environ.get(var) for var in ("CONDA_PREFIX", "MAMBA_ROOT_PREFIX", "CONDA_ROOT")]
    roots += [p.join(p.expanduser("~"), d) for d in ("miniconda3", "miniforge3", "mambaforge",
                                                      "anaconda3", "micromamba", ".local/share/mamba")]
    for root in filter(None, roots):
        candidates.append(p.join(root, "envs", envName))
        ## CONDA_PREFIX may itself be an env inside <root>/envs
        candidates.append(p.join(p.dirname(root), envName))
    return [c for c in candidates if p.basename(c.rstrip("/\\")) == envName and p.isdir(c)]

def _bin_dir(prefix: str) -> str:
    return p.join(prefix, "Library", "bin") if os.name == "nt" else p.join(prefix, "bin")

def find_tools_prefix() -> Optional[str]:
    """Locates the environment that provides AmberTools / OpenBabel (see module docstring)"""
    for var in ("DRMD_TOOLS_PREFIX", "AMBERHOME"):
        prefix = os.environ.get(var)
        if prefix and p.isdir(prefix):
            return p.abspath(prefix)
    tleap = shutil.which("tleap")
    if tleap:
        return p.dirname(p.dirname(p.realpath(tleap)))
    for candidate in _conda_env_candidates(TOOLS_ENV_NAME):
        if p.isfile(p.join(_bin_dir(candidate), "tleap")):
            return p.abspath(candidate)
    return None

def set_up_tool_paths() -> Optional[str]:
    """
    Makes the uv environment and the drMD-tools environment visible to every subprocess:
        PATH      = <python env bin> : <tools env bin> : existing PATH
        AMBERHOME = <tools env prefix>   (if not already set)
    Safe to call more than once. Returns the tools prefix (or None if not found)
    """
    pathParts = os.environ.get("PATH", "").split(os.pathsep)
    newParts: List[str] = [p.dirname(sys.executable)]       ## pdb2pqr, mempro, ... from the uv env
    toolsPrefix = find_tools_prefix()
    if toolsPrefix:
        newParts.append(_bin_dir(toolsPrefix))
        if p.isdir(p.join(toolsPrefix, "dat", "leap")):
            os.environ.setdefault("AMBERHOME", toolsPrefix)
    os.environ["PATH"] = os.pathsep.join(dict.fromkeys(newParts + pathParts))
    return toolsPrefix

def find_tools_site_packages(packageName: str) -> Optional[str]:
    """Directory of a Python package that lives inside the tools environment (e.g. packmol_memgen)"""
    toolsPrefix = find_tools_prefix()
    if toolsPrefix is None:
        return None
    patterns = [p.join(toolsPrefix, "lib", "python3*", "site-packages", packageName),
                p.join(toolsPrefix, "Lib", "site-packages", packageName)]
    for pattern in patterns:
        hits = sorted(glob.glob(pattern))
        if hits:
            return hits[-1]
    return None

#########################################################################
## CHECKS
#########################################################################
def _module_available(moduleName: str) -> bool:
    try:
        return importlib.util.find_spec(moduleName) is not None
    except (ImportError, ValueError):
        return False

def _check_cuda_platform() -> List[str]:
    try:
        import openmm
        names = [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())]
    except Exception as e:
        return [f"could not query OpenMM platforms ({e})"]
    return [] if "CUDA" in names else ["OpenMM has no CUDA platform (available: " + ", ".join(names) + ")"]

def find_memembed() -> Optional[str]:
    """MEMEMBED ships inside PACKMOL-Memgen (tools environment)"""
    exe = shutil.which("memembed")
    if exe:
        return exe
    packmolMemgenDir = find_tools_site_packages("packmol_memgen")
    if packmolMemgenDir:
        candidate = p.join(packmolMemgenDir, "lib", "memembed", "bin", "memembed")
        if p.isfile(candidate):
            return candidate
    return None

def _check_membrane_orientation(config: Optional[dict] = None) -> List[str]:
    """Checks that the requested orientation method (membraneInfo.orientation) is available"""
    orientation = "auto"
    if isinstance(config, dict) and isinstance(config.get("membraneInfo"), dict):
        orientation = config["membraneInfo"].get("orientation", "auto")
    hasMempro = shutil.which("mempro") is not None
    hasMemembed = find_memembed() is not None
    if orientation == "mempro" and not hasMempro:
        return ["orientation 'mempro' requested but MemPrO is not installed (uv sync --extra membrane)"]
    if orientation == "memembed" and not hasMemembed:
        return ["orientation 'memembed' requested but MEMEMBED was not found (comes with PACKMOL-Memgen in drMD-tools)"]
    if orientation == "auto" and not (hasMempro or hasMemembed):
        return ["no orientation method found: install MemPrO (uv sync --extra membrane) or PACKMOL-Memgen/MEMEMBED"]
    return []

def check_feature(featureName: str, config: Optional[dict] = None) -> List[str]:
    """Returns a list of human-readable problems (empty list = feature is usable)"""
    feature = FEATURES[featureName]
    problems: List[str] = []
    for module, pypiName in feature.get("python", {}).items():
        if not _module_available(module):
            problems.append(f"Python package '{pypiName}' is not installed")
    for exe in feature.get("executables", []):
        if shutil.which(exe) is None:
            problems.append(f"program '{exe}' was not found")
    extraCheck = feature.get("check")
    if extraCheck:
        checkFunction = globals()[extraCheck]
        problems += checkFunction(config) if checkFunction.__code__.co_argcount else checkFunction()
    return problems

def install_hint(featureName: str) -> str:
    feature = FEATURES[featureName]
    hints = []
    if feature.get("python") or feature.get("extra"):
        extra = feature.get("extra")
        hints.append("uv sync" + (f" --extra {extra.split()[0]}" if extra else "") +
                     (f"    # extra: {extra}" if extra and " " in extra else ""))
    if feature.get("executables"):
        hints.append(f"micromamba env create -f {TOOLS_ENV_FILE}    # or conda / mamba; creates '{TOOLS_ENV_NAME}'")
        hints.append(f"(or point drMD at an existing AmberTools/OpenBabel environment: export DRMD_TOOLS_PREFIX=/path/to/env)")
    return "\n        ".join(hints)

def required_features(config: dict) -> List[str]:
    return [name for name, feature in FEATURES.items() if feature["required"](config)]

def check_dependencies(config: dict) -> List[str]:
    """
    Checks every feature the config needs. Returns a list of formatted error blocks
    (empty = all good). Call set_up_tool_paths() first.
    """
    messages: List[str] = []
    for name in required_features(config):
        problems = check_feature(name, config)
        if problems:
            messages.append(f"Feature '{name}' ({FEATURES[name]['description']}) is missing dependencies:\n" +
                            "".join(f"    - {problem}\n" for problem in problems) +
                            f"    To fix:\n        {install_hint(name)}")
    return messages

def report() -> str:
    """Status of every registered feature, for 'python drDependencies.py'"""
    toolsPrefix = set_up_tool_paths()
    lines = [f"Python environment : {sys.prefix}",
             f"Tools environment  : {toolsPrefix or 'NOT FOUND (see drMD_tools.yaml)'}", ""]
    for name, feature in FEATURES.items():
        problems = check_feature(name)
        status = "OK" if not problems else "missing: " + "; ".join(problems)
        lines.append(f"  {name:<10} {feature['description']:<45} {status}")
    return "\n".join(lines)

if __name__ == "__main__":
    print(report())
