# Ligand Chain/Residue Mismatch in Distance Restraints

## Symptoms

The workflow was run with a minimisation configuration containing distance
restraints. Preparation and tLeap completed successfully, but the simulation
failed immediately when the first restraint was constructed:

```text
ValueError: Expected exactly two atom indices for a distance restraint.
```

The failure originates in `Surgery/drRestraints.py`, where a distance
restraint requires its selection to resolve to exactly two atoms.

The automatic methods writer also failed afterwards with:

```text
Error writing my methods section: object of type 'int' has no len()
```

That is a separate reporting bug and is not the cause of the simulation
failure.

## Diagnosis

The input PDB and the pre-solvation merged PDB use these ligand identifiers:

```text
PLM: chain B, residue 1
FAD: chain C, residue 1
```

The restraint YAML consequently selected, among others:

```text
A:171:HG  -- B:1:C2       # CYS -- PLM
B:1:O1    -- C:1:N38      # PLM -- FAD
```

The actual solvated reference PDB used by the simulation contains:

```text
PLM C2/O1/O2: chain A, residue 183
FAD N38:      chain B, residue 1
CYS HG/SG:    chain A, residue 171
ARG NH1/NE/NH2: chain A, residue 170
```

Therefore the first selection's `B:1:C2` atom does not exist in the
reference PDB. `A:171:HG` is found, but the PLM atom is not, so the selector
returns one index instead of two.

The generated preparation files are retained in `fixture/` so the mismatch
can be reproduced without rerunning tLeap:

- `fixture/input.pdb` is the workflow input.
- `fixture/whole_pre_solvation.pdb` is the merged PDB before solvation.
- `fixture/solvated_reference.pdb` is the PDB used for restraint selection.
- `fixture/fapd_00011_0_0_1_0.yaml` contains the failing restraints.
- `fixture/prep.log` documents successful tLeap preparation.
- `fixture/workflow.log` documents the restraint failure.

## Root Cause

After tLeap, `drPrep.py` calls `drFixer.reset_chains_residues()` to restore
chain and residue identifiers in the solvated PDB. The ligand implementation
in `UtilitiesCloset/drFixer.py` currently does this:

```python
for (inputChain, inputChainDf), (templateChain, templateChainDf) in zip(
    inputLigandsDf.groupby("CHAIN_ID"),
    templateLigandsDf.groupby("CHAIN_ID")
):
    for (inputRes, inputResDf), (templateRes, templateResDf) in zip(
        inputChainDf.groupby("RES_ID"),
        templateChainDf.groupby("RES_ID")
    ):
        outputDf.loc[inputDf["RES_ID"] == inputRes, "CHAIN_ID"] = templateChain
        outputDf.loc[inputDf["RES_ID"] == inputRes, "RES_ID"] = templateRes
```

This is unsafe for two reasons:

1. It pairs ligand chains by `zip()` order rather than by ligand identity.
   Group ordering in the tLeap output and template is not a reliable mapping.
2. The assignment mask matches only `RES_ID`, not the original chain and
   residue name. FAD and PLM both use residue ID `1`, so they can be confused,
   and unrelated residues with the same ID can also be modified.

## Proposed Solution

Replace positional chain pairing with ligand-name-based mapping. For each
configured ligand name:

1. Collect unique `(CHAIN_ID, RES_ID)` residues for that ligand in the input
   and template dataframes.
2. Match residues with the same `RES_NAME`.
3. For each matched residue, update only rows matching all of:
   `RES_NAME`, original `CHAIN_ID`, and original `RES_ID`.
4. Validate that the number of input and template copies is equal. Raise a
   descriptive error if a ligand cannot be mapped unambiguously.

The essential assignment mask should be equivalent to:

```python
mask = (
    (outputDf["RES_NAME"] == ligandName)
    & (outputDf["CHAIN_ID"] == inputChain)
    & (outputDf["RES_ID"] == inputRes)
)
outputDf.loc[mask, "CHAIN_ID"] = templateChain
outputDf.loc[mask, "RES_ID"] = templateRes
```

For this fixture, the expected restored identifiers are:

```text
PLM: chain B, residue 1
FAD: chain C, residue 1
```

The same chain-aware masking principle should also be applied to
`reset_chain_residues_protein()`, which currently assigns protein chains and
residue IDs using `RES_ID` alone. That is a separate latent collision bug for
multi-chain proteins or repeated residue numbers.

## Regression Coverage

Add a regression test using the files in `fixture/` that verifies:

- PLM remains chain B, residue 1 after identifier restoration.
- FAD remains chain C, residue 1 after identifier restoration.
- Protein CYS `HG` and `SG` remain chain A, residue 171.
- All five configured distance restraints resolve to exactly two atoms when
  selected from the corrected reference PDB.
- No protein rows are modified while restoring ligand identifiers.

The methods writer's integer `RES_ID` failure should be covered separately by
converting integer identifiers to strings before calling code that expects a
list-like value, or by making `identifier_list_to_str()` handle integers.
