## PIP
[] set up PIP for version 1.0.0

## tests
[] look through and see what needs to be done!

## configs
[x] frcmod --> frcmod

## error handling
[x] report what function died "__name__" dunder method

## real time updates after crash
[x] kill job upon crash

## drFirstAid
[] make 0.1 fs x 100 ps??
[x] alternate between em and 0.1 fs x 10 ps NpT // NVT (use same as sim in question??)
[x] overlap between progress decorator when firstAid is running

## drSelector
[x] "all" --> "*" as a wildcard

## drConfigTriage
[] better errors upon YAML format error
[] defaults for each simulation step?
[] better disorder handling for restraints


## ncAA support

config[miscInfo][nonCanonicalResidueNames] contains a list of ncAAs to look for. Use this to exclude from pdb triage and ligand prep steps

in pdb triage, when we flag a ncAA, check for RES_NAME.mol2 / frcmod / lib.

If not all found --> Throw error, point to drFrankenstein

If all fine, point tleap towards ncAA params

Might have problems with Protein Protonation steps [PROTON_TRANSFUSION]
Use Tleap + 



