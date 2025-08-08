## Description
This directory contains the following contents:
- The input files (SCAFIR*n*.com) and the output files (SCAFIR*n*_EQ_list.log, SCAFIR*n*_PT_list.log, SCAFIR*n*_sim.log) from the *n*th SC-AFIR calculation
- All the setup files for separating EQs and connecting substituents
- `calculationfile` package for analyzing the output files of SC-AFIR calculations
- `jobcontroller` package for automating the energy descriptor calculations

## Usage
`calculationfile` and `jobcontroller` packages can be imported normally in python programs.
```python
import calculationfile as cf
import jobcontroller as jc
```
The `calculationfile` package provides functions that parse the output files of SC-AFIR calculations and return file objects containing sufficient information, such as coordinates and energies, for constructing reaction path networks. The `save_eqs()` function saves coordinates of all EQs in a folder named "EQ".
```python
wd = "."
eq_list_path = "td_c1_193_SCAFIR6_EQ_list.log"
eq_list = cf.read_eq_list(f"{wd}/{eq_list_path}")
eq_list.save_eqs(wd)
```
The `SeparationJob` can be instantiated, and the `separate()` function separates all EQs into fragments based on their bond connectivity matrices. All separated EQs are saved in a folder named "separated_group" and numbered sequentially during this process. The mappings between EQs, separated EQs, and separated groups are exported to `separated_EQ.csv` and `separated_group.csv`. The file `separated_group_info.csv` contains information about the bond connectivity.
```python
sj = jc.SeparationJob(wd)
eq_folder = "EQ"
sj.separate(eq_folder)
```
Using the files `SP_separation.com` and `MIN_separation.com`, which contain the calculation setups for single-point energy calculation and geometry optimization, respectively, the `optimize()` function submits multiple calculations to a server. The single-point calculation is skipped when a separated group contains only one separated EQ. The `analyze()` function produces a summary of the calculations.
```python
sp_separation_path = "SP_separation.com"
min_separation_path = "MIN_separation.com"
sj.optimize(56, 24, sp_separation_path, min_separation_path)
sj.analyze()
```
Similarly, the `ConnectionJob` can be instantiated, and the `make_conformer()` function connects substituents from the specified folders to separated EQs for all combinations.
```python
wd = "."
cj = jc.ConnectionJob(wd)

note_path = "note.com"
cj.add_note(note_path)

sub_folders = ["arene", "alkene", "PG", "backbone", "pyridone"]
cj.make_conformer(sub_folders)
``` 
The information on how the substituents should be connected is provided using the `add_note()` function before calling `make_conformer()`. This information is indicated by placing a sequence `fragment atom1 atom2` to the right of the coordinate of atom0 for both separated EQs and substituents. The fragment number links the parts of the separated EQ to be swapped with a substituent, while the three atoms (atom0, atom1, and atom2) from both the separated EQ and substituent define the planes used to determine the dihedral angle between them. When generating multiple conformers, a sequence `fragment atom1 atom2 step` is provided. This generates conformers with dihedral angles of 0°, 360°/step, 360°/step * 2, ..., 360°/step * (step - 1).
```bash
Pd    -4.626032566819    -1.815194844416     0.391845605436 0
C     -2.979682622253     0.356034409536    -0.178373879507 0
O     -4.254430313631     0.012001025782    -0.278971564617 0
O     -2.605663344565     1.495671947896    -0.322729149645 0
C     -1.952608672882    -0.803497899306     0.086379059431 0 14 16 3
H     -1.326248445610    -0.500697618989     0.946218878146 0
N     -2.721517176961    -1.994025841517     0.421125519527 0
C     -2.661704302814    -3.061564955062     1.213032532607 0
O     -3.793901979686    -3.620469961627     1.388646364695 0
C     -1.414574301650    -3.597773395740     1.846856289818 0
H     -1.011626421100    -4.417399317648     1.231097999856 0
H     -0.638383130353    -2.827688294826     1.951406738999 0
H     -1.674343221793    -4.011395910408     2.831407457874 0
C     -1.038332670898    -0.992475405233    -1.167165243957 4
H     -1.708124962268    -0.992739132318    -2.046042935092 4
C     -0.051574454371     0.174544613179    -1.302796880425 4
H      0.524020145769     0.075954573425    -2.236831646384 4
H     -0.565583909098     1.143597047486    -1.304269193680 4
H      0.669107441865     0.171657985402    -0.466406419822 4
C     -0.287877131382    -2.328618254178    -1.159404190450 4
H      0.406161254912    -2.396943060016    -0.304574229789 4
H     -0.969802001740    -3.190929369059    -1.127105957038 4
H      0.316776427347    -2.421345352354    -2.075270635976 4
```
Using the files `SP_connection.com` and `MIN_connection.com`, which contain the calculation setups for single-point energy calculation and geometry optimization, respectively, the `optimize()` function submits multiple calculations to a server. The single-point calculation is skipped when only one dihedral angle is specified. The `analyze()` function produces a summary of the calculations.
```python
sp_method_path = "SP_connection.com"
sp_option = '-n 24 -m "hu01 hu02 hu03 hu04 hu05 hu06 hu07 hu08 hu09 hu10 hu11 hu12 hu13 hu14"'
min_method_path = "MIN_connection.com"
min_option = "-p1 -nhu01,hu02,hu03,hu04,hu05,hu06,hu07,hu08,hu09,hu10,hu11,hu12,hu13,hu14"
cj.optimize(56, sp_method_path, sp_option, min_method_path, min_option)
cj.analyze()
```