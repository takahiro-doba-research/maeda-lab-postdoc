## Description
This directory contains the following contents:
- `calculationfile` package for analyzing the output files of SC-AFIR calculation
- `jobcontroller` package for automating the energy descriptor calculation

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
Using the files `SP_separation.com` and `MIN_separation.com`, which contain the calculation setups for single-point energy calculation and geometry optimization, respectively, the optimize() function submits multiple calculations to a server. The single-point calculation is skipped when a separated group contains only one separated EQ. The `analyze()` function produces a summary of the calculations.
```python
sp_separation_path = "SP_separation.com"
min_separation_path = "MIN_separation.com"
sj.optimize(56, 24, sp_separation_path, min_separation_path)
sj.analyze()
```




```python
import jobcontroller as jc

wd = "."
cj = jc.ConnectionJob(wd)

note_path = "note.com"
cj.add_note(note_path)

sub_folders = ["arene", "alkene", "PG", "backbone", "pyridone"]
cj.make_conformer(sub_folders)

sp_method_path = "SP_connection.com"
sp_option = '-n 24 -m "hu01 hu02 hu03 hu04 hu05 hu06 hu07 hu08 hu09 hu10 hu11 hu12 hu13 hu14"'
min_method_path = "MIN_connection.com"
min_option = "-p1 -nhu01,hu02,hu03,hu04,hu05,hu06,hu07,hu08,hu09,hu10,hu11,hu12,hu13,hu14"
cj.optimize(56, sp_method_path, sp_option, min_method_path, min_option)

cj.analyze()
```