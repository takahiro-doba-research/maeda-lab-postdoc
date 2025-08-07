import calculationfile as cf
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
