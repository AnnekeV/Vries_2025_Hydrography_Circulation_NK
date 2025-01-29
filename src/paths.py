
from pathlib import Path

path_parent = Path.cwd().parent.parent
figpath = path_parent.joinpath(Path('Figures'))
datapath = path_parent.joinpath(Path('Data'))
path_intermediate_files = Path.cwd().parent.joinpath("intermediate_files")

path_carra = "/Users/annek/Documents/CARRA/"

# Microcat moorings 2018-2019
fpath_GF13 = str(path_parent.joinpath("Data","Moorings","20190802_SBE_GF13")) +"/"
fpath_GF10 = str(path_parent)+"/Data/Moorings/20190612_SBE_GF10/"