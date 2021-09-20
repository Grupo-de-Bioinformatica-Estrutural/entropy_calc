# Entropy calculator

## Installing
```
pip3 install -r requirements
```

## Running
To run in a dynamics, you need to create a dihedrals file, where each line is a dihedral from the small molecule. Example:
```
1,2,3,4
2,3,4,5
```
**You can also see an dihedral input example file in the repository (dihedrals_input_example.txt)**

### Script command
```
python3 entropy_ddh.py -xtc {xtx filename} -tpr {tpr filename} -i {dihedrals input file}
```

If you want to understand the parameters (or change default):
```
python3 entropy_ddh.py -h
```

### Output
| Filename | Description |
| :---: | :---: |
| entropy_for_frames.txt | Contains total system entropy in each frame |
| Entropies_DDH.png | Contains total system entropy in a plot representation |
