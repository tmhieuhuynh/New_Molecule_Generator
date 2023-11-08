# New_Molecule_Generator
This repository aims to construct new molecules based on the old ones while maintaining the chemical space. It includes two main programs:

- Molecule_To_Moiety.py: to break the provided molecules into moieties.
- Smile_Generator.py: to construct new molecules from the moieties.

The programs are designed to process in [RDkit](https://www.rdkit.org/) environment.

## Installation
### Dependent Library Installation
- RDkit environment
```bash
conda create -c conda-forge -n my-rdkit-env rdkit
conda activate my-rdkit-env
```
- numpy
```bash
pip install numpy
```
- pandas
```bash
pip install pandas
```
- pubchemy
```bash
pip install pubchempy
```
### New_Molecule_Generator Installation
```bash
git clone https://github.com/Serendipity12345/New_Molecule_Generator
cd New_Molecule_Generator
```
## Tutorial
The user needs to prepare a list of SMI(s) as presented in SMI_Samples.smi
### Molecule_To_Moiety
This program takes the SMI list as the input. Then, it breaks all the input molecules into moieties. The result is displayed in a csv file with the statistical table of the molecules and the moieties.
```bash
python Molecule_To_Moiety.py ["directory of SMI list"] ["directory for csv result"]
```
To test the program:
```bash
python Molecule_To_Moiety.py SMI_Samples.smi Test.csv
```
### New_Molecule_Generator
Before executing this program, the user needs to manually convert the moiety result into a text file as Moiety_Samples.txt, which includes the information on moieties, occurrence quantity, and binding points.
The program takes the text file of moiety information as the input and generates new molecules. The new molecules are checked with Pubchem through pubchmey to ensure the chemical feasibility.
```bash
python New_Molecule_Generator.py ["number of looping times"] ["maximum molecular weight"] ["directory of the moiety information file"] ["directory for the SMI result files"]
```
To test the program:
```bash
python New_Molecule_Generator.py 10 150 Moiety_Samples.txt New_Molecules.txt
```
