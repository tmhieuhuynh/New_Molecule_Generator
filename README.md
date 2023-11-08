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
The user need to prepare a list of SMI(s) as presented in SMI_Samples.txt
