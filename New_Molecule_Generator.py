import numpy as np
import pandas as pd
import sys
import pubchempy
from rdkit import Chem
from rdkit.Chem import Descriptors, rdmolops
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import DrawingOptions
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=True

all_moiety=pd.read_csv(sys.argv[3],sep='\t')

all_prob=np.array(all_moiety['Count'].tolist())/sum(all_moiety['Count'].tolist())
all_moiety=all_moiety[['Moiety','Binding_point']].values.tolist()
for i in range (0,len(all_moiety)):
    all_moiety[i][-1]=all_moiety[i][-1].split(',')

f=open(sys.argv[4])
results=f.readlines()
f.close()
print('Number of available SMI(s) in the result file: ',len(results))

for y in range (int(sys.argv[1])):
    print('Progress: ',y,'/',sys.argv[1])
    try:
        moiety_rand=np.random.choice(len(all_moiety),p=all_prob)
        new_molecule=Chem.MolFromSmiles(all_moiety[moiety_rand][0])

        x=0
        binding_points=all_moiety[moiety_rand][1].copy()
        while x<len(binding_points) and Descriptors.HeavyAtomMolWt(new_molecule)<int(sys.argv[2]):
            if np.random.randint(2)!=0:
                moiety_rand=np.random.choice(len(all_moiety),p=all_prob)
                molecule=Chem.CombineMols(new_molecule,Chem.MolFromSmiles(all_moiety[moiety_rand][0]))

                ed_molecule=Chem.EditableMol(molecule)
                DrawingOptions.includeAtomNumbers=True
                ed_molecule.AddBond(int(binding_points[x]),new_molecule.GetNumAtoms()+int(all_moiety[moiety_rand][1][0]),order=Chem.rdchem.BondType.SINGLE)
                if len(all_moiety[moiety_rand][1])!=1:
                    binding_points+=[str(new_molecule.GetNumAtoms()+int(i)) for i in all_moiety[moiety_rand][1][1:]]

                new_molecule = ed_molecule.GetMol()
            x+=1
        new_molecule = Chem.RemoveHs(new_molecule)
        new_molecule = Chem.AddHs(new_molecule)
        smiles=AllChem.MolToSmiles(new_molecule)
        compounds = pubchempy.get_compounds(smiles, namespace='smiles')
        match = compounds[0]

        if match.isomeric_smiles not in results:
            results.append(match.isomeric_smiles)
    except: True

print('Number of SMI(s) in the result file after running: ',len(results))
f=open(sys.argv[4],'w')
for t in results:
    f.write(t+'\n')
f.close()
print('Finished!')
