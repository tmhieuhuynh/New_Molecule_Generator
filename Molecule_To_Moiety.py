import numpy as np
import pandas as pd
import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import DrawingOptions
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=True

def generalize_unknown(smi):
    try:
        smi=smi.split('.')
        for i in range (0,len(smi)):
            smi[i]=smi[i].replace('*]','x]')
            while 'x]' in smi[i]:
                smi[i]=smi[i][:smi[i][:smi[i].index('x')].rindex('[')]+'*'+smi[i][smi[i].index('x')+2:]

            mol=Chem.MolFromSmiles(smi[i])
            if len(mol.GetRingInfo().AtomRings())>0:
                mol=Chem.DeleteSubstructs(mol,Chem.MolFromSmiles('*'))
                smi[i]=Chem.MolToSmiles(mol)

        return smi
    except: return [smi]

def break_down(smi):
    mol=Chem.MolFromSmiles(smi)

    atoms={}
    for i in mol.GetAtoms():
        atoms[i.GetIdx()]=i.GetSymbol()

    rings=[]
    for i in mol.GetRingInfo().AtomRings():
        rings.append(i)

    cut=[]
    for bond in mol.GetBonds():
        if atoms[bond.GetBeginAtomIdx()]=='C' and atoms[bond.GetEndAtomIdx()]=='C' and bond.GetBondTypeAsDouble()==1:
            status=0
            for i in rings:
                if bond.GetBeginAtomIdx() in i and bond.GetEndAtomIdx() in i:
                    status=1
            if status==0:
                cut.append(bond.GetIdx())
    
    if cut!=[]:
        frags=Chem.FragmentOnBonds(mol, cut)
    else:
        frags=mol
    return generalize_unknown(Chem.MolToSmiles(frags))

def ring_cut(smi):
    mol=Chem.MolFromSmiles(smi)

    atoms={}
    for i in mol.GetAtoms():
        atoms[i.GetIdx()]=i.GetSymbol()

    rings=[]
    for i in mol.GetRingInfo().AtomRings():
        rings.append(i)

    cut=[]
    for bond in mol.GetBonds():
        status=0
        for i in rings:
            if bond.GetBeginAtomIdx() in i and bond.GetEndAtomIdx() in i:
                status=1
                break
            if (bond.GetBeginAtomIdx() in i and bond.GetEndAtomIdx() not in i)\
                or (bond.GetBeginAtomIdx() not in i and bond.GetEndAtomIdx() in i):
                status=2
        if status==2 and bond.GetBondTypeAsDouble()==1:
            cut.append(bond.GetIdx())

    if cut!=[]:
        frags=Chem.FragmentOnBonds(mol, cut)
    else:
        frags=mol

    return generalize_unknown(Chem.MolToSmiles(frags))

f=open(sys.argv[1])
smi=f.readlines()
f.close()
smi=[x[:-1] for x in smi]

moiety=[]
for i in smi:
    for x in break_down(i):
        if x not in moiety:
            moiety.append(x)

group=[]
ring=[]
ring_tail=[]

for i in moiety:
    mol=Chem.MolFromSmiles(i)
    if len(mol.GetRingInfo().AtomRings())>0:
        for x in ring_cut(i):
            sub=Chem.MolFromSmiles(x)
            if len(sub.GetRingInfo().AtomRings())>0:
                if x not in ring:
                    ring.append(x)
            elif x not in ring_tail:
                ring_tail.append(x)
    else:
        group.append(i)

data=[]
for i in smi:
    sub=[]
    for x in break_down(i):
        mol=Chem.MolFromSmiles(x)
        if len(mol.GetRingInfo().AtomRings())>0:
            sub+=ring_cut(x)
        else:
            sub.append(x)

    part=[i]
    for x in group+ring+ring_tail:
        part.append(sub.count(x))
    data.append(part)

data=pd.DataFrame(data,columns=['Smile']+group+ring+ring_tail)
data.to_csv(sys.argv[2],index=False)

        