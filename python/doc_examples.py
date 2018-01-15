from rdkit import Chem
from rdkit.Chem import AllChem
import rdkitjson

mols = []

smi='c1c(C=CC)cccc1O/C=C\\[C@H]([NH3+])Cl'
m = Chem.MolFromSmiles(smi)
m.SetProp("_Name","example 1")
m.SetIntProp("prop1",1)
m.SetDoubleProp("prop2",3.14)
m.SetProp("prop3","foo")
mols.append(m)

m = Chem.AddHs(Chem.MolFromSmiles('O[C@H]([35Cl])F'))
m.SetProp("_Name","example 2")
AllChem.ComputeGasteigerCharges(m)
AllChem.Compute2DCoords(m)
AllChem.EmbedMolecule(m,clearConfs=False)
mols.append(m)

m = Chem.AddHs(Chem.MolFromSmiles('O[C@H](Cl)/C=C/C'))
m.SetProp("_Name","example 3")
AllChem.ComputeGasteigerCharges(m)
AllChem.Compute2DCoords(m)
AllChem.EmbedMolecule(m,clearConfs=False)
mols.append(m)

for m in mols:
    mjson = rdkitjson.molstojson([m],includePartialCharges=True)
    print(mjson.replace('}, ','},\n ').replace('], "','], \n"'))
    newm = rdkitjson.jsontomols(mjson)[0]
    assert(Chem.MolToSmiles(newm)==Chem.MolToSmiles(m))
