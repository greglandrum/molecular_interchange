from rdkit import Chem
from rdkit.Chem import AllChem
import rdkitjson

mols = []

smi='c1c(C=CC)cccc1O/C=C\\[C@H]([NH3+])Cl'
m = Chem.MolFromSmiles(smi)
m.SetProp("_Name","example 1")
mols.append(m)

m = Chem.AddHs(Chem.MolFromSmiles('O[C@H](Cl)F'))
m.SetProp("_Name","example 2")
AllChem.ComputeGasteigerCharges(m)
AllChem.Compute2DCoords(m)
AllChem.EmbedMolecule(m,clearConfs=False)
mols.append(m)

for m in mols:
    mjson = rdkitjson.moltojson(m)
    print(mjson.replace('}, ','},\n ').replace('], "','], \n"'))
    newm = rdkitjson.jsontomol(mjson)
    assert(Chem.MolToSmiles(newm)==Chem.MolToSmiles(m))
