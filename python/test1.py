import os.path
import gzip,unittest
import rdkitjson
from rdkit import Chem,RDConfig

class TestCase(unittest.TestCase):
    @classmethod
    def readShortFile(cls):
        with gzip.open(os.path.join(RDConfig.RDBaseDir,"Regress","Data","mols.1000.txt.gz")) as inf:
            for i,line in enumerate(inf):
                nm,smi = line.split()
                print(i,smi)
                m = Chem.MolFromSmiles(smi)
                m.SetProp("_Name",nm)
                #print(Chem.MolToSmiles(m,True))
                yield m
    @classmethod
    def readLongFile(cls):
        with gzip.open(os.path.join(RDConfig.RDBaseDir,"Regress","Data","znp.50k.smi.gz")) as inf:
            for i,line in enumerate(inf):
                smi,nm = line.split()
                print(i,smi)
                m = Chem.MolFromSmiles(smi)
                m.SetProp("_Name",nm)
                #print(Chem.MolToSmiles(m,True))
                yield m
    def test01(self):
        for m in self.readShortFile():
            mjson = rdkitjson.moltojson(m)
            newm = rdkitjson.jsontomol(mjson)
            #Chem.SanitizeMol(newm)
            self.assertEqual(Chem.MolToSmiles(m,isomericSmiles=True),Chem.MolToSmiles(newm,isomericSmiles=True))

    def test02(self):
        if not doLong:
            return
        for m in self.readLongFile():
            mjson = rdkitjson.moltojson(m)
            newm = rdkitjson.jsontomol(mjson)
            #Chem.SanitizeMol(newm)
            self.assertEqual(Chem.MolToSmiles(m,isomericSmiles=True),Chem.MolToSmiles(newm,isomericSmiles=True))

    def testzbo(self):
        m = Chem.RWMol(Chem.MolFromSmiles('C=O.[Fe]'))
        m.AddBond(1,2,Chem.BondType.ZERO)
        mjson = rdkitjson.moltojson(m)
        print(mjson)
        newm = rdkitjson.jsontomol(mjson)
        self.assertTrue(newm.GetBondBetweenAtoms(1,2))
        self.assertEqual(newm.GetBondBetweenAtoms(1,2).GetBondType(),Chem.BondType.ZERO)


if __name__ == '__main__':  # pragma: nocover
  import argparse
  import sys
  global doLong
  parser = argparse.ArgumentParser()
  parser.add_argument('-l', default=False, action='store_true', dest='doLong')
  args = parser.parse_args()
  doLong = args.doLong

  # Remove the -l flag if present so that it isn't interpreted by unittest.main()
  if '-l' in sys.argv:
    sys.argv.remove('-l')
  unittest.main()
