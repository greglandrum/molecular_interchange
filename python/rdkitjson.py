#
# Copyright (C) 2017 Greg Landrum
# All Rights Reserved
#
import json
from rdkit import rdBase,Chem
from rdkit.Chem import AllChem

def moltojson(m,includePartialCharges=True):
    from collections import OrderedDict
    # copy it since we're going to kekulize
    m = Chem.Mol(m)
    if includePartialCharges:
        AllChem.ComputeGasteigerCharges(m)
    Chem.Kekulize(m)
    obj = {}
    nm = "no name"
    if m.HasProp("_Name"):
        nm = m.GetProp("_Name")
    obj_type = OrderedDict
    res = obj_type()
    res['header'] = obj_type(version=10,name=nm)
    res["atomDefaults"] = obj_type(Z=6,impHs=0,chg=0,stereo="unspecified",nrad=0)
    res["bondDefaults"] = obj_type(bo=1,stereo="unspecified",stereoAtoms=[])
    res["atoms"] = []
    for i,at in enumerate(m.GetAtoms()):
        obj = obj_type(Z=at.GetAtomicNum())
        if at.GetTotalNumHs():
            obj["impHs"]=at.GetTotalNumHs()
        if at.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            chi = {Chem.ChiralType.CHI_TETRAHEDRAL_CW:"cw",Chem.ChiralType.CHI_TETRAHEDRAL_CCW:"ccw",Chem.ChiralType.CHI_OTHER:"other"}[at.GetChiralTag()]
            obj['stereo']=chi
        if at.GetFormalCharge():
            obj['chg'] = at.GetFormalCharge()
        if at.GetNumRadicalElectrons():
            obj['nRad'] = at.GetNumRadicalElectrons()
        res["atoms"].append(obj)
    res["atomProperties"] = []
    if includePartialCharges and m.GetAtomWithIdx(0).HasProp("_GasteigerCharge"):
        obj = obj_type(type="partialcharges",method="rdkit-gasteiter")
        obj["values"] = [float('%.3f'%float(x.GetProp("_GasteigerCharge"))) for x in m.GetAtoms()]

    res["bonds"] = []
    for i,bnd in enumerate(m.GetBonds()):
        bo = {Chem.BondType.SINGLE:1,Chem.BondType.DOUBLE:2,Chem.BondType.TRIPLE:3}[bnd.GetBondType()]
        obj = obj_type(atoms=[bnd.GetBeginAtomIdx(),bnd.GetEndAtomIdx()])
        if bo != 1:
            obj["bo"] = bo
        if bnd.GetStereo() in (Chem.BondStereo.STEREOE,Chem.BondStereo.STEREOZ,Chem.BondStereo.STEREOCIS,Chem.BondStereo.STEREOTRANS):
            obj['stereoAtoms']=list(bnd.GetStereoAtoms())
            if bnd.GetStereo() in (Chem.BondStereo.STEREOCIS,Chem.BondStereo.STEREOZ):
                obj['stereo'] = 'cis'
            elif bnd.GetStereo() in (Chem.BondStereo.STEREOTRANS,Chem.BondStereo.STEREOE):
                obj['stereo'] = 'trans'
        elif bnd.GetStereo() == Chem.BondStereo.STEREOANY:
            obj['stereo'] = 'either'
        res["bonds"].append(obj)

    if m.GetNumConformers():
        res["conformers"] = []
        for conf in m.GetConformers():
            dim=3
            if not conf.Is3D():
                dim=2
            obj = obj_type(dim=dim)
            coords=[]
            for i in range(m.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                coord = [pos.x,pos.y]
                if dim==3:
                    coord.append(pos.z)
                coords.append([float('%.4f'%x) for x in coord])
            obj['coords']=coords
            res["conformers"].append(obj)

    if m.GetAtomWithIdx(0).GetPDBResidueInfo() is not None:
        residues = {}
        for a in m.GetAtoms():
            ri = a.GetPDBResidueInfo()
            if ri is None:
                continue
            num = ri.GetResidueNumber()
            code = ri.GetInsertionCode()
            key = (num,code)
            if key not in residues:
                d = {}
                d['num'] = num
                d['insertioncode'] = code
                d['name'] = ri.GetResidueName()
                d['atoms'] = []
                d['atomNames'] = []
                d['serialnumbers'] = []
                # include serial number?
                residues[key] = d
            residue = residues[key]
            if ri.GetIsHeteroAtom():
                residue['containsHetatms'] = True
            residue['atoms'].append(a.GetIdx())
            residue['atomNames'].append(ri.GetName())
            residue['serialnumbers'].append(ri.GetSerialNumber())

        res['residues'] = [residues[x] for x in sorted(residues)]


    obj = obj_type(toolkit="RDKit",toolkit_version=rdBase.rdkitVersion,format_version=1)
    obj["aromaticAtoms"] = [x.GetIdx() for x in m.GetAtoms() if x.GetIsAromatic()]
    obj["aromaticBonds"] = [x.GetIdx() for x in m.GetBonds() if x.GetIsAromatic()]
    if(m.GetAtomWithIdx(0).HasProp('_CIPRank')):
        obj["cipRanks"] = [int(x.GetProp("_CIPRank")) for x in m.GetAtoms()]
    obj["cipCodes"] = [[x.GetIdx(),x.GetProp("_CIPCode")] for x in m.GetAtoms() if x.HasProp("_CIPCode")]

    obj["atomRings"] = [list(x) for x in m.GetRingInfo().AtomRings()]
    res["representations"] = [obj]

    return json.dumps(res)

def jsontomol(text,strict=True):
    obj = json.loads(text)
    if obj['header']['version'] != 10:
        raise ValueError('bad version %s'%obj['header']['version'])
    nm = obj['header']['name']
    m = Chem.RWMol()
    m.SetProp('_Name',nm)
    if 'atomDefaults' in obj:
        atomDefaults = obj['atomDefaults']
    else:
        atomDefaults={}
    if 'bondDefaults' in obj:
        bondDefaults = obj['bondDefaults']
    else:
        bondDefaults={}
    # ---------------------------------
    #      Atoms
    for entry in obj['atoms']:
        atm = Chem.Atom(entry.get('Z',atomDefaults.get('Z',6)))
        atm.SetNoImplicit(True)
        atm.SetNumExplicitHs(entry.get('impHs',atomDefaults.get('impHs',0)))
        atm.SetFormalCharge(entry.get('chg',atomDefaults.get('chg',0)))
        tags = {'unspecified':Chem.ChiralType.CHI_UNSPECIFIED,'ccw':Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
              'cw':Chem.ChiralType.CHI_TETRAHEDRAL_CW,'other':Chem.ChiralType.CHI_OTHER}
        atm.SetChiralTag(tags[entry.get('stereo',atomDefaults.get('stereo','unspecified'))])
        atm.SetNumRadicalElectrons(entry.get('nRad',atomDefaults.get('nRad',0)))
        m.AddAtom(atm)
    # ---------------------------------
    #      Atom Properties
    for entry in obj['atomProperties']:
        # we ignore these for the moment
        pass

    # ---------------------------------
    #      Bonds
    # at the moment we can't set bond stereo directly because all atoms need to be there, so hold
    # that info for a bit
    bondStereos={}
    for entry in obj['bonds']:
        bos = {1:Chem.BondType.SINGLE,2:Chem.BondType.DOUBLE,3:Chem.BondType.TRIPLE}
        bo = bos[entry.get('bo',bondDefaults.get('bo',Chem.BondType.SINGLE))]
        nbs = m.AddBond(entry['atoms'][0],entry['atoms'][1],bo)
        bnd = m.GetBondWithIdx(nbs-1)
        tags = {'cis':Chem.BondStereo.STEREOCIS,'trans':Chem.BondStereo.STEREOTRANS,
               'either':Chem.BondStereo.STEREOANY,'unspecified':Chem.BondStereo.STEREONONE}
        stereo = tags[entry.get('stereo',bondDefaults.get('stereo','unspecified'))]
        if 'stereoAtoms' in entry:
            bondStereos[bnd.GetIdx()]=(entry['stereoAtoms'],stereo)
        elif stereo in (Chem.BondStereo.STEREOCIS,Chem.BondStereo.STEREOTRANS):
            raise ValueError("bond stereo set, but stereoatoms not provided")
    for idx,(vs,stereo) in bondStereos.items():
        bnd = m.GetBondWithIdx(idx)
        bnd.SetStereoAtoms(vs[0],vs[1])
        bnd.SetStereo(stereo)

    # ---------------------------------
    #      Conformers
    for entry in obj.get('conformers',[]):
        conf = Chem.Conformer(m.GetNumAtoms())
        dim = entry.get('dim',3)
        if dim==3:
            conf.Set3D(True)
        else:
            conf.Set3D(False)
        for i in range(m.GetNumAtoms()):
            coord = entry['coords'][i]
            if dim != 3:
                coord.append(0.)
            conf.SetAtomPosition(i,Chem.rdGeometry.Point3D(coord[0],coord[1],coord[2]))
        m.AddConformer(conf,assignId=True)

    # ---------------------------------
    #      representation
    for entry in obj.get('representations'):
        if entry['toolkit'] == 'RDKit':
            if entry['format_version'] != 1:
                raise ValueError("bad format_version %s"%entry['format_version'])
            aromAtoms = entry.get('aromaticAtoms',[])
            for idx in aromAtoms:
                m.GetAtomWithIdx(idx).SetIsAromatic(True)
            aromBonds = entry.get('aromaticBonds',[])
            for idx in aromBonds:
                bnd = m.GetBondWithIdx(idx)
                bnd.SetIsAromatic(True)
                bnd.SetBondType(Chem.BondType.AROMATIC)
            atomRings = entry.get('atomRings',[])
            for ring in atomRings:
                ringBonds = []
                alist = ring+[ring[0]]
                for i in range(len(ring)):
                    ringBonds.append(m.GetBondBetweenAtoms(alist[i],alist[i+1]).GetIdx())
                m.GetRingInfo().AddRing(ring,ringBonds)
            for i,x in enumerate(entry.get('cipRanks',[])):
                m.GetAtomWithIdx(i).SetProp('_CIPRank',str(x))
            for i,x in entry.get('cipCodes',[]):
                m.GetAtomWithIdx(i).SetProp('_CIPCode',x)
                #m.GetAtomWithIdx(i).SetIntProp('_ChiralityPossible',1)
            break
    m.UpdatePropertyCache(strict=strict)
    m.SetIntProp("_StereochemDone",1)
    return m

if(__name__=='__main__'):
    from rdkit.Chem import AllChem
    if 0:
        smi='c1c(C=CC)cccc1O/C=C\\[C@H]([NH3+])Cl'
        # smi ='Cc1nnc(SCC2CS[C@@H]3[C@H](NC(=O)Cn4cnnn4)C(=O)N3C=2C(=O)[O-])s1.[Na+]'
        # smi ='Cc1nnc(SCC2CS[C@@]3(C)[C@](C)(NC(=O)Cn4cnnn4)C(=O)N3C=2C(=O)[O-])s1.[Na+]'
        # smi ='N[C@H]1[C@H]2SCC=C(N2C1=O)C([O-])=O'
        # smi = 'Cl.CCN(CCN1c2cccc3c(C)c(C)n(c32)CC1=O)CC'
        m = Chem.MolFromSmiles(smi)
        m.SetProp("_Name","example 1")
    else:
        from rdkit import RDConfig
        import os
        m = Chem.MolFromPDBFile(os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers','test_data','1CRN.pdb'))
    mjson = moltojson(m)
    print(mjson)
    newm = jsontomol(mjson)
    assert(Chem.MolToSmiles(newm)==Chem.MolToSmiles(m))
