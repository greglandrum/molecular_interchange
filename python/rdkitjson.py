#
# Copyright (C) 2017 Greg Landrum
# All Rights Reserved
#
import json
from rdkit import rdBase,Chem
from rdkit.Chem import AllChem

def moltojson(m,includePartialCharges=True):
    from collections import OrderedDict
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
    res["atomDefaults"] = obj_type(Z=6,chiral=False,impHs=0,chg=0,stereo="unspecified",nrad=0)
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

    obj = obj_type(toolkit="RDKit",toolkit_version=rdBase.rdkitVersion,format_version=1)
    obj["aromaticAtoms"] = [x.GetIdx() for x in m.GetAtoms() if x.GetIsAromatic()]
    obj["aromaticBonds"] = [x.GetIdx() for x in m.GetBonds() if x.GetIsAromatic()]
    obj["bondRings"] = [list(x) for x in m.GetRingInfo().BondRings()]
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
                m.GetBondWithIdx(idx).SetIsAromatic(True)
            # bondRings is also there, but we can't do anything with that from Python.
            break

    m.UpdatePropertyCache(strict=strict)
    return m

if(__name__=='__main__'):
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles('c1ccccc1O/C=C\\[C@H]([NH3+])Cl')
    AllChem.Compute2DCoords(m)
    m.SetProp("_Name","example 1")
    mjson = moltojson(m)
    print(mjson)
    newm = jsontomol(mjson)
    print(Chem.MolToSmiles(newm))
    print(Chem.MolToMolBlock(newm))
