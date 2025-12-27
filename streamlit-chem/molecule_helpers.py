import re
import numpy as np

# lazy import helper (keeps Streamlit fast even if deps missing)
def _maybe_imports():
    mods = {}
    try:
        import pubchempy as pcp
        mods["pcp"] = pcp
    except Exception:
        mods["pcp"] = None
    try:
        import py3Dmol
        mods["py3Dmol"] = py3Dmol
    except Exception:
        mods["py3Dmol"] = None
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Draw
        mods["Chem"] = Chem; mods["AllChem"] = AllChem; mods["Draw"] = Draw
    except Exception:
        mods["Chem"] = mods["AllChem"] = mods["Draw"] = None
    try:
        from stmol import showmol
        mods["showmol"] = showmol
    except Exception:
        mods["showmol"] = None
    return mods

def handle_polymer_name(chemical_name: str):
    name_l = (chemical_name or "").lower()
    if name_l.startswith("poly") and "saccharide" not in name_l:
        monomer_name = re.sub(r"^poly", "", chemical_name, flags=re.IGNORECASE).strip().strip("()")
        mods = _maybe_imports()
        pcp, Chem = mods["pcp"], mods["Chem"]
        if not pcp or not Chem:
            return None, None, None
        monomer = get_compound_from_name(monomer_name)
        if monomer:
            mon_smiles = monomer.smiles
            if "C=C" in mon_smiles:
                repeat = mon_smiles.replace("C=C", "C-C")
                return monomer_name, monomer, f"*{repeat}*"
    return None, None, None

def get_compound_from_name(chemical_name: str):
    mods = _maybe_imports(); pcp = mods["pcp"]
    if not pcp or not chemical_name.strip(): return None
    try:
        cs = pcp.get_compounds(chemical_name, "name")
        return cs[0] if cs else None
    except Exception:
        return None

def get_compound_from_iupac_web(chemical_name: str):
    mods = _maybe_imports(); pcp, Chem = mods["pcp"], mods["Chem"]
    if not pcp or not Chem: return None
    import requests
    try:
        enc = requests.utils.quote(chemical_name)
        r = requests.get(f"https://opsin.ch.cam.ac.uk/opsin/{enc}.smi",
                         headers={"Accept": "chemical/x-smiles"}, timeout=5)
        if r.status_code == 200:
            smiles = r.text
            if smiles and Chem.MolFromSmiles(smiles):
                return pcp.get_compounds(smiles, "smiles")[0]
    except Exception:
        pass
    return None

def get_chemical_properties(smiles: str):
    mods = _maybe_imports(); pcp = mods["pcp"]
    if not pcp: return {}
    try:
        c = pcp.get_compounds(smiles, "smiles")[0]
        return {
            "Molecular Formula": f"`{c.molecular_formula}`",
            "Molecular Weight": f"{(c.molecular_weight or 0):.2f} g/mol",
            "LogP": f"{(c.xlogp or 0):.2f}" if getattr(c, "xlogp", None) is not None else "N/A",
            "H-Bond Donors": getattr(c, "h_bond_donor_count", "N/A"),
            "H-Bond Acceptors": getattr(c, "h_bond_acceptor_count", "N/A"),
        }
    except Exception:
        return {}

def generate_3d_mol(smiles: str):
    mods = _maybe_imports(); Chem, AllChem = mods["Chem"], mods["AllChem"]
    if not Chem or not AllChem: return None
    m = Chem.MolFromSmiles(smiles)
    if not m: return None
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m, randomSeed=1)
    AllChem.UFFOptimizeMolecule(m)
    return m

def smiles_to_2d_image(smiles: str, legend=""):
    mods = _maybe_imports(); Chem, AllChem, Draw = mods["Chem"], mods["AllChem"], mods["Draw"]
    if not Chem or not AllChem or not Draw: return None
    m = Chem.MolFromSmiles(smiles)
    if not m: return None
    AllChem.Compute2DCoords(m)
    return Draw.MolToImage(m, size=(400, 400), legend=legend)

def get_functional_group_library():
    return {
        "Alkene (C=C)": "[CX3]=[CX3]", "Benzene Ring": "c1ccccc1", "Hydroxyl (-OH)": "[OX2H]",
        "Carbonyl (C=O)": "[CX3]=O", "Carboxylic Acid (-COOH)": "[CX3](=O)[OX2H1]",
        "Amine (-NH2, -NHR, -NR2)": "[NX3;H2,H1,H0;!$(NC=O)]", "Ester (-COO-)": "[#6][CX3](=O)[OX2H0][#6]",
        "Halogen (F,Cl,Br,I)": "[#9,#17,#35,#53]",
    }

def add_functional_group_labels(view, mol):
    mods = _maybe_imports(); Chem = mods["Chem"]
    if not Chem: return view
    patt = {n: Chem.MolFromSmarts(s) for n, s in get_functional_group_library().items()}
    conf = mol.GetConformer()
    labeled = set()
    for name, p in patt.items():
        if not p: continue
        for match in mol.GetSubstructMatches(p):
            if all(i in labeled for i in match): continue
            coords = [conf.GetAtomPosition(i) for i in match]
            cx, cy, cz = np.mean([(p.x, p.y, p.z) for p in coords], axis=0)
            view.addLabel(name, {'position': {'x': cx, 'y': cy, 'z': cz},
                                 'backgroundColor': 'white', 'fontColor': 'black',
                                 'backgroundOpacity': 0.7, 'borderColor': 'black', 'fontSize': 12})
            labeled.update(match)
    return view