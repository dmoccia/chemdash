
from rdkit.Chem import Descriptors, Mol, rdMolDescriptors, Crippen, GraphDescriptors, Lipinski
from chemdash.SA_Score import sascorer
from chemdash.NP_Score import npscorer


round_digs = 2

def get_logp_(mol: Mol) -> float:
    return round(Crippen.MolLogP(mol), round_digs)


def get_qed_(mol: Mol) -> float:
    return round(Descriptors.qed(mol), round_digs)


def get_tpsa_(mol: Mol) -> float:
    return round(rdMolDescriptors.CalcTPSA(mol), round_digs)


def get_bertz_(mol: Mol) -> float:
    return round(GraphDescriptors.BertzCT(mol), round_digs)


def get_mol_weight_(mol: Mol) -> float:
    return round(Descriptors.MolWt(mol), round_digs)


def get_num_h_donors_(mol: Mol) -> int:
    return Lipinski.NumHDonors(mol)


def get_num_h_acceptors_(mol: Mol) -> int:
    return Lipinski.NumHAcceptors(mol)


def get_num_rotatable_bonds_(mol: Mol) -> int:
    return Lipinski.NumRotatableBonds(mol)


def get_num_rings_(mol: Mol) -> int:
    return rdMolDescriptors.CalcNumRings(mol)


def get_num_aromatic_rings_(mol: Mol) -> int:
    return rdMolDescriptors.CalcNumAromaticRings(mol)


def get_num_heavy_atoms_(mol: Mol) -> int:
    """Returns the number of heavy atoms in the molecule"""
    return Lipinski.HeavyAtomCount(mol)

def get_fraction_csp3_(mol: Mol) -> float:
    """Returns fraction of SP3 Carbon atoms"""
    return round(rdMolDescriptors.CalcFractionCSP3(mol), round_digs)


def get_synth_accessibility_(mol: Mol) -> float:
    """
    Computes RDKit's Synthetic Accessibility score
    """
    return round(sascorer.calculateScore(mol), round_digs)


def get_nat_prod_likeness_(mol: Mol) -> float:
    """
    Computes RDKit's Natural Product-likeness score
    """
    return round(npscorer.scoreMol(mol), round_digs)


def get_mpo_score_(mol: Mol) -> float:
    """
    Computes MPO score
    """
    #TODO implement MPO Score...

    return 1.0




def get_descriptor_type():
    return      {'MolecularWeight': 'continuous',
                 'HeavyAtomCount': 'discrete',
                 'HydrogenBondAcceptors': 'discrete',
                 'HydrogenBondDonors': 'discrete',
                 'RotatableBonds': 'discrete',
                 'RingCount': 'discrete',
                 'FractionCSP3': 'continuous',
                 'LogP': 'continuous',
                 'TPSA': 'continuous',
                 'NumAromaticRings': 'discrete',
                 'QED': 'continuous',
                 'SyntheticAccessibility': 'continuous',
                 'NaturalProductLikeness': 'continuous',
                 #'MPOScore': 'continuous',
                 }