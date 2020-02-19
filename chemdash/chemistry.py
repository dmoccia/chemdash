import pandas as pd
from typing import List

from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Chem import Draw


def df_canonicalize_from_smiles(df, smiles_col: str, include_stereocenters=True)->pd.Series:
    """
    Canonicalize the SMILES strings with RDKit.
    Args:
        df: dataframe
        smiles_col: column name in df
        include_stereocenters: whether to keep the stereochemical information in the canonical SMILES string

    Returns:
        Canonicalized SMILES string, None if the molecule is invalid.
    """

    def df_canonicalize_from_smiles_inner(row):
        mol = Chem.MolFromSmiles(row[smiles_col])
        if mol is not None:
            return Chem.MolToSmiles(mol, isomericSmiles=include_stereocenters)
        else:
            return None

    return df.apply(lambda x: pd.Series(df_canonicalize_from_smiles_inner(x)), axis=1)



def df_smiles_to_rdkit_mol(df, smiles_col: str) ->pd.Series:
    """
    Converts a SMILES string to a RDKit molecule.

    Args:
        df: dataframe
        smiles_col: column name in df

    Returns:
        RDKit Mol, None if the SMILES string is invalid
    """

    def df_smiles_to_rdkit_mol_inner(row):
        mol = Chem.MolFromSmiles(row[smiles_col])
        #  Sanitization check (detects invalid valence)
        if mol is not None:
            try:
                Chem.SanitizeMol(mol)
            except ValueError:
                return None
        return mol

    return df.apply(lambda x: pd.Series(df_smiles_to_rdkit_mol_inner(x)), axis=1)


def df_get_image_file_url(df, mol_col: str, id_col: str)->pd.Series:
    """
    Takes an rdkti molecule and generates a markdown url and saves file locally

    Args:
        df: dataframe
        mol_columns: mol column name in df

    Returns:
        markdown string for file location (used in the Dash Datatable)
    """

    def df_get_image_file_url_inner(row):
        mol = row[mol_col]
        id = row[id_col]
        img_file_url = None
        if mol is not None:
            try:
                img_file_name = 'assets/mol_images/' + id + '.png'
                Chem.Draw.MolToFile(mol, img_file_name, size=(200, 200))
                img_file_url = '![' + id + '](' + img_file_name + ')'
            except ValueError:
                return None
        return img_file_url

    return df.apply(lambda x: pd.Series(df_get_image_file_url_inner(x)), axis=1)

def df_create_image_file(df, smiles_col: str, id_col: str):

    """
    As the above df_get_image_file_url but simply generates the image file, not currently used
    but could be useful for on the fly generation of structures
    """

    def df_create_image_file_inner(row):
        smiles = row[smiles_col]
        id = row[id_col]
        img_file_url = None
        if smiles is not None:
            try:
                img_file_name = 'assets/mol_images/' + id + '.png'
                Chem.Draw.MolToFile(Chem.MolFromSmiles(smiles), img_file_name, size=(200, 200))
                img_file_url = '![' + id + '](' + img_file_name + ')'
            except ValueError:
                return None

    df.apply(lambda x: pd.Series(df_create_image_file_inner(x)), axis=1)



def df_add_maccs_fps(df, mol_column: str) -> (pd.DataFrame, List, str):

    num_bits = 167
    prefix_bit_cols = 'MACCS_'

    bit_columns = []
    for i in range(num_bits):
        bit_columns.append(prefix_bit_cols + str(i))

    def split(word):
        return list(word)

    def bit_vector_expand(row):
        x = row[mol_column]
        fp = MACCSkeys.GenMACCSKeys(x)
        vector_split = split(fp)
        return vector_split

    lambda_func = lambda x: pd.Series(bit_vector_expand(x))

    df[bit_columns] = df.apply(lambda_func, axis=1)

    return df, bit_columns, prefix_bit_cols


def df_add_ecfp_1024_4_fps(df, mol_column: str) -> (pd.DataFrame, List, str):


    num_bits = 1024
    prefix_bit_cols = 'Morgan_4_1024_'

    bit_columns = []
    for i in range(num_bits):
        bit_columns.append(prefix_bit_cols + str(i))

    def split(word):
        return list(word)

    def bit_vector_expand(row):
        x = row[mol_column]
        fp = AllChem.GetMorganFingerprintAsBitVect(x, 2, 1024)
        vector_split = split(fp)
        return vector_split

    lambda_func = lambda x: pd.Series(bit_vector_expand(x))

    df[bit_columns] = df.apply(lambda_func, axis=1)

    return df, bit_columns, prefix_bit_cols

