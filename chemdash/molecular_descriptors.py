from chemdash import descriptors as _descriptors
import pandas as pd
from rdkit import Chem


# will look to move this to a file in the future but should work for now


class MolecularDescriptors(object):
    """
    Takes a dataframe with RDKit Structures and generates descriptors
    based on input list.  There is a default list as well.  Returns a dataframe with
    descriptors
    """

    default_descriptors = ['MolecularWeight', 'HydrogenBondAcceptors', 'HydrogenBondDonors',
                           'RotatableBonds', 'LogP', 'TPSA', 'HeavyAtomCount', 'RingCount',
                           'FractionCSP3', 'NumAromaticRings', 'QED', 'SyntheticAccessibility',
                           'NaturalProductLikeness']

    def __init__(self,

                 descriptors=default_descriptors,
                 mol_field='structure',
                 smiles_field=None,
                 id_field=None):

        """Extracts data from input as Pandas data frame"""
        if not isinstance(descriptors, list):
            raise ValueError("descriptors must be a list.")

        # get descriptor dictionary...
        self.descriptors_dict= self.get_descriptor_dict()
        #print(self.descriptors_dict)

        if all(i in self.descriptors_dict.keys() for i in descriptors):
            self.descriptors = descriptors

        else:
            raise ValueError(
                "Descriptor list contains unknown descriptor type, please refer to supported descriptors in documentation ")
        if (mol_field is None and smiles_field is None):
            raise ValueError("must define a mol_field or smiles_field which calculations will be based on")
        self.mol_field = mol_field
        self.smiles_field = smiles_field
        if id_field is None:
            self.id_field = smiles_field
        else:
            self.id_field = id_field



    def generate_descriptors(self, df):


        def inner_lambda(row):
            calcs = []

            if self.mol_field is not None:
                mol = row[self.mol_field]

            else:
                print('No Mol')
                mol = Chem.MolFromSmiles(row[self.smiles_field])
            if mol is not None:
                for elem in self.descriptors:
                    # Lookup the function in your dictionary
                    func_to_call = self.descriptors_dict[elem]
                    # And call it!
                    calcs.append(func_to_call(mol))

            return calcs

        lambda_func = lambda x: pd.Series(inner_lambda(x))
        df[self.descriptors] = df.apply(lambda_func, axis=1)

        return df


    def get_descriptor_dict(self):
        # mapping input names to functions, probably move this to a file, don't want to change function names at this time
        return  {'MolecularWeight': _descriptors.get_mol_weight_,
                                 'HeavyAtomCount': _descriptors.get_num_heavy_atoms_,
                                 'HydrogenBondAcceptors': _descriptors.get_num_h_acceptors_,
                                 'HydrogenBondDonors': _descriptors.get_num_h_donors_,
                                 'RotatableBonds': _descriptors.get_num_rotatable_bonds_,
                                 'RingCount': _descriptors.get_num_rings_,
                                 'FractionCSP3': _descriptors.get_fraction_csp3_,
                                 'LogP': _descriptors.get_logp_,
                                 'TPSA': _descriptors.get_tpsa_,
                                 'NumAromaticRings': _descriptors.get_num_aromatic_rings_,
                                 'QED': _descriptors.get_qed_,
                                 'SyntheticAccessibility': _descriptors.get_synth_accessibility_,
                                 'NaturalProductLikeness': _descriptors.get_nat_prod_likeness_,
                                 #MPOScore': _descriptors.get_mpo_score_,
                                 }

