import json
import sys

from chembl_webresource_client.new_client import new_client

from rdkit import Chem
from rdkit.Chem import Descriptors


def get_chembl_id(molecule_name):
    """
    Searches for a molecule by its name/synonym and retrieves its ChEMBL ID.

    Args:
        molecule_name (str): The name or synonym of the molecule for which the ChEMBL ID is to be retrieved.

    Returns:
        str or None: The ChEMBL ID of the molecule, or None if not found.

    Example:
        chembl_id = get_chembl_id('aspirin')
    """
    # Access ChEMBL client to search for the molecule by name/synonym
    molecule = new_client.molecule
    # Perform a search and retrieve the ChEMBL ID if found
    result = molecule.search(molecule_name)
    if result:
        return result[0]['molecule_chembl_id']
    # Return None if the molecule is not found
    return None


def get_compound_details(chembl_id):
    """
    Retrieves compound details including synonyms, structure, and properties for a given ChEMBL ID.

    Args:
        chembl_id (str): The ChEMBL ID of the compound for which details are to be retrieved.

    Returns:
        dict: A dictionary containing the compound details.

    Example:
        compound_details = get_compound_details('CHEMBL123456')
    """
    # Access ChEMBL client to retrieve compound details
    molecule = new_client.molecule
    # Return the details of the compound for the given ChEMBL ID
    return molecule.get(chembl_id)


def get_activity_data(chembl_id):
    """
    Retrieves the activity data for a given compound, including target information and references from ChEMBL database.

    Args:
        chembl_id (str): The ChEMBL ID of the compound for which activity data is to be retrieved.

    Returns:
        list: A list containing the activity data for the given compound.

    Example:
        activity_data = get_activity_data('CHEMBL123456')
    """
    # Access ChEMBL client to retrieve activity data
    activities = new_client.activity
    # Return the filtered activity data for the given ChEMBL ID
    return activities.filter(molecule_chembl_id=chembl_id)


def find_target_by_chembl_id(data, chembl_id):
    """
    Find dictionaries in the data list where the targets list contains the specified chembl_id.
    
    This method is required since the target field in the component-target mapping is not 
    searchable via the API call directly.

    Parameters:
    - data: List of dictionaries to search through.
    - chembl_id: The target_chembl_id to search for.

    Returns:
    - List of dictionaries where the target_chembl_id is found.
    """
    return [
        entry for entry in data
        if any(target.get('target_chembl_id') == chembl_id for target in entry.get('targets', []))
    ]


def fill_data_model(molecule_name):
    """
    Fills the data model with details obtained from ChEMBL database for a given molecule.

    Args:
        molecule_name (str): The name of the molecule for which data model is to be filled.

    Returns:
        dict: A dictionary containing the data model for the given molecule.

    Example:
        data_model = fill_data_model("aspirin")
    """
    # Retrieve the ChEMBL ID for the molecule
    chembl_id = get_chembl_id(molecule_name)
    if chembl_id:
        references = []
        activity_details = []
        # If we find a match, retrieve activity data for the molecule
        activity_data = get_activity_data(chembl_id)
        # Extract compound details
        compound_name, inchi_key, synonyms, compound_structure, properties = process_compound_details(chembl_id)
        # Extract activity data and populate references and activity_details lists
        for activity in activity_data:
            process_activity_data(references, activity_details, activity)
        # Now we can constructing the data in the data model:
        data_model = {
            inchi_key: {
                "pref_name": compound_name,
                "inchi_key": inchi_key,
                "synonyms": synonyms,
                "structure": compound_structure,
                "properties": properties,
                "targets": activity_details,
                "references": references
            }
        }
        return data_model


def process_activity_data(references, activity_details, activity):
    """
    Processes activity data for a given compound and target, and populates the references and activity_details lists.

    Args:
        references (list): A list of references to which reference details will be appended.
        activity_details (list): A list of activity details to which activity details will be appended.
        activity (dict): A dictionary containing details of the activity for a compound and target.

    Returns:
        None

    Raises:
        SystemExit: If unexpected conditions are encountered during the data processing.

    Example:
        references = []
        activity_details = []
        activity_data = {...}  # Provide the activity data dictionary
        process_activity_data(references, activity_details, activity_data)

    """
    target_info = None
    targets = []
    # Retrieve the target information based on the target_chembl_id
    target_chembl_id = activity['target_chembl_id']
    target_details = new_client.target.filter(target_chembl_id=target_chembl_id)    
    
    # Filter and process the target details
    for target in target_details:
        print(activity['target_chembl_id'], target['target_type'])
        if target['target_type'] == "SINGLE PROTEIN":
            targets.append(target)
    # if we've found results, then we can decide how we go forward
    if targets:
        target_info = targets[0]
    elif len(targets) > 1:
        print(f"Since we're only looking for single proteins, we do not expect to see {len(targets)} targets for {target_chembl_id}!")
        print(activity)
        sys.exit(0)

    # If target_info is available and is a single protein, continue processing
    if target_info and target_info['target_type'] == "SINGLE PROTEIN":
        target_component_results = new_client.target_component.filter(tax_id=target_info['tax_id'])
        # Find matching target components based on target_chembl_id, given that the target_component class 
        # isn't searchable by target_chembl_id, we need to iterate through the results to find the specific match
        matched_targets = find_target_by_chembl_id(target_component_results, target_chembl_id)
        # We really should only find a single instance since we're looking for "SINGLE PROTEIN"
        # if we don't let's drop out of the code and explore what the results look like 
        # so we can decide what to do.
        if len(matched_targets) > 1:
            print(f"TOO MANY TARGETS for {target_chembl_id}! ({len(matched_targets)})")
            print(matched_targets)
            sys.exit(0)
        
         # If matched_targets exist, process the reference and activity details
        if matched_targets:
            target_component = matched_targets[0]
            # Handle cases where target type is not a single protein
            if target_info['target_type'] != "SINGLE PROTEIN":
                print("somthing's going wonky...")
                sys.exit(0)
            # Retrieve reference details and ensure uniqueness
            reference = new_client.document.get(activity['document_chembl_id'])
            reference_details = {
                        "title": reference['title'],
                        "authors": reference['authors'],
                        "journal": reference['journal'],
                        "year": reference['year'],
                        "doi": reference['doi']
                    }
                    # let's make sure the reference list is unique:
            if reference_details not in references:
                references.append(reference_details)
            # Populate activity details and ensure uniqueness
            activity_detail = {
                        "target_name": target_info['pref_name'],
                        "uniprot_accession": target_component['accession'],
                        "taxon": target_info['tax_id'],
                        "activity_type": activity['standard_type'],
                        "activity_value": activity['standard_value'],
                        "unit": activity['standard_units'],
                        "source": activity['src_id'],
                        "source_id": activity['document_chembl_id'],
                        'assay_type': activity['assay_type'],
                        'assay_description': activity['assay_description']
                    }
                    # and make sure the activity detail is also unique (no duplicates!)
            if activity_detail not in activity_details:
                activity_details.append(activity_detail)


def process_compound_details(chembl_id):
    """
    Retrieves and processes details of a compound from ChEMBL database using the ChEMBL ID.

    Args:
        chembl_id (str): The ChEMBL ID of the compound.

    Returns:
        tuple: A tuple containing the compound name, InChI Key, list of synonyms, compound structure (SMILES notation), and compound properties.

    Raises:
        Exception: If there is an error in retrieving or processing the compound details.

    Example:
        compound_name, inchi_key, synonyms, compound_structure, properties = process_compound_details('CHEMBL123456')
    """
    # Retrieve compound details from ChEMBL database using ChEMBL ID
    compound_details = get_compound_details(chembl_id)
    # Extract compound details
    compound_name = compound_details['pref_name']
    inchi_key = compound_details['molecule_structures']['standard_inchi_key']
    # Extract synonyms from compound details
    synonyms = []
    if 'molecule_synonyms' in compound_details:
        synonyms = [{"synonym": syn['molecule_synonym'], "source": 'ChEMBL', 
                     "type": syn['syn_type']} for syn in compound_details['molecule_synonyms']]
    compound_structure = compound_details['molecule_structures']['canonical_smiles']
    smiles = compound_details['molecule_structures']['canonical_smiles']
    mol = Chem.MolFromSmiles(smiles)  # calculate the RDKit Molecule object for generating the calculated descriptors
    # Generate compound properties
    properties = {
            "smiles": smiles,
            "inchi": compound_details['molecule_structures']['standard_inchi'],
            "molecular_weight": compound_details['molecule_properties']['full_mwt'],
            "formula": compound_details['molecule_properties']['full_molformula'],
            "h_bond_donors": compound_details['molecule_properties']['hbd'],
            "h_bond_acceptors": compound_details['molecule_properties']['hba'],
            "rotatable_bonds": compound_details['molecule_properties']['rtb'],
            "solubility": compound_details['molecule_properties']['alogp'],
            "logP": Descriptors.MolLogP(mol),            
            "pKa": {'acidic': compound_details['molecule_properties']['cx_most_apka'], 
                    'basic':  compound_details['molecule_properties']['cx_most_bpka']},
            "tpsa": round(Descriptors.TPSA(mol), 2),
        }    
    # return the results:
    return compound_name, inchi_key, synonyms, compound_structure, properties


def main():  # Example usage
    names_list = ["aspirin", "paracetamol", "omeprazole"]  # define the list of molecules.
    model_data = {}  # TODO, modify this to be a dictionary
    # process the given list of molecule names
    for molecule_name in names_list:
        print(f"Processing {molecule_name}")
        # gather the data from ChEMBL for the current molecule_name
        molecule_data = fill_data_model(molecule_name)
        # add the results to the final list only if we've not seen the result before
        # Can you think of a more elegant way of doing this given that we have a unique id of InChI-Key?
        for inchikey, annotations in molecule_data.items():
            if inchikey not in model_data:
                model_data[inchikey] = annotations

    # write the results to a json file
    with open("compounds.json", "w") as outfile:
        json.dump(model_data, outfile)


if __name__ == '__main__':
    main()
