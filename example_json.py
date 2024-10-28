import json
import os


def is_json_file(file_path):
    """
    Checks if a file contains correctly formatted JSON.

    Args:
        file_path (str): The path to the file to be checked.

    Returns:
        bool: True if the file contains correctly formatted JSON, False otherwise.

    Example:
        is_valid = is_json_file('file.json')
    """
    try:
        with open(file_path, 'r') as file:
            json.load(file)
        return True
    except (ValueError, IOError):
        return False


def combine_json_files(input_dir, output_file):
    """
    Combines multiple JSON files into a single file, ensuring unique top-level keys and checking for conflicting data.

    Args:
        input_dir (str): The path to the directory containing the input JSON files.
        output_file (str): The path to the output file for the combined JSON data.

    Raises:
        ValueError: If conflicting data is found for a key, indicating a manual check is required.

    Example:
        combine_json_files('input_directory', 'output_file.json')
    """
    combined_data = {}
    # Iterate through each file in the input directory
    for file_name in os.listdir(input_dir):
        # build the full file path from the input_dir and file_name
        file_path = os.path.join(input_dir, file_name)
        # only iterate through files that are correctly formatted json!
        if is_json_file(file_path):
            with open(file_path, 'r') as file:
                # get the json data, we assume the data is in the format {inchi-key: {results}} as opposed to a list of entries
                data = json.load(file)  
                # we need to check to see if the base files have a duplicate key (InChI-Key)
                for key in data:
                    # if we've seen the key before
                    if key in combined_data:
                        # check to see if it's got the same data in it.
                        # if not, raise an error, we need to validate the input
                        if combined_data[key] != data[key]:
                            raise ValueError(f"Conflicting data found for key '{key}'. Manual check required.")
                    else:
                        combined_data[key] = data[key]
        else:
            print(f"The input file {file_path} is incorrectly formatted json")

    # Write the combined data to a single output file
    with open(output_file, 'w') as outfile:
        json.dump(combined_data, outfile, indent=4)


# Example usage
input_directory = 'path_to_input_directory'  # Replace with the actual input directory path
file_path = os.path.join(input_directory, "file_name.json")
if is_json_file(file_path):
    print(f"The file {file_path} contains correctly formatted JSON.")
else:
    print(f"The file {file_path} does not contain correctly formatted JSON.")

output_file_path = 'path_to_output_file.json'  # Replace with the desired output file path
combine_json_files(input_directory, output_file_path)
