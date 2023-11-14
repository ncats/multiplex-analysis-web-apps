def get_foundry_dataset(alias='input'):
    """Create a dataset object.
    This should be fast.
    """
    import foundry.transforms
    return foundry.transforms.Dataset.get(alias)

def get_file_objects_from_dataset(dataset):
    """Get a list of file objects in a dataset.
    This is slow.
    """
    return list(dataset.files())

def list_files_in_dataset(dataset_file_objects):
    """Get a list of strings of the filenames in a dataset.
    This should be fast.
    """
    return [x.path for x in dataset_file_objects]

def get_dataset_file_object(dataset_file_objects, selected_filename='sample_txt_file.txt'):
    """Get the Foundry file object associated with a file.
    This should be fast.
    """
    dataset_file_list = list_files_in_dataset(dataset_file_objects)
    return dataset_file_objects[dataset_file_list.index(selected_filename)]

def download_file_from_dataset(dataset_file_object):
    """Download a file from a dataset to somewhere local.
    This never overwrites existing files because the local download path (/tmp/data/RANDOM-STRING/filename) is always different.
    This returns a string of the local download path.
    This is slow.
    """
    return dataset_file_object.download()

def upload_file_to_dataset(dataset, selected_filepath='/home/user/repo/bleh.txt'):
    """Upload a local file to a dataset.
    This overwrites existing files.
    This trivially returns a string of the uploaded filename.
    This should be slow.
    """
    return dataset.upload_file(selected_filepath)
