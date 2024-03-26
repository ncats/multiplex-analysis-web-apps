def get_foundry_dataset(alias='input'):
    """Create a dataset object.
    This should be fast.
    This is consistent with Code Workspaces snippets on 3/10/24.
    """
    from foundry.transforms import Dataset
    return Dataset.get(alias)

def upload_file_to_dataset(dataset, selected_filepath='/home/user/repo/bleh.txt'):
    """Upload a local file to a dataset.
    This overwrites existing files.
    This trivially returns a string of the uploaded filename.
    This is consistent with Code Workspaces snippets on 3/10/24.
    This should be slow.
    """
    return dataset.upload_file(selected_filepath)

def upload_dir_to_dataset(dataset, path_to_dir_to_upload='../junk_files'):
    """Upload a local directory to a dataset.
    This overwrites existing files.
    This returns a dictionary where each key is the path to the local file in path_to_dir_to_upload and each value is the name of the file in the dataset on NIDAP, where the name includes separators (e.g., "/") and is equal to the local file path without the prefix path_to_dir_to_upload (plus following "/"), probably as you would expect on Amazon S3. E.g., the return value could be:
        {'../junk_files/junk-200mb-20': 'junk-200mb-20',
         '../junk_files/junk-200mb-09': 'junk-200mb-09',
         '../junk_files/junk-200mb-39': 'junk-200mb-39',
         '../junk_files/subdir/tmp1.txt': 'subdir/tmp1.txt',
         '../junk_files/subdir/tmp2.txt': 'subdir/tmp2.txt',
         '../junk_files/subdir/subdir2/junk-200mb-2': 'subdir/subdir2/junk-200mb-2',
         '../junk_files/subdir/subdir2/junk-200mb-4': 'subdir/subdir2/junk-200mb-4',
         '../junk_files/subdir/subdir2/junk-200mb-9': 'subdir/subdir2/junk-200mb-9',
         '../junk_files/subdir/subdir2/junk-200mb-0': 'subdir/subdir2/junk-200mb-0'}
    Using this function on 3/[9-10]/24, I get about 70-120 MB/s upload speed.
    Note there is at least a single-file upload limit of about 2000 MB, which is higher than I reported in an old Issue to Palantir.
    This should be slow.
    """
    return dataset.upload_directory(path_to_dir_to_upload)

def download_files_from_dataset(dataset, dataset_filter_func=lambda f: f.path.startswith("junk-200mb"), limit=15):
    """
    Download files from a dataset.
    This is what Palantir recommends per https://nidap.nih.gov/workspace/issues-app/issue/ri.issues.main.issue.4a72e98a-5233-47ef-bafd-01d9f5f0966c on 1/26/24.
    next_downloaded_files updates on each iteration through the `while` loop, which stops if this dictionary is empty. It is a dictionary with limit elements, where the keys are the names of the files (probably including "/"'s) in the NIDAP dataset and the values are the local paths to the downloaded files. The loop as written below is possibly inefficient because I believe the return of .download() are the just-downloaded files, not the next ones to be downloaded, so the loop runs one extra time. But this is what they recommend doing and may really be the best way to do it, and regardless it is likely trivially slower if at all.
    This function should likely supplant download_file_from_dataset(), get_dataset_file_object(), list_files_in_dataset(), and get_file_objects_from_dataset().
    This returns a dictionary of *all* downloaded files as described above.
    I'm pretty sure this does not overwrite already-downloaded files, so if it's been run once with the same parameters, it will be fast on subsequent runs.
    This is otherwise slow.
    """
    # limit=15 seems to be the best value for downloading 60 200MB files. It's unclear to me exactly what this limit is doing. But in this situation, I get an overall download speed of about 800 MB/s!
    # If I had to guess, without the limit keyword I believe all matching files are downloaded in a single batch (so e.g. the loop below is iterated once and is unnecessary), and with the limit keyword, the files are downloaded in batches of size limit (so each batch has limit files). I believe that each batch is downloaded sequentially, but within each batch, multiple CPUs are used to download files in the batch in parallel. So I think you want to have at least as many files in each batch (i.e., limit) as there are CPUs available to download files in parallel. It's a bit unclear why a single batch with all the files isn't fastest because I'd think it'd parallelize the file downloads efficiently, but e.g. limit=15 was faster than limit=20, which was faster than larger limits. Likewise, smaller limit values (than 15) were slower.
    specific_files = dataset.files().filter(dataset_filter_func)
    all_downloaded_files = dict()
    while True:
        next_downloaded_files = specific_files.download(limit=limit)
        all_downloaded_files.update(next_downloaded_files)
        if not next_downloaded_files:
            break
        # Process next_downloaded_files
    return all_downloaded_files


# ---- These are likely supplanted by download_files_from_dataset() ----------------------------------------------------------------
# Actually, not necessarily. For downloads/uploads, yes, but for ultimately getting file listings, no.
def get_file_objects_from_dataset(dataset):
    """Get a list of file objects in a dataset.
    This is likely supplanted by download_files_from_dataset().
    This is slow.
    """
    return list(dataset.files())
def list_files_in_dataset(dataset_file_objects):
    """Get a list of strings of the filenames in a dataset.
    This is likely supplanted by download_files_from_dataset().
    This should be fast.
    """
    return [x.path for x in dataset_file_objects]
def get_dataset_file_object(dataset_file_objects, selected_filename='sample_txt_file.txt'):
    """Get the Foundry file object associated with a file.
    This is likely supplanted by download_files_from_dataset().
    This should be fast.
    """
    dataset_file_list = list_files_in_dataset(dataset_file_objects)
    return dataset_file_objects[dataset_file_list.index(selected_filename)]
def download_file_from_dataset(dataset_file_object):
    """Download a file from a dataset to somewhere local.
    This never overwrites existing files because the local download path (/tmp/data/RANDOM-STRING/filename) is always different.
    This returns a string of the local download path.
    This is likely supplanted by download_files_from_dataset().
    This is slow.
    """
    return dataset_file_object.download()
    # dataset_files = dataset.files().download()  # <-- This downloads all files in the dataset per the Code Workspaces documentation, consistent with Code Workspaces snippets on 3/10/24
# ----------------------------------------------------------------------------------------------------------------------------------
