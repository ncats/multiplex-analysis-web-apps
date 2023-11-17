import io
import os
import pandas as pd

class foundry_IO_lib:
    """
    A class for handling the input and output SDK commands for the Palantir Foundry library

    Attributes:
        dataset (obj): object created by the Palantir package The attributes and methods attached
                       to this object will differ depending on if you are operating on Foundry
                       or if you are operating off of Foundry.
        onNIDAP (bool): Boolean statement to whether or not the FIOL obj recognizes that its on 
                        Foundry (NIDAP) or not.
    """
    def __init__(self):
        """
        Initializes a FIOL (Foundry IO Lib) Object
        """

        hostN = os.environ.get('FOUNDRY_HOSTNAME', 'Not found')
        token = os.environ.get('FOUNDRY_TOKEN', 'Not found')
        if (hostN == 'Not found') | (token == 'Not found'):
            # Import SDK handling library
            from palantir.datasets import dataset
            self.dataset = dataset
            # Inform on working environment
            print('Not Operating on NIDAP')
            self.onNIDAP = False
            # Set SDK secrets
            self.setEnvVar()
        else:
            # Import SDK handling library
            from foundry.transforms import Dataset
            self.dataset = Dataset
            # Inform on working environment
            print('Operating on NIDAP')
            self.onNIDAP = True

    def setEnvVar(self):
        """
        Called when the Palantir SDK HOSTNAME and TOKEN cannot be found 
        in the enviornment variables. This function reads in stored 
        values from a txt file and sets them as the approrpiately named 
        Environmental variables. 
        """

        APIFile = os.path.dirname(__file__) + '/SDK_SECRETS.txt'
        f = open(APIFile, 'r')
        client_HOSTNAME = str(f.readline()).strip()
        client_TOKEN = str(f.readline()).strip()
        f.close()

        print('Setting SDK EnvVar')
        os.environ['PALANTIR_HOSTNAME'] = str(client_HOSTNAME)
        os.environ['PALANTIR_TOKEN'] = str(client_TOKEN)

    def load_listofFiles(self, projectPath):
        """
        Returns a dictionary of the file paths and names for each of the
        files listed in an unstructured dataset
        """

        if self.onNIDAP:
            filesDict = self.dataset.get(projectPath).files().download()
        else:
            files = self.dataset(projectPath).list_files()
            fileList = list(files)
            filePaths = [file.path for file in fileList]
            filesDict = dict(zip(filePaths, fileList))
        return filesDict

    def load_dataset(self, dataset_path, files_dict, file_path, loadCompass=False):
        """
        Load dataset from Foundry Compass or from local directory as a PANDAS dataframe(df)
        """
        if loadCompass:
            if self.onNIDAP:
                # Unstructred Dataset
                if files_dict is not None:
                    # Load the datafile into a Pandas dataframe
                    file = files_dict[file_path]

                    return pd.read_csv(file)
                # Structured Dataset
                else:
                    return self.dataset.get(file_path).read_table(format="arrow").to_pandas()
            else:
                # Unstructred Dataset
                if files_dict is not None:
                    # Load the datafile into a Pandas dataframe
                    file = files_dict[file_path]
                    rawfile = file.read().read()

                    return pd.read_csv(io.BytesIO(rawfile))
                # Structured Dataset
                else:
                    fullpath = dataset_path + file_path
                    return self.dataset(fullpath).read_pandas()
        # Local Import
        else:
            return pd.read_csv(dataset_path, sep=',')
        
    def export_results_dataset(self, df, path, filename, saveCompass=False, type = 'S', create = True):
        """
        Write a dataset and phenotype assignments to disk.

        Args:
            df (Pandas dataframe): Dataframe containing some data
            path (str): Name of the CSV file
            saveCompass (bool, optional): Whether or not we are saving to NIDAP or not
            type (string, 'S' or 'U'): (S)tructured or (U)nstructured dataset export
        """
        # Save to Compass
        if saveCompass:
            if self.onNIDAP:
                # Unstructured Dataset Export
                if type == 'U':
                    # First save it to a temp file
                    foundryObj = self.dataset.get(path)
                    fileNameFull = filename + '.csv'
                    df.to_csv(fileNameFull, index=False)
                    foundryObj.upload_file(fileNameFull)
                    os.remove(fileNameFull)

                # Structured Dataset Export
                elif type == 'S':
                    full_path = path + filename
                    ds = self.dataset(full_path, create=create)
                    ds.write_pandas(df)
            else:
                # Unstructured Dataset Export
                if type == 'U':
                    # First save it to a temp file
                    df.to_csv('temp.csv', index=False)
                    # Next open it as binary, read it, then close
                    in_file = open('temp.csv', "rb") # opening for [r]eading as [b]inary
                    csvData = in_file.read()
                    in_file.close()
                    # Finally, ship it to NIDAP
                    self.dataset(path) \
                        .file(filename + '.csv') \
                        .write(content = csvData)
                    os.remove('temp.csv')
                # Structured Dataset Export
                elif type == 'S':
                    full_path = path + filename
                    ds = self.dataset(full_path, create=create)
                    ds.write_pandas(df)
        # Local csv Export
        else:
            full_path = path + filename + '.csv'
            df.to_csv(full_path, index=False)
        print(f'Uploaded {filename}')

    def save_png_dataset(self, datafile, pngFileName, pltFig):
        """
        Save png to a NIDAP dataset
        """

        # File Name
        fileNameFull = pngFileName + '.png'
        # Save as a png in the local directory using the Matplotlib 'savefig' method
        pltFig.savefig(fileNameFull)

        if self.onNIDAP:
            foundryObj = self.dataset.get(datafile)
            foundryObj.upload_file(fileNameFull)

        else:
            in_file = open(fileNameFull, 'rb') # opening for [r]eading as [b]inary
            pngData = in_file.read()
            in_file.close()

            self.dataset(datafile) \
                .file(fileNameFull) \
                .write(content = pngData)
        
        # Once Upload is complete, delete the local file
        os.remove(fileNameFull)
        print(f'Uploaded {fileNameFull}')

    def export_file_dataset(self, dataset, filename):
        if self.onNIDAP:
            foundryObj = self.dataset.get(dataset)
            foundryObj.upload_file(f'output/{filename}')

        else:
            in_file = open(f'output/{filename}', 'rb') # opening for [r]eading as [b]inary
            file_data = in_file.read()
            in_file.close()

            self.dataset(dataset) \
                .file(filename) \
                .write(content = file_data)

