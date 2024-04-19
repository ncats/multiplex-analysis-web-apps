# Import relevant libraries
import numpy as np
import pandas as pd


# This should be applied per image
def sample_index(index, n_samples):
    """
    Randomly samples `n_samples` indices from the given `index` array without replacement, returning the sorted result.

    Parameters:
        index (array-like): The array of indices to sample from.
        n_samples (int): The number of indices to sample.

    Returns:
        numpy.array: An array of randomly sampled indices.
    """
    return np.sort(np.random.choice(index, n_samples, replace=False))


# Main function
def main():

    # Parameters
    smallest_image_size_frac = 0.1
    num_analysis_subsets = 10

    # Junk load the training data
    df_train = pd.read_csv('data/train.csv')

    # Group the training data grouped by image
    df_train_by_image = df_train.groupby('Slide ID')
    
    # Get the desired number of samples per image
    smallest_image_size = df_train_by_image.size().min()
    num_samples_per_image = int(smallest_image_size * smallest_image_size_frac)

    # Get a sampling of df_train that samples rows from each image equally
    train_indices = df_train_by_image.apply(lambda df_image: sample_index(df_image.index, num_samples_per_image)).explode()  # goes into df_train. This is a series whose index is the image ID and the values are the indices of the samples. This should be used for fitting the UMAP model.

    # Do the same for num_analysis_subsets subsets but first ensuring that the training indices are not included in the possible analysis indices
    subset_indices_holder = [train_indices]
    df_train_only_analysis_indices_by_image = df_train.drop(train_indices).groupby('Slide ID')
    for _ in range(num_analysis_subsets):
        subset_indices_holder.append(df_train_only_analysis_indices_by_image.apply(lambda df_image: sample_index(df_image.index, num_samples_per_image)).explode())

    # Get one big dataframe with all the indices, both "train" and analysis
    df_subset_indices = pd.concat(subset_indices_holder, axis='columns', columns=['train_subset'] + [f'analysis_subset_{i}' for i in range(num_analysis_subsets)])

    # Ensure that the concatenation hasn't added any rows
    assert len(df_subset_indices) == len(train_indices), f'The concatenation of the indices has added rows, going from {len(train_indices)} to {len(df_subset_indices)}'


# Call the main function
if __name__ == '__main__':
    main()
