import utils
import new_phenotyping_lib

class TIMECellInteraction:
    '''
    Instantiation of this class mainly loads Consolidata_data.txt into a Pandas dataframe (or reads in a simulated one in the case of simulated data) and performs some preprocessing on it

    It will create a pickle file (initial_data.pkl) of the read-in and preprocessed data, unless the file already exists, in which case this step is skipped
    '''

    def __init__(self, dataset_obj, project_dir, allow_compound_species, nslices=1, thickness_new=4, n_neighs=6, radius_instead_of_knn=True, simulate_data=False, refine_plotting_map_using_mapping_dict=False, flatten=True, use_analytical_significance=True, decimate_top_species=False, **kwargs):  # incorporating new dataset object

        # Import relevant module
        import os

        # Extract attributes, some not actually needed, for compatibility with original way the SIP library worked before there was a dataset object
        coord_units_in_microns, min_coord_spacing, input_data_filename, _, mapping_dict, phenotype_identification_tsv_file = dataset_obj.extract_original_parameters()

        # Set local variables
        thickness = thickness_new / coord_units_in_microns  # we want "thickness" to be in the same units as the coordinates in the input file. Essentially, we're just making sure the units match from the get-go
        dataset_name = 'slices_{}x{}'.format(nslices, thickness_new)
        pickle_dir = os.path.join('.', 'output', 'checkpoints')
        if not os.path.exists(pickle_dir):
            os.makedirs(pickle_dir)
        webpage_dir = os.path.join('.', 'output', 'images')
        if not os.path.exists(webpage_dir):
            os.makedirs(webpage_dir)

        # Set attributes
        self.dataset_name = dataset_name
        self.project_dir = project_dir
        self.input_data_filename = input_data_filename  # if this is a dataframe then that's because dataset_formats was called with an already existing dataframe, not a string path
        self.webpage_dir = webpage_dir
        self.mapping_dict = mapping_dict
        self.nslices = nslices
        self.thickness = thickness
        self.n_neighs = n_neighs
        self.radius_instead_of_knn = radius_instead_of_knn
        self.coord_units_in_microns = coord_units_in_microns
        self.min_coord_spacing = min_coord_spacing
        self.thickness_new = thickness_new
        self.use_analytical_significance = use_analytical_significance

        # These next block isn't technically needed but it helps to set these here to help for linting purposes
        # These are set in this method but not saved in the traditional way (instead, using make_pickle_dict())
        self.pickle_dir = pickle_dir  # directory for storing the processed data, i.e., pickle files
        self.unique_species = []
        self.doubling_type = None
        self.unique_slides = []
        self.is_real_data = None
        self.compound_species_allowed = None
        self.plotting_map = []
        self.num_colors = None

        # These are set in other functions in this class but not saved in the traditional way (instead, using make_pickle_dict())
        self.data_by_slide = []
        self.dr = None
        self.k_max = None
        self.min_nvalid_centers = None

        # Assign local variables that aren't the same as those inputted in order to save them later using make_pickle_dict()
        is_real_data = not simulate_data
        compound_species_allowed = allow_compound_species

        # Constant
        pickle_file = 'initial_data.pkl'

        # If the pickle file doesn't exist...
        if not os.path.exists(os.path.join(pickle_dir, pickle_file)):

            # If requesting simulated data...
            if simulate_data:

                # ...generate simulated data
                midpoints = kwargs['midpoints']
                max_real_area = kwargs['max_real_area']
                mult = kwargs['mult']
                self.doubling_type = kwargs['doubling_type']
                self.data = get_simulated_data(kwargs['doubling_type'], midpoints, max_real_area, min_coord_spacing, mult)

            else:

                # ...otherwise, read in the data from the CSV file
                doubling_type = None
                self.data = dataset_obj.data

            # Phenotype as desired            
            if phenotype_identification_tsv_file is not None:
                pheno_method = 'Custom'
            else:
                if allow_compound_species:
                    pheno_method = 'Species'
                else:
                    pheno_method = 'Marker'
            self.data, self.phenotypes, species_int_to_pheno_name = new_phenotyping_lib.apply_phenotyping(self.data, pheno_method, phenotype_identification_tsv_file, species_int_colname='Species int')

            #### ADD IN species_int_to_pheno_name TO THE FUNCTION BELOW!!!!

            # Get the plotting map, number of colors, unique species, and the unique slides
            plotting_map, num_colors, unique_species, unique_slides, self.data = get_dataframe_info(self.data, self.phenotypes, species_int_to_pheno_name=species_int_to_pheno_name)

            # Save the case/slide/ROI "table of contents", including the ROI widths and heights
            self.df_data_by_roi = generate_case_slide_roi_contents(self.data)

            # Print and save the recommended radius as 10% of the minimum average ROI size
            self.recommended_radius_ = calculate_recommended_radius(self.df_data_by_roi)

            # Reduce dataset here in order to test independence of significance methodology of non-involved species
            if decimate_top_species:
                print('WARNING: Decimating top species by 90%!!')
                top_species = plotting_map[0][0]
                data_sample = self.data[self.data['Species int'] == top_species].sample(frac=0.9, random_state=42)  # didn't yet test addition of ", random_state=42"
                self.data = self.data.drop(data_sample.index).reset_index()
                plotting_map[0][2] = plotting_map[0][2] - len(data_sample)

            # Set some properties that haven't already been defined
            self.num_colors = num_colors
            self.plotting_map = plotting_map
            self.unique_species = unique_species

            # Add to the used-to-be-called-contents dataframe the rest of the ROI-specific data
            if flatten:
                df_data_by_roi = self.flatten_roi_plotting_data()

            # Save the data to a pickle file
            data = self.data
            self.make_pickle_dict(['pickle_dir', 'is_real_data', 'compound_species_allowed', 'doubling_type', 'data', 'phenotypes', 'plotting_map', 'num_colors', 'unique_species', 'unique_slides', 'df_data_by_roi'], locals(), pickle_file)

        else:

            # Load the data from the pickle file if it already exists
            self.load_pickle_dict(pickle_file, pickle_dir=pickle_dir)

            # However, overwrite the pickle and webpage directories as we should be able to load these same pickle files on different systems
            self.pickle_dir = pickle_dir
            self.webpage_dir = webpage_dir


    def calculate_metrics(self, nworkers=1, do_logging=True, use_multiprocessing=True, delete_intermediate_pkl_files=True, keep_unnecessary_calculations=False):
        '''
        Calculate the P values (and Z scores) from the coordinates of the species in every ROI in every slide.

        The assumed distributions are Poisson for the densities and binomial for the PMFs (see summary_of_formulas.lyx) for the formula and physical notebook notes on 12-16-20 for the derivations.

        As seen in the Lyx/PDF file, the precise random variable that is Poisson-distributed (S) is the total number of neighbors of a particular species over all n samples (i.e., over all n centers of a particular species). This is converted to "density" by dividing by n, but then note this is not technically Poisson distributed; e.g., the Poisson distribution applies to only discrete random variables. The corresponding "average number of neighbors around the n centers" is what is typically referred to as "density," which is more technically thought of as the average number of neighbors DIVIDED BY an area. Note that <S> = n * lambda, where lambda is basically the number of neighbors (of a certain species) in the ROI, divided by the ROI area, times the area of the slice under consideration. Note that S runs from 0, 1, 2, ..., infinity.

        The precise random variable that is binomial-distributed (Z_j) is the number of centers (of a particular species) having j neighbors (of a particular species). This (divided by n) is the probability mass function (PMF) of the random variable Y, which is the number of neighbors (of a particular species) around a center (of a particular species). (Note that Y is Poisson-distributed with mean lambda and S = Sum[Y_i], where i goes from 1, 2, ..., n, and Y_i is the number of neighbors around center i.) Note that <Z_j> = n * lambda^j*exp(-lambda)/j!. Note that Z_j runs from 0, 1, 2, ..., n. n is, as above, the total number of centers of a particular species.

        Assuming these distributions, the null hypotheses are that means are as defined above. The "left" or "less" alternative hypotheses are that the means are "less" than (i.e., to the "left" of) these stated values. The "right" or "greater" alternative hypotheses are that the means are "greater" than (i.e., to the "right" of) these stated values. The P values are then the probabilities of observing values of S or Z_j as extreme as (in the directions of the alternative hypotheses) the observed instances of S or Z_j under the assumptions of the null hypotheses.

        Notes:

          * Took ~47 (later: 65) minutes on laptop
          * Units here should be in the same units provided in Consolidated_data.txt (which were originally half-microns, i.e., for dr=8, each slice is 4 microns thick)
          * nworkers should probably be the number of CPUs allocated by SLURM less 1 just to be safe
        '''

        # Import relevant libraries
        import os
        import subprocess

        # Set variables already defined as attributes
        unique_slides = self.unique_slides
        pickle_dir = self.pickle_dir
        plotting_map = self.plotting_map
        nslices = self.nslices
        thickness = self.thickness
        n_neighs = self.n_neighs
        radius_instead_of_knn = self.radius_instead_of_knn
        min_coord_spacing = self.min_coord_spacing
        df_data_by_roi = self.df_data_by_roi
        use_analytical_significance = self.use_analytical_significance

        # Set some attributes from the method parameters
        self.k_max = nslices
        self.dr = thickness

        # Constants
        pickle_file = 'calculated_metrics.pkl'

        # If the pickle file doesn't already exist...
        if not os.path.exists(os.path.join(pickle_dir, pickle_file)):

            # Print what we're doing
            print('Calculating metrics...')

            # Experiment-wide variables
            all_species_list = [x[0] for x in plotting_map]
            nall_species = len(all_species_list)
            nrois = len(df_data_by_roi)
            log_file = 'calculated_metrics.log'

            # ---- Calculate the metrics for all the ROIs that haven't already been calculated, saving the results in individual pickle files

            # Determine the ROIs whose metrics need to be calculated, i.e., those whose corresponding pickle files are not present on the filesystem
            retval = subprocess.run(['ls {}/calculated_metrics-roi_index_*.pkl'.format(pickle_dir)], shell=True, capture_output=True)
            roi_ids_not_present = set(range(nrois)) - set([int(x.split('.pkl')[0].split('_')[-1]) for x in retval.stdout.decode().split('\n')[:-1]])

            # Generate a list of tuple arguments each of which is inputted into calculate_metrics_for_roi() to be run by a single worker
            constant_tuple = (pickle_dir, nslices, thickness, n_neighs, radius_instead_of_knn, min_coord_spacing, all_species_list, nall_species, do_logging, use_analytical_significance, df_data_by_roi, keep_unnecessary_calculations, nworkers)
            # list_of_tuple_arguments = [constant_tuple + (x,) for x in range(nrois)]  # doing it this lazy way potentially messes up the multiprocessing module, causing too many unnecessary-to-be-calculated ROIs to be sent into the Pool, causing only a single worker to actually be used
            list_of_tuple_arguments = [constant_tuple + (x,) for x in roi_ids_not_present]

            # Since Squidpy commandeers multiprocessing, completely disable it if Squidpy has been requested
            if use_analytical_significance:
                # Farm out the metrics calculations to the worker CPUs. This ensures that a pickle file gets created for each ROI
                utils.execute_data_parallelism_potentially(function=calculate_metrics_for_roi, list_of_tuple_arguments=list_of_tuple_arguments, nworkers=(0 if not use_multiprocessing else nworkers), task_description='calculation of ROI metrics (Poisson)')
            else:  # Squidpy is being requested
                print('Running {} function calls using 1 worker WITHOUT the multiprocessing module because Squidpy is being employed, which commandeers threads'.format(len(list_of_tuple_arguments)))
                utils.execute_data_parallelism_potentially(function=calculate_metrics_for_roi, list_of_tuple_arguments=list_of_tuple_arguments, nworkers=0, task_description='calculation of ROI metrics (permutation test)')

            # ---- Load the resulting individual pickle files into a new, single pickle file called calculated_metrics.pkl

            # For each slide...
            data_by_slide = []
            roi_pickle_files = []  # store filenames to delete later
            roi_log_files = []
            for uslide in unique_slides:
                print('Reading slide ' + uslide + '...')

                # Get the unique ROIs in the current slide
                unique_rois = df_data_by_roi[df_data_by_roi['unique_slide'] == uslide]['unique_roi'].unique()  # note the .unique() is likely unneeded

                # For each ROI in the slide...
                data_by_roi = []
                for uroi in unique_rois:
                    print('  Reading ROI ' + uroi + '...')

                    # Load the appropriate pickle file
                    roi_index = df_data_by_roi.loc[df_data_by_roi['unique_roi'] == uroi, :].index[0]
                    roi_pickle_file = 'calculated_metrics-roi_index_{:06}.pkl'.format(roi_index)
                    roi_pickle_files.append(roi_pickle_file)
                    roi_log_files.append('calculated_metrics-roi_index_{:06}.log'.format(roi_index))
                    roi_data_item = load_pickle(pickle_dir, roi_pickle_file)

                    # Save the loaded data
                    data_by_roi.append(roi_data_item)
                data_by_slide.append([uslide, unique_rois, data_by_roi])  # save the current slide data and the inputted parameters

            # Create the single pickle file saving all the data
            make_pickle(data_by_slide, pickle_dir, pickle_file)

            # Concatenate all metrics calculation log files (one per ROI) into a single log file
            logs_dir = os.path.join('.', 'output', 'logs')
            if not os.path.exists(logs_dir):
                os.makedirs(logs_dir)
            with open(os.path.join(logs_dir, log_file), 'w') as outfile:
                for roi_log_file in roi_log_files:
                    with open(os.path.join(pickle_dir, roi_log_file)) as infile:
                        outfile.write(infile.read())

            # If the overall pickle file was successfully created, delete all intermediate pickle files for the ROIs
            if delete_intermediate_pkl_files:
                if os.path.exists(os.path.join(pickle_dir, pickle_file)):
                    for roi_pickle_file in roi_pickle_files:
                        os.remove(os.path.join(pickle_dir, roi_pickle_file))

            # If the overall log file was successfully created, delete all intermediate log files for the ROIs
            if os.path.exists(os.path.join(logs_dir, log_file)):
                for roi_log_file in roi_log_files:
                    os.remove(os.path.join(pickle_dir, roi_log_file))

        # If the pickle file already exists, load it
        else:
            data_by_slide = load_pickle(pickle_dir, pickle_file)

        # Save the calculated data as a property of the class object
        self.metrics = data_by_slide

        # Create a density P value Pandas dataframe from the just-calculated metrics
        self.flatten_density_pvals()


    def save_figs_and_corresp_data(self, plot_real_data=True, pval_figsize=(8, 12), log_pval_range=(-40, 0), calculate_empty_bin_pvals=True, max_nbins=None, roi_figsize=(6, 4), marker_size_step=0.80, pval_dpi=150, alpha=1, roi_dpi=200, square=True, yticklabels=2, pickle_dir=None, save_individual_pval_plots=True, edgecolors='k', default_marker_size_fac=1, yaxis_dir=1, nworkers=1):
        '''

        NOTE: This is likely no longer used!

        Create and save all the figures and the data corresponding to the figures (i.e., the left and right P values for both the densities and PMFs), i.e., the actual data that is plotted. This is a version of save_figs_and_corresp_data() that uses an arbitrary number of CPUs.

        This way, I can see exactly what the read-in data actually look like and therefore trust them more.

        * alpha is the transparency for the circles in the ROI plots (0 is fully transparent, 1 is fully opaque)
        * marker_size_step=0.80 means the radius should be 80% larger for cells plotted behind other cells
        '''

        # Antonio's fix to enable plot generation in SLURM's batch mode
        import matplotlib
        matplotlib.use('Agg')

        # Import relevant libraries
        import os
        import pandas as pd
        import subprocess
        import multiprocessing as mp

        # Set a default value for this; note the input parameter to this call can be something like max_nbins=np.max([slices1.max_nbins_over_exp, slices2.max_nbins_over_exp])
        if max_nbins is None:
            max_nbins = self.max_nbins_over_exp
        self.max_nbins_used = max_nbins

        # Set variables already defined as attributes
        plotting_map = self.plotting_map
        metrics = self.metrics
        num_colors = self.num_colors
        if pickle_dir is None:
            webpage_dir = self.webpage_dir
            pickle_dir = self.pickle_dir
        else:
            webpage_dir = pickle_dir
        mapping_dict = self.mapping_dict
        metrics = self.metrics
        contents = self.contents
        coord_units_in_microns = self.coord_units_in_microns

        # Define the directory holding all the images for the webpage and the filename of the file holding all the corresponding figure data
        webpage_dir = os.path.join(webpage_dir, ('real' if plot_real_data else 'simulated'))
        pickle_file = 'figure_data-{}.pkl'.format(('real' if plot_real_data else 'simulated'))

        # If the pickle file doesn't already exist...
        pickle_pathname = os.path.join(pickle_dir, pickle_file)
        if not os.path.exists(pickle_pathname):

            # Experiment-wide variables
            all_species_list = [x[0] for x in plotting_map]
            nall_species = len(all_species_list)
            if mapping_dict is not None:  # species_names = [x[1] for x in plotting_map]
                species_names = [get_descriptive_cell_label(x[1], mapping_dict)[0] for x in plotting_map]
            else:
                species_names = [phenotypes_to_string(x[1]) for x in plotting_map]

            # Extract the correct number of colors from the default color palette
            ielem = 0
            colors = []
            for elem in matplotlib.rcParams['axes.prop_cycle']():
                color = elem['color']
                colors.append(color)
                ielem = ielem + 1
                if ielem == num_colors:
                    break
            default_marker_size = matplotlib.rcParams['lines.markersize'] * default_marker_size_fac

            # Since at this point we're probably saving some images, ensure their parent directory exists
            os.makedirs(webpage_dir, exist_ok=True)

            # Define the experiment-wide metadata to save and initialize the filedata
            metadata = {
                'webpage_dir': webpage_dir, 'pickle_pathname': pickle_pathname, 'num_slides': len(metrics), 'plot_real_data': plot_real_data, 'all_species_ids': all_species_list, 'all_species_names': species_names, 'nall_species': nall_species,
                'roi_dpi': roi_dpi, 'roi_figsize': roi_figsize,
                'pval_dpi': pval_dpi, 'pval_figsize': pval_figsize, 'log_pval_range': log_pval_range, 'calculate_empty_bin_pvals': calculate_empty_bin_pvals, 'max_nbins': max_nbins
            }

            # This block simply flattens out all the metrics data to the level of the ROIs, instead of the slide-then-ROI heirarchy
            # Note that the ROI order is the same as in contents, which is out of order due to the decompounding process
            roi_index_holder = []
            roi_name_holder = []
            roi_data_holder = []
            slide_name_holder = []
            num_rois_in_slide_holder = []
            for slide_data in metrics:  # for each slide of metrics data...
                num_rois_in_slide = len(slide_data[1])  # get the number of ROIs in the current slide
                for roi_name, roi_data in zip(slide_data[1], slide_data[2]):  # for each ROI of metrics data...
                    roi_index = contents[contents['unique_roi'] == roi_name].index[0]  # use the "contents" dataframe to determine the contents index of the current ROI
                    slide_name = contents.loc[roi_index, 'unique_slide']  # get the corresponding slide name (not currently used)
                    roi_index_holder.append(roi_index)
                    roi_name_holder.append(roi_name)
                    roi_data_holder.append(roi_data)
                    slide_name_holder.append(slide_name)
                    num_rois_in_slide_holder.append(num_rois_in_slide)
            df = pd.DataFrame({'roi_name': roi_name_holder, 'roi_data': roi_data_holder, 'slide_name': slide_name_holder, 'num_rois_in_slide': num_rois_in_slide_holder}, index=roi_index_holder).sort_index()  # create a Pandas dataframe holding all of the flattened data in the same order as in the contents dataframe

            # Determine the ROIs whose figures need to be calculated and corresponding data saved, based only on the data, not the figures, i.e., those whose corresponding pickle files are not present on the filesystem
            retval = subprocess.run(['ls {}/figure_data-{}-roi_index_*.pkl'.format(pickle_dir, ('real' if plot_real_data else 'simulated'))], shell=True, capture_output=True)
            roi_ids_not_present = set(range(len(contents))) - set([int(x.split('.pkl')[0].split('_')[-1]) for x in retval.stdout.decode().split('\n')[:-1]])

            # Generate a list of tuple arguments each of which is inputted into save_figs_and_corresp_data_for_roi() to be run by a single worker
            constant_tuple = (pickle_dir, plotting_map, roi_figsize, pval_figsize, colors, marker_size_step, default_marker_size, roi_dpi, mapping_dict, coord_units_in_microns, alpha, edgecolors, yaxis_dir, webpage_dir, plot_real_data, nall_species, species_names, log_pval_range, calculate_empty_bin_pvals, max_nbins, square, yticklabels, all_species_list, save_individual_pval_plots, pval_dpi, df)
            list_of_tuple_arguments = [constant_tuple + (x,) for x in roi_ids_not_present]

            # Farm out the figure creation and corresponding data saving to the worker CPUs. This ensures that a pickle file gets created for each ROI
            print('Running {} function calls using {} workers'.format(len(list_of_tuple_arguments), nworkers))
            with mp.get_context('spawn').Pool(nworkers) as pool:
                pool.map(save_figs_and_corresp_data_for_roi, list_of_tuple_arguments)

            # Lists to contain data to save and intermediate files to delete
            filedata = []  # this is the main data we want to save
            roi_pickle_files = []  # these intermediate files should later be deleted

            # Read in all valid center-neighbor pairs for each ROI in each slide; this is in the same hierarchical slide-ROI order, as opposed to the order in contents. Order shouldn't matter, but still
            for islide, slide_data in enumerate(metrics):  # for each slide of metrics data...
                print('Reading slide {} of {}...'.format(islide + 1, len(metrics)))
                for roi_name in slide_data[1]:  # for each ROI of metrics data...

                    # Determine and print the current ROI index and pickle file
                    roi_index = contents[contents['unique_roi'] == roi_name].index[0]  # use the "contents" dataframe to determine the contents index of the current ROI
                    roi_pickle_file = 'figure_data-{}-roi_index_{:06}.pkl'.format(('real' if plot_real_data else 'simulated'), roi_index)
                    print('  Reading data for ROI {} from file {}...'.format(roi_name, roi_pickle_file))

                    # Save the ROI pickle filename to later delete
                    roi_pickle_files.append(roi_pickle_file)

                    # Load the ROI data
                    filedata_for_roi = load_pickle(pickle_dir, roi_pickle_file)

                    # For each center-neighbor species pair in the current ROI, append the data to an overall list filedata
                    for pair_data in filedata_for_roi:
                        filedata.append(pair_data)

            # Load the individual pickle files into a new, single pickle file called figure_data-{real,simulated}.pkl
            make_pickle((metadata, filedata), pickle_dir, pickle_file)

            # If the overall pickle file was successfully created, delete all intermediate pickle files for the ROIs
            if os.path.exists(os.path.join(pickle_dir, pickle_file)):
                for roi_pickle_file in roi_pickle_files:
                    os.remove(os.path.join(pickle_dir, roi_pickle_file))

        # If the overall pickle file already exists...
        else:

            # Read in the figure data from disk
            (metadata, filedata) = load_pickle(pickle_dir, pickle_file)

        # Save the calculated data as a property of the class object
        self.metadata = metadata
        self.filedata = filedata


    def plot_rois(self, plot_real_data=True, roi_figsize=(6, 4), marker_size_step=0.80, alpha=0.5, roi_dpi=150, edgecolors='k', default_marker_size_fac=1, yaxis_dir=-1, nworkers=1, use_multiprocessing=True):
        '''
        asdf

        * alpha is the transparency for the circles in the ROI plots (0 is fully transparent, 1 is fully opaque)
        * marker_size_step=0.80 means the radius should be 80% larger for cells plotted behind other cells
        '''

        # # Antonio's fix to enable plot generation in SLURM's batch mode
        # import matplotlib
        # matplotlib.use('Agg')

        # Import relevant libraries
        import os
        import subprocess
        import matplotlib

        # Set variables already defined as attributes
        plotting_map = self.plotting_map
        num_colors = self.num_colors
        webpage_dir = self.webpage_dir
        mapping_dict = self.mapping_dict
        df_data_by_roi = self.df_data_by_roi
        coord_units_in_microns = self.coord_units_in_microns

        # Define the directory holding all the images for the webpage and the filename of the file holding all the corresponding figure data
        savedir = os.path.join(webpage_dir, 'roi_plots')
        if not os.path.exists(savedir):
            os.makedirs(savedir)

        # Extract the correct number of colors from the default color palette
        ielem = 0
        colors = []
        for elem in matplotlib.rcParams['axes.prop_cycle']():
            color = elem['color']
            colors.append(color)
            ielem = ielem + 1
            if ielem == num_colors:
                break
        default_marker_size = matplotlib.rcParams['lines.markersize'] * default_marker_size_fac

        # Determine the ROIs whose figures need to be plotted based on those whose corresponding image files are not present on the filesystem
        retval = subprocess.run(['ls {}/roi_plot_*.png'.format(savedir)], shell=True, capture_output=True)
        roi_ids_not_present = set(range(len(df_data_by_roi))) - set([int(x.split('.png')[0].split('_')[-1]) for x in retval.stdout.decode().split('\n')[:-1]])

        # Generate a list of tuple arguments each of which is inputted into plot_single_roi() to be run by a single worker
        constant_tuple = (plotting_map, roi_figsize, colors, marker_size_step, default_marker_size, roi_dpi, mapping_dict, coord_units_in_microns, alpha, edgecolors, yaxis_dir, savedir, df_data_by_roi)
        list_of_tuple_arguments = [constant_tuple + (x,) for x in roi_ids_not_present]

        # Farm out the figure creation to the worker CPUs. This ensures that an image file gets (or already is) created for each ROI
        if 'roi_plot_fig_pathname' not in df_data_by_roi:
            df_data_by_roi['roi_plot_fig_pathname'] = ''
        utils.execute_data_parallelism_potentially(function=plot_single_roi, list_of_tuple_arguments=list_of_tuple_arguments, nworkers=(0 if not use_multiprocessing else nworkers), task_description='plotting of the ROIs')

    def calculate_max_nbins(self):
        '''
        Calculate the maximum number of possible bins (of numbers neighbors) over the entire experiment, for both the real and simulated data.
        '''

        # Set variables already defined as attributes
        plotting_map = self.plotting_map
        metrics = self.metrics

        # Experiment-wide variables
        nall_species = len(plotting_map)

        # Initialize the variable of interest
        max_nbins_over_exp = 0

        # For each slide...
        for slide_data in metrics:

            # For each ROI in the slide...
            for roi_data in slide_data[2]:  # get the metrics data for the current ROI of the current slide

                # Determine which dataset to plot, each of which is (nall_species, nall_species, nslices)
                # The elements of the dataset ndarray are either a tuple (if both the center and neighbor species are in the ROI) or None
                for data_to_plot in roi_data[1:3]:

                    # Determine whether both the center and neighbor species are in the dataset
                    center_and_neighbor_species_in_dataset = ((data_to_plot != None).sum(axis=2)) > 0  # (nall_species, nall_species)... since we want to do element-wise comparisons here, don't listen to linting when it says the right way to do the comparison is "data_to_plot is not None"

                    # For all combinations of centers and neighbors...
                    for icenter_spec in range(nall_species):
                        for ineighbor_spec in range(nall_species):

                            # If data exist for at least one of the slices for the current center/neighbor combination...
                            if center_and_neighbor_species_in_dataset[icenter_spec, ineighbor_spec]:

                                # Create the P value figure containing four heatmaps and return the figure object
                                # ^---This must be an incorrect comment right?!
                                max_nbins_over_exp = get_max_nbins_for_center_neighbor_pair(data_to_plot[icenter_spec, ineighbor_spec, :], max_nbins_over_exp)

        # Set the maximum number of bins as a property of the object itself and print this result
        self.max_nbins_over_exp = max_nbins_over_exp
        print('Calculated maximum number of bins over the entire experiment: {}'.format(max_nbins_over_exp))

        return(max_nbins_over_exp)


    def load_pickle_class(self, pickle_file, pickle_dir=None):
        '''
        Load some data from a pickle file ("class" just refers to this function being part of the TIMECellInteraction class)
        '''
        if pickle_dir is None:
            pickle_dir = self.pickle_dir
        return(load_pickle(pickle_dir, pickle_file))


    def load_pickle_dict(self, pickle_file, pickle_dir=None):
        '''
        Load a bunch of values to the self object from a pickle file by way of a dictionary
        '''
        dict2load = self.load_pickle_class(pickle_file, pickle_dir=pickle_dir)
        for key in dict2load:
            val = dict2load[key]
            setattr(self, key, val)


    def make_pickle_dict(self, vars2save, local_dict, pickle_file):
        '''
        Make a pickle file of a dictionary of data
        '''
        dict2save = {}
        for key in vars2save:
            if key in local_dict:
                val = local_dict[key]
                setattr(self, key, val)
            else:
                val = getattr(self, key)
            dict2save.update({key: val})
        make_pickle(dict2save, self.pickle_dir, pickle_file)


    def preprocess_dataframe(self, allow_compound_species):
        '''
        Preprocess the initial Pandas dataframe from Consolidata_data.txt (or a simulated one for simulated data) by creating another column (Species int) specifying a unique integer identifying the cell type
        If requested, remove compound species, and return the list of single-protein "phenotypes" contained in the data
        '''

        # Import relevant module
        import numpy as np

        # Due to my pre-Pandas-learned days and thus implementing functionality below suboptimally, we should ensure the dataset has a sorted index
        self.data = self.data.reset_index(drop=True)

        # Preprocess the pandas dataframe in various ways
        data_phenotypes = self.data.filter(regex='^Phenotype ')  # get just the "Phenotype " columns
        # data_phenotypes = data_phenotypes.reset_index(drop=True)  # this is bad if the input data is not already sorted and should be commented out! However, note that was likely already here because of the non-optimal way we drop trivial objects below using the drop() method. Thus, we just need to ensure that the index is reset at the beginning of this function, as we do above
        phenotype_cols = list(data_phenotypes.columns)  # get a list of those column names
        phenotypes = np.array([x.replace('Phenotype ', '') for x in phenotype_cols])  # extract just the phenotypes from that list
        n_phenotypes = len(phenotypes)  # get the number of possible phenotypes in the datafile
        self.data['Species string'] = data_phenotypes.map(lambda x: '1' if (str(x)[-1] == '+') else '0').apply(lambda x: ''.join(list(x)), axis=1)  # add a column to the original data that tells us the unique "binary" species string of the species corresponding to that row/cell
        check_roi_order(self.data, 'prior to dropping invalid objects')
        self.data = self.data.drop(np.nonzero((self.data['Species string'] == ('0' * n_phenotypes)).to_numpy())[0])  # delete rows that all have '...-' as the phenotype or are blank. As mentioned above, this is a non-optimal way of performing this; a better way is likely something like:

        # df_phenotypes = df.filter(regex='^Phenotype ')
        # num_phenotypes = df_phenotypes.shape[1]
        # df_pared = df[df_phenotypes.apply(lambda x: ''.join(x), axis='columns') != '-' * num_phenotypes]
        check_roi_order(self.data, 'after dropping invalid objects')

        self.data = self.data.reset_index(drop=True)  # reset the indices
        self.data['Species int'] = self.data['Species string'].apply(lambda x: int(x, base=2))  # add an INTEGER species column

        # Remove compound species if requested
        if not allow_compound_species:
            self.remove_compound_species()
            self.remove_compound_species()  # ensure nothing happens

        return(phenotypes)


    def remove_compound_species(self):
        '''
        For each compound species ('Species int' not just a plain power of two), add each individual phenotype to the end of the dataframe individually and then delete the original compound entry
        '''

        # Import relevant modules
        import numpy as np
        import pandas as pd

        # Get the species IDs
        x = np.array(self.data['Species int'])
        print('Data size:', len(self.data))

        # Determine which are not powers of 2, i.e., are compound species
        powers = np.log2(x)
        compound_loc = np.nonzero(powers != np.round(powers))[0]
        ncompound = len(compound_loc)
        print('  Compound species found:', ncompound)

        check_roi_order(self.data, 'prior to decompounding')

        # If compound species exist...
        if ncompound > 0:

            print('  Removing compound species from the dataframe...')

            # Get a list of tuples each corresponding to compound species, the first element of which is the row of the compound species, and the second of which is the species IDs of the pure phenotypes that make up the compound species
            compound_entries = [(cl, 2**np.nonzero([int(y) for y in bin(x[cl])[2:]][-1::-1])[0]) for cl in compound_loc]

            # For each compound species...
            data_to_add = []
            for index, subspecies in compound_entries:

                # For each pure phenotype making up the compound species...
                for subspec in subspecies:

                    # You have to put this here instead of outside this loop for some weird reason! Even though you can see the correct change made to series and EVEN TO DATA_TO_ADD by looking at series and data_to_add[-1] below, for some Godforsaken reason the actual data_to_add list does not get updated with the change to 'Species int' when you print data_to_add at the end of both these loops, and therefore the data that gets added to the data dataframe contains all the same 'Species string' values, namely the last one assigned. Thus, we are actually adding the SAME species to multiple (usually 2) spatial points, so that the even-spacing problem arises.
                    series = self.data.iloc[index].copy()

                    # Note the only field I'm updating here is "Species int" and NOT the phenotype columns nor the "Species string" column. "Species int" is the only one that matters. So don't get confused when I print the dataframe and see apparently incorrect other columns... I am copying their values!
                    series['Species int'] = subspec  # set the species ID of the series data to that of the current phenotype
                    data_to_add.append(series)  # add the data to a running list

            # Add all the data in the list to the dataframe
            # self.data = self.data.append(data_to_add, ignore_index=True)
            self.data = pd.concat([self.data, pd.DataFrame(data_to_add)], ignore_index=True)  # needed because Pandas 2.0 deprecates .append()
            print('  Added rows:', len(data_to_add))
            print('  Data size:', len(self.data))

            # Delete the original compound species entries
            self.data = self.data.drop(compound_loc)
            self.data = self.data.reset_index(drop=True)
            print('  Deleted rows:', len(compound_loc))
            print('  Data size:', len(self.data))

        check_roi_order(self.data, 'post decompounding')


    def average_over_rois(self, plot_real_data=True, log_pval_range=(-200, 0), figsize=(14, 4), dpi=150, img_file_suffix='', regular_pval_figsize=(8, 12), square=True, yticklabels=2, pval_dpi=150, plot_summary_dens_pvals=True, plot_main_pvals=True, write_csv_files=True, write_pickle_datafile=True, start_slide=None, num_slides=None):
        '''
        Perform a weighted geometric mean of the density and PMF P values over the valid ROIs, saving the corresponding figures as PNG files and data as CSV files.

        Old:

        Read all the data corresponding to the P value plots into a Pandas dataframe, using this structure to write an ndarray holding the left and right density (not PMF) P values (for a single-slice dataset) for every patient (slide) and center/neighbor combination in the entire experiment, performing a weighted average of the P values over all the ROIs for each patient.

        The resulting array can then be read back in to other functions (such as create_summary_plots(), below) to plot these results in heatmaps or to write them to CSV files.

        From def create_summary_plots():
          Read in the summary ndarray from create_summary_array(), above, to plot and save heatmaps of all these data.
          Call like, e.g., create_summary_plots(dataset_name, project_dir, plot_real_data, summary_arr, unique_slides, metadata['all_species_names']).

        From def create_summary_csv_files():
          Write the left and right summary P values to a CSV file.
          Call like:
          tci.create_summary_csv_files(dataset_name, project_dir, plot_real_data, summary_arr, unique_slides, metadata['all_species_names'])
        '''

        # Import relevant libraries
        import pandas as pd
        import numpy as np
        import seaborn as sns
        import matplotlib.pyplot as plt
        import os

        # Define variables already defined as attributes
        filedata = self.filedata
        nslices = self.nslices
        max_nbins = self.max_nbins_used
        pickle_dir = self.pickle_dir
        webpage_dir = self.webpage_dir
        all_species_names = self.all_species_names
        nall_species = self.nall_species
        all_species_ids = self.all_species_ids

        # Create the Pandas dataframe from the list of dictionaries that were outputted by save_figs_and_corresp_data()
        df = pd.DataFrame()
        df = df.append(filedata)  # may need to fix due to deprecation of .append(); use pd.concat() instead

        # Add a slide column to the dataframe and get a list of the unique slides while maintaining order
        df['slide_name'] = df['roi_name'].apply(lambda x: x.split('_')[0])
        unique_slides = df['slide_name'].unique()
        nunique_slides = len(unique_slides)

        # Determine which unique slides to plot in the main plotting functionality below
        if start_slide is None:
            start_slide = 0
            num_slides = nunique_slides
        stop_slide = start_slide + num_slides
        unique_slides_to_plot = unique_slides[start_slide:stop_slide]

        # Define the directory holding all the images of the averaged data for the webpage and the filename of the file holding all the corresponding data
        webpage_dir = os.path.join(webpage_dir, ('real' if plot_real_data else 'simulated'))
        pickle_file = 'averaged_data-{}.pkl'.format(('real' if plot_real_data else 'simulated'))

        # If the pickle file doesn't already exist...
        if not os.path.exists(os.path.join(pickle_dir, pickle_file)):

            # Initialize the arrays of interest
            nmatches_holder    = np.zeros((nunique_slides, nall_species, nall_species, nslices))
            log_dens_pvals_avg = np.zeros((nunique_slides, nall_species, nall_species, 2, nslices))
            log_pmf_pvals_avg  = np.zeros((nunique_slides, nall_species, nall_species, max_nbins, 2, nslices))

            # For every unique slide...
            for islide, slide in enumerate(unique_slides[:]):

                # For every possible center species...
                for icenter_spec, center_spec in enumerate(all_species_ids[:]):

                    # For every possible neighbor species...
                    for ineighbor_spec, neighbor_spec in enumerate(all_species_ids[:]):

                        # Get a temporary dataframe containing the current slide/center/neighbor combination (this will usually be five rows, one per ROI)
                        df_tmp = df[(((df['slide_name'] == slide) & (df['center_species_id'] == center_spec) & (df['neighbor_species_id'] == neighbor_spec)))]
                        nrows = len(df_tmp)

                        # For every slice...
                        for islice in range(nslices):

                            # Initialize and populate the three temporary arrays
                            nvalid_centers_holder = np.zeros((nrows,))
                            log_dens_pvals = np.zeros((nrows, 2))
                            log_pmf_pvals = np.zeros((nrows, max_nbins, 2))
                            for irow, (nvalid_centers_per_slice, left_log_dens_pvals, right_log_dens_pvals, left_log_pmf_pvals, right_log_pmf_pvals) in enumerate(zip(df_tmp['nvalid_centers_per_slice'], df_tmp['left_log_dens_pvals'], df_tmp['right_log_dens_pvals'], df_tmp['left_log_pmf_pvals'], df_tmp['right_log_pmf_pvals'])):
                                nvalid_centers_holder[irow] = nvalid_centers_per_slice[islice]
                                log_dens_pvals[irow, :] = np.array([left_log_dens_pvals[0, islice], right_log_dens_pvals[0, islice]])  # (2,)
                                log_pmf_pvals[irow, :, :] = np.c_[left_log_pmf_pvals[:, islice], right_log_pmf_pvals[:, islice]]  # (max_nbins,2)

                            # Determine the rows in the temporary dataframe that have at least 1 valid center
                            matches = nvalid_centers_holder >= 1  # (nrows,)
                            nmatches = matches.sum()
                            nmatches_holder[islide, icenter_spec, ineighbor_spec, islice] = nmatches

                            # Perform the weighted averaging over the ROIs
                            if nmatches >= 1:
                                nvalid_centers_tmp = nvalid_centers_holder[matches][:, np.newaxis]  # (nmatches, 1)
                                log_dens_pvals_avg[islide, icenter_spec, ineighbor_spec, :, islice] = (nvalid_centers_tmp * log_dens_pvals[matches, :]).sum(axis=0) / nvalid_centers_tmp.sum(axis=0)  # (2,)
                                nvalid_centers_tmp = nvalid_centers_tmp[:, :, np.newaxis]  # (nmatches, 1, 1)
                                log_pmf_pvals_avg[islide, icenter_spec, ineighbor_spec, :, :, islice] = (nvalid_centers_tmp * log_pmf_pvals[matches, :, :]).sum(axis=0) / nvalid_centers_tmp.sum(axis=0)  # (max_nbins, 2)
                            else:
                                log_dens_pvals_avg[islide, icenter_spec, ineighbor_spec, :, islice] = None
                                log_pmf_pvals_avg[islide, icenter_spec, ineighbor_spec, :, :, islice] = None

            # Set the log of the P value range for the color plotting
            vmin = log_pval_range[0]
            vmax = log_pval_range[1]

            # Set the (negative) infinite values to the darkest color (or else they won't be plotted, as inf values are not plotted)
            log_dens_pvals_avg[np.isneginf(log_dens_pvals_avg)] = vmin
            log_pmf_pvals_avg[np.isneginf(log_pmf_pvals_avg)] = vmin

            # Initialize the figure and axes
            fig, ax = plt.subplots(nrows=1, ncols=2, figsize=figsize)

            if plot_summary_dens_pvals:

                # Plot the average density P values for every slice
                for islice in range(nslices):

                    # For every unique slide...
                    for islide in range(nunique_slides):

                        # Determine the filename of the figure
                        filename = os.path.join(webpage_dir, 'average_density_pvals-{}-{}-slice_{:02d}_of_{:02d}{}.png'.format(('real' if plot_real_data else 'simulated'), unique_slides[islide], islice + 1, nslices, img_file_suffix))

                        # Reset the figure/axes
                        fig.clf()
                        ax = fig.subplots(nrows=1, ncols=2)

                        # Plot the log of the left/less P values
                        sns.heatmap(log_dens_pvals_avg[islide, :, :, 0, islice], vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[0], cbar=True, xticklabels=all_species_names, yticklabels=all_species_names, square=True)
                        ax[0].set_title('log10(\"less\" density pvals)')
                        ax[0].set_xlabel('Neighbor species')
                        ax[0].set_ylabel('Center species')

                        # Plot the log of the right/greater P values
                        sns.heatmap(log_dens_pvals_avg[islide, :, :, 1, islice], vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[1], cbar=True, xticklabels=all_species_names, yticklabels=all_species_names, square=True)
                        ax[1].set_title('log10(\"greater\" density pvals)')
                        ax[1].set_xlabel('Neighbor species')
                        ax[1].set_ylabel('Center species')

                        # Set the figure title to the slide title and ensure the facecolor is white
                        fig.suptitle('Average density P values - {} - {} data - slice {} of {}'.format(unique_slides[islide], ('real' if plot_real_data else 'simulated'), islice + 1, nslices))
                        fig.patch.set_facecolor('white')

                        # Save the figure to disk
                        fig.savefig(filename, dpi=dpi, bbox_inches='tight')

            # Initialize the P value figure
            fig = plt.subplots(nrows=2, ncols=2, figsize=regular_pval_figsize)[0]

            if plot_main_pvals:

                # For every slide, center, and neighbor species...
                # for islide, slide_name in enumerate(unique_slides):
                for islide2, slide_name in enumerate(unique_slides_to_plot):
                    islide = islide2 + start_slide
                    for icenter, center_name in enumerate(all_species_names):
                        for ineighbor, neighbor_name in enumerate(all_species_names):

                            # Initialize the current figure by clearing it and all its axes
                            fig.clf()
                            ax = fig.subplots(nrows=2, ncols=2)

                            # Create the four-axis figure, where the top row has the left and right density P values and the bottom row has the left and right PMF P values

                            # Plot the log10 of the left density P values
                            sns.heatmap(log_dens_pvals_avg[islide, icenter, ineighbor, 0, :][np.newaxis, :], vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[0, 0], cbar=True, yticklabels=True, square=True)
                            ax[0, 0].set_title('log10(\"less\" density pvals)')
                            ax[0, 0].set_xlabel('Slice')

                            # Plot the log10 of the right density P values
                            sns.heatmap(log_dens_pvals_avg[islide, icenter, ineighbor, 1, :][np.newaxis, :], vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[0, 1], cbar=True, yticklabels=True, square=True)
                            ax[0, 1].set_title('log10(\"greater\" density pvals)')
                            ax[0, 1].set_xlabel('Slice')

                            # Plot the log10 of the left PMF P values
                            sns.heatmap(log_pmf_pvals_avg[islide, icenter, ineighbor, :, 0, :], vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[1, 0], cbar=True, yticklabels=yticklabels, square=square)
                            ax[1, 0].set_title('log10(\"less\" PMF pvals)')
                            ax[1, 0].set_xlabel('Slice')
                            ax[1, 0].set_ylabel('Number of neighbors')

                            # Plot the log10 of the right PMF P values
                            sns.heatmap(log_pmf_pvals_avg[islide, icenter, ineighbor, :, 1, :], vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[1, 1], cbar=True, yticklabels=yticklabels, square=square)
                            ax[1, 1].set_title('log10(\"greater\" PMF pvals)')
                            ax[1, 1].set_xlabel('Slice')
                            ax[1, 1].set_ylabel('Number of neighbors')

                            # Place a descriptive title on the figure
                            figure_title = 'Averaged-over-ROI, {} data\nSlide: {}\ncenter={}, neighbor={}'.format(('real' if plot_real_data else 'simulated'), slide_name, center_name, neighbor_name)
                            fig.suptitle(figure_title)

                            # Save the figure
                            pvals_fig_filename = 'averaged_pvals_{}_center-{}_neighbor-{}.png'.format(slide_name, all_species_ids[icenter], all_species_ids[ineighbor])
                            pvals_fig_dirname = os.path.join(webpage_dir, slide_name)
                            pvals_fig_pathname = os.path.join(pvals_fig_dirname, pvals_fig_filename)
                            os.makedirs(pvals_fig_dirname, exist_ok=True)
                            fig.savefig(pvals_fig_pathname, dpi=pval_dpi, bbox_inches='tight')

            # Determine the filename of the CSV files we intend to write
            density_csv_pathname = os.path.join(pickle_dir, 'average_density_pvals-{}.csv'.format(('real' if plot_real_data else 'simulated')))
            pmf_csv_pathname = os.path.join(pickle_dir, 'average_pmf_pvals-{}.csv'.format(('real' if plot_real_data else 'simulated')))

            # Part the slide names into the individual patients and drug treatment
            patient = [int(x.split('-')[0][:-1]) for x in unique_slides]
            pre_or_post = [x.split('-')[0][-1:] for x in unique_slides]
            upatient = []
            upre_or_post = []
            for curr_patient, curr_pre_or_post in zip(patient, pre_or_post):
                if curr_patient not in upatient:
                    upatient.append(curr_patient)
                if curr_pre_or_post not in upre_or_post:
                    upre_or_post.append(curr_pre_or_post)

            # Reshape the average holder for the main arrays to the individual patients and drug treatment, filling in "blanks" if necessary
            impossible_val = 444.4
            log_dens_pvals_avg2 = np.ones((len(upatient), len(upre_or_post), nall_species, nall_species, 2, nslices)) * impossible_val
            log_pmf_pvals_avg2 = np.ones((len(upatient), len(upre_or_post), nall_species, nall_species, max_nbins, 2, nslices)) * impossible_val
            for iunique_slide, (curr_patient, curr_pre_or_post) in enumerate(zip(patient, pre_or_post)):
                upatient_idx = upatient.index(curr_patient)
                upre_or_post_idx = upre_or_post.index(curr_pre_or_post)
                log_dens_pvals_avg2[upatient_idx, upre_or_post_idx, :, :, :, :] = log_dens_pvals_avg[iunique_slide, :, :, :, :]
                log_pmf_pvals_avg2[upatient_idx, upre_or_post_idx, :, :, :, :, :] = log_pmf_pvals_avg[iunique_slide, :, :, :, :, :]
            log_dens_pvals_avg = log_dens_pvals_avg2
            log_pmf_pvals_avg = log_pmf_pvals_avg2

            # Create the Pandas indexes
            index_density = pd.MultiIndex.from_product([upre_or_post, all_species_names, all_species_names, ['left', 'right'], np.arange(nslices) + 1], names=['drug_status', 'center_species', 'neighbor_species', 'pval_type', 'slice_number'])
            index_pmf = pd.MultiIndex.from_product([upre_or_post, all_species_names, all_species_names, np.arange(max_nbins), ['left', 'right'], np.arange(nslices) + 1], names=['drug_status', 'center_species', 'neighbor_species', 'poss_num_neighbors', 'pval_type', 'slice_number'])

            # For the density, Create the Pandas dataframes from the indexes
            df_log_dens_pvals_avg = pd.DataFrame(data=log_dens_pvals_avg.reshape((len(upatient), -1)), index=upatient, columns=index_density).rename_axis('subject')
            df_log_pmf_pvals_avg = pd.DataFrame(data=log_pmf_pvals_avg.reshape((len(upatient), -1)), index=upatient, columns=index_pmf).rename_axis('subject')

            # Write the Pandas dataframes to disk
            if write_csv_files:
                df_log_dens_pvals_avg.to_csv(density_csv_pathname)
                df_log_pmf_pvals_avg.to_csv(pmf_csv_pathname)

            # Save the averaged data to disk
            if write_pickle_datafile:
                make_pickle((nmatches_holder, log_dens_pvals_avg, log_pmf_pvals_avg, df_log_dens_pvals_avg, df_log_pmf_pvals_avg), pickle_dir, pickle_file)

            # Close the figure
            plt.close(fig)

        else:

            # Read in the averaged data from disk
            (nmatches_holder, log_dens_pvals_avg, log_pmf_pvals_avg, df_log_dens_pvals_avg, df_log_pmf_pvals_avg) = load_pickle(pickle_dir, pickle_file)

        return(nmatches_holder, log_dens_pvals_avg, log_pmf_pvals_avg, unique_slides, df_log_dens_pvals_avg, df_log_pmf_pvals_avg)


    # def average_over_rois2(self, plot_real_data=True, log_pval_range=(-200, 0), figsize=(14, 4), dpi=150, img_file_suffix='', regular_pval_figsize=(8, 12), square=True, yticklabels=2, pval_dpi=150, plot_summary_dens_pvals=True, plot_main_pvals=True, write_csv_files=True, write_pickle_datafile=True, start_slide=None, num_slides=None, min_num_valid_centers=1):
    def average_over_rois2(self, plot_real_data=True, log_pval_range=(-200, 0), figsize=(14, 4), dpi=150, img_file_suffix='', plot_summary_dens_pvals=True, write_pickle_datafile=True, start_slide=None, min_num_valid_centers=1, weight_rois_by_num_valid_centers=False):
        '''
        Perform a weighted geometric mean of the density and PMF P values over the valid ROIs, saving the corresponding figures as PNG files and data as CSV files.

        Old:

        Read all the data corresponding to the P value plots into a Pandas dataframe, using this structure to write an ndarray holding the left and right density (not PMF) P values (for a single-slice dataset) for every patient (slide) and center/neighbor combination in the entire experiment, performing a weighted average of the P values over all the ROIs for each patient.

        The resulting array can then be read back in to other functions (such as create_summary_plots(), below) to plot these results in heatmaps or to write them to CSV files.

        From def create_summary_plots():
          Read in the summary ndarray from create_summary_array(), above, to plot and save heatmaps of all these data.
          Call like, e.g., create_summary_plots(dataset_name, project_dir, plot_real_data, summary_arr, unique_slides, metadata['all_species_names']).

        From def create_summary_csv_files():
          Write the left and right summary P values to a CSV file.
          Call like:
          tci.create_summary_csv_files(dataset_name, project_dir, plot_real_data, summary_arr, unique_slides, metadata['all_species_names'])
        '''

        # Import relevant libraries
        # import pandas as pd
        import numpy as np
        import seaborn as sns
        import matplotlib.pyplot as plt
        import os

        # Define variables already defined as attributes
        filedata = self.df_density_pvals
        nslices = self.nslices
        # max_nbins = self.max_nbins_used
        pickle_dir = self.pickle_dir
        webpage_dir = self.webpage_dir
        all_species_names = self.all_species_names
        nall_species = self.nall_species
        all_species_ids = self.all_species_ids

        # Create the Pandas dataframe from the list of dictionaries that were outputted by save_figs_and_corresp_data()
        # df = pd.DataFrame()
        # df = df.append(filedata)
        df = filedata

        # Add a slide column to the dataframe and get a list of the unique slides while maintaining order
        df['slide_name'] = df['roi_name'].apply(lambda x: x.split('_')[0])
        unique_slides = df['slide_name'].unique()
        nunique_slides = len(unique_slides)

        # Determine which unique slides to plot in the main plotting functionality below
        if start_slide is None:
            start_slide = 0
            # num_slides = nunique_slides
        # stop_slide = start_slide + num_slides
        # unique_slides_to_plot = unique_slides[start_slide:stop_slide]

        # Define the directory holding all the images of the averaged data for the webpage and the filename of the file holding all the corresponding data
        webpage_dir = os.path.join(webpage_dir, ('real' if plot_real_data else 'simulated'))
        pickle_file = 'averaged_data-{}.pkl'.format(('real' if plot_real_data else 'simulated'))

        # If the pickle file doesn't already exist...
        if not os.path.exists(os.path.join(pickle_dir, pickle_file)):

            # Initialize the arrays of interest
            nmatches_holder    = np.zeros((nunique_slides, nall_species, nall_species, nslices))
            log_dens_pvals_avg = np.zeros((nunique_slides, nall_species, nall_species, 2, nslices))
            # log_pmf_pvals_avg  = np.zeros((nunique_slides, nall_species, nall_species, max_nbins, 2, nslices))

            # Set the log of the P value range for the color plotting
            vmin = log_pval_range[0]
            vmax = log_pval_range[1]

            # For every unique slide...
            for islide, slide in enumerate(unique_slides[:]):

                print('Averaging slide {}'.format(slide))

                # For every possible center species...
                for icenter_spec, center_spec in enumerate(all_species_ids[:]):

                    # For every possible neighbor species...
                    for ineighbor_spec, neighbor_spec in enumerate(all_species_ids[:]):

                        # Get a temporary dataframe containing the current slide/center/neighbor combination (this will usually be five rows, one per ROI)
                        df_tmp = df[(((df['slide_name'] == slide) & (df['center_species_id'] == center_spec) & (df['neighbor_species_id'] == neighbor_spec)))]
                        nrows = len(df_tmp)

                        # For every slice...
                        for islice in range(nslices):

                            # Initialize and populate the three temporary arrays
                            nvalid_centers_holder = np.zeros((nrows,))
                            log_dens_pvals = np.zeros((nrows, 2))
                            # log_pmf_pvals = np.zeros((nrows, max_nbins, 2))
                            # for irow, (nvalid_centers_per_slice, left_log_dens_pvals, right_log_dens_pvals, left_log_pmf_pvals, right_log_pmf_pvals) in enumerate(zip(df_tmp['nvalid_centers_per_slice'], df_tmp['left_log_dens_pvals'], df_tmp['right_log_dens_pvals'], df_tmp['left_log_pmf_pvals'], df_tmp['right_log_pmf_pvals'])):
                            for irow, (nvalid_centers_per_slice, left_log_dens_pvals, right_log_dens_pvals) in enumerate(zip(df_tmp['nvalid_centers_per_slice'], df_tmp['left_log_dens_pvals'], df_tmp['right_log_dens_pvals'])):
                                nvalid_centers_holder[irow] = nvalid_centers_per_slice[islice]
                                log_dens_pvals[irow, :] = np.array([left_log_dens_pvals[0][islice], right_log_dens_pvals[0][islice]])  # (2,)
                                # log_pmf_pvals[irow, :, :] = np.c_[left_log_pmf_pvals[:, islice], right_log_pmf_pvals[:, islice]]  # (max_nbins,2)

                            # Determine the rows in the temporary dataframe that have at least 1 valid center
                            matches = nvalid_centers_holder >= min_num_valid_centers  # (nrows,)
                            nmatches = matches.sum()
                            nmatches_holder[islide, icenter_spec, ineighbor_spec, islice] = nmatches

                            # Perform the weighted averaging over the ROIs
                            if nmatches >= 1:
                                log_dens_pvals_tmp = log_dens_pvals[matches, :]
                                log_dens_pvals_tmp[np.isneginf(log_dens_pvals_tmp)] = vmin
                                if weight_rois_by_num_valid_centers:
                                    nvalid_centers_tmp = nvalid_centers_holder[matches][:, np.newaxis]  # (nmatches, 1)
                                    log_dens_pvals_avg[islide, icenter_spec, ineighbor_spec, :, islice] = (nvalid_centers_tmp * log_dens_pvals_tmp).sum(axis=0) / nvalid_centers_tmp.sum(axis=0)  # (2,)
                                else:
                                    log_dens_pvals_avg[islide, icenter_spec, ineighbor_spec, :, islice] = (log_dens_pvals_tmp).sum(axis=0) / nmatches  # (2,)
                                # nvalid_centers_tmp = nvalid_centers_tmp[:, :, np.newaxis]  # (nmatches, 1, 1)
                                # log_pmf_pvals_avg[islide, icenter_spec, ineighbor_spec, :, :, islice] = (nvalid_centers_tmp * log_pmf_pvals[matches, :, :]).sum(axis=0) / nvalid_centers_tmp.sum(axis=0)  # (max_nbins, 2)
                            else:
                                log_dens_pvals_avg[islide, icenter_spec, ineighbor_spec, :, islice] = None
                                # log_pmf_pvals_avg[islide, icenter_spec, ineighbor_spec, :, :, islice] = None

            # Set the (negative) infinite values to the darkest color (or else they won't be plotted, as inf values are not plotted)
            # log_dens_pvals_avg[np.isneginf(log_dens_pvals_avg)] = vmin
            # log_pmf_pvals_avg[np.isneginf(log_pmf_pvals_avg)] = vmin

            # Initialize the figure and axes
            fig, ax = plt.subplots(nrows=1, ncols=2, figsize=figsize)

            if plot_summary_dens_pvals:

                # Plot the average density P values for every slice
                for islice in range(nslices):

                    # For every unique slide...
                    for islide in range(nunique_slides):

                        print('Plotting slide {}'.format(islide + 1))

                        # Determine the filename of the figure
                        filename = os.path.join(webpage_dir, 'average_density_pvals-{}-{}-slice_{:02d}_of_{:02d}{}.png'.format(('real' if plot_real_data else 'simulated'), unique_slides[islide], islice + 1, nslices, img_file_suffix))

                        # Reset the figure/axes
                        fig.clf()
                        ax = fig.subplots(nrows=1, ncols=2)

                        # Plot the log of the left/less P values
                        sns.heatmap(log_dens_pvals_avg[islide, :, :, 0, islice], vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[0], cbar=True, xticklabels=all_species_names, yticklabels=all_species_names, square=True)
                        ax[0].set_title('log10(\"less\" density pvals)')
                        ax[0].set_xlabel('Neighbor species')
                        ax[0].set_ylabel('Center species')

                        # Plot the log of the right/greater P values
                        sns.heatmap(log_dens_pvals_avg[islide, :, :, 1, islice], vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[1], cbar=True, xticklabels=all_species_names, yticklabels=all_species_names, square=True)
                        ax[1].set_title('log10(\"greater\" density pvals)')
                        ax[1].set_xlabel('Neighbor species')
                        ax[1].set_ylabel('Center species')

                        # Set the figure title to the slide title and ensure the facecolor is white
                        fig.suptitle('Average density P values - {} - {} data - slice {} of {}'.format(unique_slides[islide], ('real' if plot_real_data else 'simulated'), islice + 1, nslices))
                        fig.patch.set_facecolor('white')

                        # Save the figure to disk
                        fig.savefig(filename, dpi=dpi, bbox_inches='tight')

            # # Initialize the P value figure
            # fig = plt.subplots(nrows=2, ncols=2, figsize=regular_pval_figsize)[0]

            # if plot_main_pvals:

            #     # For every slide, center, and neighbor species...
            #     # for islide, slide_name in enumerate(unique_slides):
            #     for islide2, slide_name in enumerate(unique_slides_to_plot):
            #         islide = islide2 + start_slide
            #         for icenter, center_name in enumerate(all_species_names):
            #             for ineighbor, neighbor_name in enumerate(all_species_names):

            #                 # Initialize the current figure by clearing it and all its axes
            #                 fig.clf()
            #                 ax = fig.subplots(nrows=2, ncols=2)

            #                 # Create the four-axis figure, where the top row has the left and right density P values and the bottom row has the left and right PMF P values

            #                 # Plot the log10 of the left density P values
            #                 sns.heatmap(log_dens_pvals_avg[islide, icenter, ineighbor, 0, :][np.newaxis, :], vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[0, 0], cbar=True, yticklabels=True, square=True)
            #                 ax[0, 0].set_title('log10(\"less\" density pvals)')
            #                 ax[0, 0].set_xlabel('Slice')

            #                 # Plot the log10 of the right density P values
            #                 sns.heatmap(log_dens_pvals_avg[islide, icenter, ineighbor, 1, :][np.newaxis, :], vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[0, 1], cbar=True, yticklabels=True, square=True)
            #                 ax[0, 1].set_title('log10(\"greater\" density pvals)')
            #                 ax[0, 1].set_xlabel('Slice')

            #                 # Plot the log10 of the left PMF P values
            #                 sns.heatmap(log_pmf_pvals_avg[islide, icenter, ineighbor, :, 0, :], vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[1, 0], cbar=True, yticklabels=yticklabels, square=square)
            #                 ax[1, 0].set_title('log10(\"less\" PMF pvals)')
            #                 ax[1, 0].set_xlabel('Slice')
            #                 ax[1, 0].set_ylabel('Number of neighbors')

            #                 # Plot the log10 of the right PMF P values
            #                 sns.heatmap(log_pmf_pvals_avg[islide, icenter, ineighbor, :, 1, :], vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[1, 1], cbar=True, yticklabels=yticklabels, square=square)
            #                 ax[1, 1].set_title('log10(\"greater\" PMF pvals)')
            #                 ax[1, 1].set_xlabel('Slice')
            #                 ax[1, 1].set_ylabel('Number of neighbors')

            #                 # Place a descriptive title on the figure
            #                 figure_title = 'Averaged-over-ROI, {} data\nSlide: {}\ncenter={}, neighbor={}'.format(('real' if plot_real_data else 'simulated'), slide_name, center_name, neighbor_name)
            #                 fig.suptitle(figure_title)

            #                 # Save the figure
            #                 pvals_fig_filename = 'averaged_pvals_{}_center-{}_neighbor-{}.png'.format(slide_name, all_species_ids[icenter], all_species_ids[ineighbor])
            #                 pvals_fig_dirname = os.path.join(webpage_dir, slide_name)
            #                 pvals_fig_pathname = os.path.join(pvals_fig_dirname, pvals_fig_filename)
            #                 os.makedirs(pvals_fig_dirname, exist_ok=True)
            #                 fig.savefig(pvals_fig_pathname, dpi=pval_dpi, bbox_inches='tight')

            # # Determine the filename of the CSV files we intend to write
            # density_csv_pathname = os.path.join(pickle_dir, 'average_density_pvals-{}.csv'.format(('real' if plot_real_data else 'simulated')))
            # pmf_csv_pathname = os.path.join(pickle_dir, 'average_pmf_pvals-{}.csv'.format(('real' if plot_real_data else 'simulated')))

            # # Part the slide names into the individual patients and drug treatment
            # patient = [int(x.split('-')[0][:-1]) for x in unique_slides]
            # pre_or_post = [x.split('-')[0][-1:] for x in unique_slides]
            # upatient = []
            # upre_or_post = []
            # for curr_patient, curr_pre_or_post in zip(patient, pre_or_post):
            #     if curr_patient not in upatient:
            #         upatient.append(curr_patient)
            #     if curr_pre_or_post not in upre_or_post:
            #         upre_or_post.append(curr_pre_or_post)

            # # Reshape the average holder for the main arrays to the individual patients and drug treatment, filling in "blanks" if necessary
            # impossible_val = 444.4
            # log_dens_pvals_avg2 = np.ones((len(upatient), len(upre_or_post), nall_species, nall_species, 2, nslices)) * impossible_val
            # log_pmf_pvals_avg2 = np.ones((len(upatient), len(upre_or_post), nall_species, nall_species, max_nbins, 2, nslices)) * impossible_val
            # for iunique_slide, (curr_patient, curr_pre_or_post) in enumerate(zip(patient, pre_or_post)):
            #     upatient_idx = upatient.index(curr_patient)
            #     upre_or_post_idx = upre_or_post.index(curr_pre_or_post)
            #     log_dens_pvals_avg2[upatient_idx, upre_or_post_idx, :, :, :, :] = log_dens_pvals_avg[iunique_slide, :, :, :, :]
            #     log_pmf_pvals_avg2[upatient_idx, upre_or_post_idx, :, :, :, :, :] = log_pmf_pvals_avg[iunique_slide, :, :, :, :, :]
            # log_dens_pvals_avg = log_dens_pvals_avg2
            # log_pmf_pvals_avg = log_pmf_pvals_avg2

            # # Create the Pandas indexes
            # index_density = pd.MultiIndex.from_product([upre_or_post, all_species_names, all_species_names, ['left', 'right'], np.arange(nslices) + 1], names=['drug_status', 'center_species', 'neighbor_species', 'pval_type', 'slice_number'])
            # index_pmf = pd.MultiIndex.from_product([upre_or_post, all_species_names, all_species_names, np.arange(max_nbins), ['left', 'right'], np.arange(nslices) + 1], names=['drug_status', 'center_species', 'neighbor_species', 'poss_num_neighbors', 'pval_type', 'slice_number'])

            # # For the density, Create the Pandas dataframes from the indexes
            # df_log_dens_pvals_avg = pd.DataFrame(data=log_dens_pvals_avg.reshape((len(upatient), -1)), index=upatient, columns=index_density).rename_axis('subject')
            # df_log_pmf_pvals_avg = pd.DataFrame(data=log_pmf_pvals_avg.reshape((len(upatient), -1)), index=upatient, columns=index_pmf).rename_axis('subject')

            # # Write the Pandas dataframes to disk
            # if write_csv_files:
            #     df_log_dens_pvals_avg.to_csv(density_csv_pathname)
            #     df_log_pmf_pvals_avg.to_csv(pmf_csv_pathname)

            # Save the averaged data to disk
            if write_pickle_datafile:
                # make_pickle((nmatches_holder, log_dens_pvals_avg, log_pmf_pvals_avg, df_log_dens_pvals_avg, df_log_pmf_pvals_avg), pickle_dir, pickle_file)
                make_pickle((nmatches_holder, log_dens_pvals_avg), pickle_dir, pickle_file)

            # Close the figure
            plt.close(fig)

        else:

            # Read in the averaged data from disk
            (nmatches_holder, log_dens_pvals_avg) = load_pickle(pickle_dir, pickle_file)

        # return(nmatches_holder, log_dens_pvals_avg, unique_slides)


    def plot_density_pvals_over_slides(self, roi_figsize=(15, 10), marker_size_step=0.80, alpha=0.1, default_marker_size_fac=1, roi_dpi=200, yaxis_dir=-1, edgecolors=None, filename_suffix2='', nworkers=1, use_multiprocessing=True, plot_real_data=True):
        # Possible reasons there could be no data for a ROI/center/neighbor combo: not having any of the center species present, having too few center species present, or to the ROI having a zero area

        # Import relevant libraries
        import numpy as np
        import os

        # Define variables already defined as attributes
        df = self.data  # get a shortcut for the main data dataframe
        plotting_map = self.plotting_map
        num_colors = self.num_colors
        webpage_dir = self.webpage_dir
        mapping_dict = self.mapping_dict
        coord_units_in_microns = self.coord_units_in_microns
        all_species_names = self.all_species_names
        log_pval_range = self.log_pval_range
        df_data_by_roi = self.df_data_by_roi
        unique_slides = self.unique_slides
        df_density_pvals_arrays = self.df_density_pvals_arrays

        # Constant
        zero = 1e-8

        # If the directory holding the density P values plotted over the slides does not exist, create and populate it
        if not os.path.exists(os.path.join(webpage_dir, 'density_pvals_over_slide_spatial_plot')):

            # # Create a directory for the images created herein
            # webpage_dir = os.path.join(webpage_dir, ('real' if plot_real_data else 'simulated'), 'density_pvals_over_slide_spatial_plot')
            # os.makedirs(webpage_dir, exist_ok=True)

            # Set the log of the P value range for the color plotting
            vmin = log_pval_range[0]
            vmax = log_pval_range[1]

            # Get rid of any overlap created by the patches so that the same cells aren't plotted multiple times, showing a gridded pattern
            depatched, _ = undo_patching_overlaps_and_decompounding(df)

            # Initialize the list that will hold one dataframe per slide whose rows will contain the coordinates of all the ROIs in that slide
            rectangle_data = []

            # For each slide...
            num_rois = 0  # just as a check for later
            for unique_slide in unique_slides:

                # Get the ROI data for the current slide
                curr_rectangle_data = df_data_by_roi[df_data_by_roi['unique_slide'] == unique_slide].copy()

                # Add a column containing a tuple of the minimum x- and y-coordinates, just as patches.Rectangle() requires
                curr_rectangle_data['xy'] = curr_rectangle_data[['x_roi', 'y_roi']].map(min).apply(tuple, axis='columns')

                # Calculate the area of every ROI
                curr_rectangle_data['area'] = curr_rectangle_data[['width', 'height']].apply(np.product, axis='columns')

                # Reformat the rectangle data dataframe
                curr_rectangle_data = curr_rectangle_data[['unique_roi', 'width', 'height', 'xy', 'area']].rename({'unique_roi': 'tag'}, axis='columns').set_index('tag').sort_index()

                # Drop ROIs of zero area
                len_with_zero_areas = len(curr_rectangle_data)
                curr_rectangle_data = curr_rectangle_data[curr_rectangle_data['area'] > zero]

                # Note whether ROIs of zero area exist
                num_dropped = len_with_zero_areas - len(curr_rectangle_data)
                if num_dropped != 0:
                    print('Note: {} ROIs of zero areas were found and excluded'.format(num_dropped))

                # Update the total number of ROIs
                num_rois = num_rois + len(curr_rectangle_data)

                # Save this ROI coordinate dataframe to the master list
                rectangle_data.append(curr_rectangle_data)

            # Print the total number of ROIs as a check
            print('A total of {} ROIs of non-zero area has been found in the plot_density_pvals_over_slides() method'.format(num_rois))

            # Generate the arguments to run in the parallel function by looping over all possible combinations of centers and neighbors
            rectangle_data_orig = rectangle_data.copy()
            list_of_tuple_arguments = []
            for icenter, center_name in enumerate(all_species_names):
                for ineighbor, neighbor_name in enumerate(all_species_names):
                    list_of_tuple_arguments.append((rectangle_data_orig, unique_slides, df_density_pvals_arrays, filename_suffix2, depatched, plotting_map, num_colors, webpage_dir, mapping_dict, coord_units_in_microns, roi_figsize, marker_size_step, alpha, edgecolors, default_marker_size_fac, roi_dpi, yaxis_dir, vmin, vmax, icenter, ineighbor, center_name, neighbor_name, 'density_pvals_over_slide_spatial_plot'))

            # Potentially farm out the P values plotted over the slides to the worker CPUs
            utils.execute_data_parallelism_potentially(function=plot_pvals_over_slides_for_center_neighbor_pair, list_of_tuple_arguments=list_of_tuple_arguments, nworkers=(0 if not use_multiprocessing else nworkers), task_description='plotting of the density P values at the ROI locations on top of all slides')

        # Print that we're not running this component
        else:
            print('Directory {} already exists; not plotting the density P values over the slides now'.format(os.path.join(webpage_dir, 'density_pvals_over_slide_spatial_plot')))


    def plot_dens_pvals_per_roi(self, nworkers=1, figsize=(10, 4), dpi=100, img_file_suffix='', plot_real_data=True, use_multiprocessing=True):

        # Import relevant libraries
        import os
        import subprocess

        # Define variables already defined as attributes
        webpage_dir = self.webpage_dir
        all_species_names = self.all_species_names
        # df_data_by_roi = self.df_data_by_roi
        log_pval_range = self.log_pval_range
        df_density_pvals_arrays = self.df_density_pvals_arrays
        num_valid_centers_minimum = self.num_valid_centers_minimum

        # Get the ROIs to plot here
        roi_indexes_with_any_heatmap_data_set = set(df_density_pvals_arrays.index)

        # Print the number of ROIs found in the data
        print('{} ROIs found with any heatmap data'.format(len(roi_indexes_with_any_heatmap_data_set)))

        # Define the directory to hold the plots of the density P values for each ROI and ensure this directory exists
        webpage_dir_for_dens_pvals_per_roi = os.path.join(webpage_dir, 'dens_pvals_per_roi')
        if not os.path.exists(webpage_dir_for_dens_pvals_per_roi):
            os.makedirs(webpage_dir_for_dens_pvals_per_roi)

        # Get a reasonable title suffix
        title_suffix = ' - min # valid centers: {}'.format(num_valid_centers_minimum)

        # Determine the ROIs whose density P values need to be plotted, i.e., those whose corresponding PNG files are not present on the filesystem
        retval = subprocess.run(['ls {}/density_pvals-*-slice_*-roi_index_*.png'.format(webpage_dir_for_dens_pvals_per_roi)], shell=True, capture_output=True)
        roi_ids_not_present = roi_indexes_with_any_heatmap_data_set - set([int(x.split('.png')[0].split('_')[-1]) for x in retval.stdout.decode().split('\n')[:-1]])

        # Generate a list of tuple arguments each of which is inputted into plot_single_density_pvals() to be run by a single worker
        # I am sending in df_data_by_roi now, which was totally unnecessary, but only to update the image filenames!! If things get weird or don't work, try removing this!!!!
        constant_tuple = (log_pval_range, figsize, webpage_dir_for_dens_pvals_per_roi, plot_real_data, img_file_suffix, all_species_names, dpi, df_density_pvals_arrays, title_suffix)
        list_of_tuple_arguments = [constant_tuple + (x,) for x in roi_ids_not_present]

        # Potentially farm out the metrics calculations to the worker CPUs. This ensures that a pickle file gets created for each ROI
        utils.execute_data_parallelism_potentially(function=plot_single_density_pvals, list_of_tuple_arguments=list_of_tuple_arguments, nworkers=(0 if not use_multiprocessing else nworkers), task_description='plotting of the density heatmaps')

    def incorporate_phenotype_identifications(self, tsv_file):
        """Take hyphen-separated phenotypes ("reduced_marker_set" below) specified by biologists in Excel for sets of surface markers exported from the plotting map and incorporate these identifications into the main data object

        Workflow:
          (1) Run tci.TIMECellInteraction() as usual
          (2) Export the plotting map using utils.export_plotting_map()
          (3) Ask biologists to fill out the resulting CSV file in Excel in order to identify the phenotypes
          (4) Delete the .pkl files created so far
          (5) Run the entire pipeline again from scratch, this time adding, e.g., "phenotype_identification_tsv_file='gmb_phenotype_ids_as_of_2022-06-06.tsv'" in the arguments to tci.preprocess_dataset()

        Originally taken from utilities/utils.py.

        Args:
            tsv_file (str): Name of the TSV file generated by pasting into a text editor columns from Excel in the order specified in the read_csv() method below
        """

        # Define variable already defined as attributes
        df_objects_orig = self.data

        # Import relevant libraries
        import pandas as pd
        import numpy as np

        # Define a hardcoded constant that in the future I can put in the same place as phenotype_identification_tsv_file='gmb_phenotype_ids_as_of_2022-06-06.tsv'
        species_id_col = 'Species int'

        # Read in the data from the filled-out-by-collaborators TSV file specifying the phenotype identifications from the surface markers, setting the 'species_id' as both the index and keeping the values in a column, dropping the unnecessary "species_perc" column
        df_phenotypes = pd.read_csv(tsv_file, sep='\t', names=['species_id', 'species_count', 'species_perc', 'positive_markers', 'reduced_marker_set']).set_index('species_id', drop=False).drop('species_perc', axis='columns')

        # Initialize a dataframe containing just the final set of species, using as initial data the "first" rows of the phenotype data grouped by species
        df_reduced_species = df_phenotypes[df_phenotypes['reduced_marker_set'].apply(lambda x: '-' not in x)].groupby(by='reduced_marker_set').first()

        # Add a column to the phenotypes dataframe with the species IDs to which we want to copy the object data
        df_phenotypes['ids_to_map'] = df_phenotypes['reduced_marker_set'].apply(lambda x: [df_reduced_species['species_id'].loc[y] for y in x.split('-')])

        # Clean up the dataframe of the final phenotypes, dropping the "positive_markers" column (since it only contains one of multiple possible sets), initializing the species counts to zero, and setting the index to the species_id, which is unique though corresponds to the single set of positive markers, which again is only one of multiple possible sets
        df_reduced_species = df_reduced_species.drop('positive_markers', axis='columns')
        df_reduced_species['species_count'] = 0
        df_reduced_species = df_reduced_species.reset_index().set_index('species_id')

        # Store the value counts of the original set of species in the original dataset
        orig_species_value_counts = df_objects_orig[species_id_col].value_counts()

        # Check that the total number of objects in the original dataset is consistent using multiple measures. Note that this can get flagged if the dataset is inconsistent with the TSV file used for phenotype identification, which hardcodes some properties of the original dataset, like the total number of cells of each naiive phenotype.
        tot_number_of_objects_orig = [orig_species_value_counts.sum(), df_phenotypes['species_count'].sum(), len(df_objects_orig)]
        if len(set(tot_number_of_objects_orig)) != 1:
            print('WARNING: Inconsistent total number of original objects! ({})'.format(tot_number_of_objects_orig))
            # exit()

        # Initialize some variables
        size_of_orig_data = 0  # stores total number of objects in original dataset
        size_of_new_data = 0  # stores total number of objects in new dataset
        list_objects_new = []  # stores dataframes to be combined into the new dataset

        # For each row in the original set of phenotypes...
        for row in df_phenotypes.iterrows():

            # Get the current species ID and phenotype data for that ID, in particular, the species IDs to which to copy the current set of object data
            species_id, row_data = row
            ids_to_map = row_data['ids_to_map']

            # Print out what we're doing
            print('Now copying data for objects with positive surface markers {} to the phenotype(s) {}'.format(row_data['positive_markers'], row_data['reduced_marker_set']))

            # Get a boolean series of the original object data that correspond to the current species ID
            object_is_species = df_objects_orig[species_id_col] == species_id

            # Determine the total number of matches for the current species ID
            curr_num_objects = object_is_species.sum()

            # Check that the current number of object matches is consistent with what was previously calculated using the value_counts() method
            if species_id not in orig_species_value_counts.index.to_list():
                num_objects_in_value_counts = 0
            else:
                num_objects_in_value_counts = orig_species_value_counts[species_id]
            if curr_num_objects != num_objects_in_value_counts:
                print('ERROR: Current number of species matches is inconsistent! ({}, {})'.format(curr_num_objects, orig_species_value_counts[species_id]))
                exit()

            # Store a copy of the object data for the current species ID... copy is indeed needed or else I'd get a typical Pandas warning!
            curr_object_data = df_objects_orig[object_is_species].copy()

            # Update the total number of original and new objects
            size_of_orig_data = size_of_orig_data + curr_num_objects
            size_of_new_data = size_of_new_data + curr_num_objects * len(ids_to_map)

            # For each species ID to which to map the current object data...
            for id_to_map in ids_to_map:

                # Print out what we're doing
                print('  surface markers {} (ID {}) --> phenotype {} (ID {}) '.format(row_data['positive_markers'], species_id, df_phenotypes.loc[id_to_map, 'reduced_marker_set'], id_to_map))

                # Update the total numbers of the new, mapped species
                df_reduced_species.loc[id_to_map, 'species_count'] = df_reduced_species.loc[id_to_map, 'species_count'] + curr_num_objects

                # Overwrite the species IDs for the current data dataframe with the new species ID
                curr_object_data[species_id_col] = id_to_map

                # Save the current object data with updated species ID to the dataframe holder that we will later combine into the final dataset... note the .copy() is absolutely needed or else when there are multiple phenotypes per set of surface markers, copies of the dataframes are added to list_objects_new instead of ones with different species IDs!!
                list_objects_new.append(curr_object_data.copy())

        # Convert the Case ID column to integer so it's sorted correctly
        df_objects_new = pd.concat(list_objects_new)
        df_objects_new['Case ID'] = df_objects_new['Case ID'].astype(int)

        # Combine the list holding all the new object data into a single new dataframe
        df_objects_new = df_objects_new.sort_values(by=['Case ID', 'Slide ID', 'tag', 'Cell X Position', 'Cell Y Position', species_id_col]).reset_index(drop=True)

        # Sort the new species holder by decreasing frequency in the entire dataset
        df_reduced_species = df_reduced_species.sort_values(by='species_count', ascending=False)

        # Run checks on the sizes of the old and new datasets
        if size_of_orig_data != len(df_objects_orig):
            print('ERROR: Inconsistent total number of original objects (second check)! ({}, {})'.format(size_of_orig_data, len(df_objects_orig)))
            exit()
        if size_of_new_data != len(df_objects_new):
            print('ERROR: Inconsistent total number of new objects! ({}, {})'.format(size_of_new_data, len(df_objects_new)))
            exit()

        # Visually compare the species counts calculated two different ways
        print('NOTE: Visually ensure the following two structures are the same (automate this in the future)!')
        print(df_reduced_species)
        print(df_objects_new[species_id_col].value_counts())

        # Set the number of colors for plotting the ROIs in this case as the total number of final phenotypes
        num_colors = len(df_reduced_species)

        # Create the plotting map from df_reduced_species
        df_plotting_map = df_reduced_species.reset_index(drop=False).set_index('species_id', drop=False)
        df_plotting_map['reduced_marker_set'] = df_plotting_map['reduced_marker_set'].apply(lambda x: [x])
        df_plotting_map['marker_ids'] = range(num_colors)
        # df_plotting_map['circle_sizes'] = [x + 1 for x in range(num_colors)]
        df_plotting_map['circle_sizes'] = 1
        plotting_map = df_plotting_map.to_numpy()

        # Store the unique species IDs as a numpy array
        unique_species = df_plotting_map['species_id'].to_numpy(dtype=np.int64)

        # Save the calculated data as properties of the class object; actually as this doesn't seem to be working, just return these four variables and overwrite them in the class instantiation; I may have run into this issue before
        # self.data = df_objects_new
        # self.num_colors = num_colors
        # self.plotting_map = plotting_map
        # self.unique_species = unique_species
        return df_objects_new, num_colors, plotting_map, unique_species


    def flatten_density_pvals(self):
        """Return a Pandas dataframe of the heatmap data (the log of the density P values) flattened to a single dataframe, just like is obtained from slices.filedata. This can be directly compared with slices.filedata using:

        pd.DataFrame(data=slices.filedata).drop(['roi_fig_pathname', 'pvals_fig_pathname', 'left_log_pmf_pvals', 'right_log_pmf_pvals'], axis='columns').equals(slices.flatten_density_pvals())

        Also, output a log describing what data are *not* present in slices.metrics, which is useful for debugging.

        See utils.create_filedata_from_metrics_log() for some more details of the comments below that describe the two different situations when None values are obtains (i.e., when (1) either the centers or neighbors do not exist in the ROI, or (2) they both exist but there are no *valid* centers).

        This function adds rows to df_density_pvals if both a center and neighbor exist and there is at least 1 *valid* center.

        Returns:
            Pandas dataframe: Dataframe holding the bare minimum heatmap data.
        """

        # Import relevant libraries
        import numpy as np
        import pandas as pd
        import os

        # Define variables already defined as attributes
        metrics = self.metrics
        plotting_map = self.plotting_map
        mapping_dict = self.mapping_dict

        # Print what we're doing
        print('Flattening the calculated metrics into a single dataframe...')

        # Get the full list of phenotypes in the dataset in decreasing frequency order
        plotting_map_species = [x[0] for x in plotting_map]  # same as all_species_list previously
        num_all_species = len(plotting_map_species)  # same as nall_species previously
        range_num_all_species = range(num_all_species)

        # Get an array of the species names in decreasing frequency order
        # all_species_list = [x[0] for x in plotting_map]
        # nall_species = len(all_species_list)
        if mapping_dict is not None:
            species_names = [get_descriptive_cell_label(x[1], mapping_dict)[0] for x in plotting_map]
        else:
            species_names = [phenotypes_to_string(x[1]) for x in plotting_map]

        # Initialize the main file data holder that will be concatenated into a Pandas dataframe
        filedata_list = []

        # Open a log file for writing notes about the data as we add it to the filedata structure
        logs_dir = os.path.join('.', 'output', 'logs')
        if not os.path.exists(logs_dir):
            os.makedirs(logs_dir)
        with open(os.path.join(logs_dir, 'metrics_check.log'), 'wt') as f:

            num_slides = len(metrics)

            # For every slide...
            progress_int_old = 0
            for islide, (curr_slide, unique_rois, data_by_roi) in enumerate(metrics):  # for each slide of metrics data...

                # Determine the number of ROIs in the current slide
                nrois_in_slide = len(unique_rois)  # get the number of ROIs in the current slide; I checked that this len always equals that of data_by_roi but I can always add an assert statement to make sure

                # For every ROI in the current slide...
                for iroi, (roi_name, roi_data) in enumerate(zip(unique_rois, data_by_roi)):

                    # Don't do tuple unpacking in the loop definition but do it here so we can explicitly define roi_data to save for later
                    roi_min_coord_spacing = roi_data['roi_min_coord_spacing']
                    unique_species_in_roi = roi_data['unique_species_in_roi']
                    real_data = roi_data['real_data']
                    all_species_list = roi_data['all_species_list']
                    roi_x_range = roi_data['roi_x_range']
                    roi_y_range = roi_data['roi_y_range']

                    # For every center-neighbor combination in the current ROI...
                    for icenter_spec in range_num_all_species:
                        for ineighbor_spec in range_num_all_species:

                            # Get an ID string for the current slide/ROI/center/neighbor
                            center_species_id = all_species_list[icenter_spec]
                            neighbor_species_id = all_species_list[ineighbor_spec]
                            current_id = 'slide {}, roi {}, center_spec {}, neighbor_spec {}'.format(curr_slide, roi_name, center_species_id, neighbor_species_id)

                            # Get the current dataset
                            curr_real_data = real_data[icenter_spec, ineighbor_spec, 0]

                            # If there is no current dataset (i.e., when either the center or neighbor species does not exist in the ROI), say so
                            if curr_real_data is None:
                                f.write('NOTE: No real data exist for {}\n'.format(current_id))

                            # If there currently is a real dataset (i.e., if both the center and neighbor species are in the ROI)...
                            else:

                                # Extract the individual pieces from the current real dataset
                                density_metrics_real, pmf_metrics_real, nexpected_real, nvalid_centers_real, coords_centers_real, coords_neighbors_real, valid_centers_real, edges_real, npossible_neighbors_real, roi_area_used_real, slice_area_used_real, center_species, neighbor_species, small_rad, large_rad, islice = curr_real_data

                                # If there is no current set of density metrics (i.e., there are no *valid* centers), say so
                                if density_metrics_real is None:
                                    f.write('NOTE: Real data exist but the density metrics do not for {}\n'.format(current_id))

                                # If there currently is a set of density metrics (i.e., there are *some* valid centers)...
                                else:

                                    # Extract the individual pieces from the current density metrics
                                    z_score, left_p_val, right_p_val = density_metrics_real

                                    # Print a warning if either the left or right P values do not exist (this never occurs)
                                    if (left_p_val is None) or (right_p_val is None):
                                        f.write('NOTE: Density metrics exist but the left ({}) and/or right ({}) P values are None for {}\n'.format(left_p_val, right_p_val, current_id))

                                    # If both the left and right P values exist, which at this point they always do, add all relevant data to the filedata structure
                                    else:
                                        filedata_list.append({'nrois_in_slide': nrois_in_slide, 'roi_name': roi_name, 'unique_species_in_roi': list(filter(lambda x: x in unique_species_in_roi, plotting_map_species)), 'roi_x_range': np.array(roi_x_range), 'roi_y_range': np.array(roi_y_range), 'roi_spacing': roi_min_coord_spacing, 'center_species_id': center_species_id, 'neighbor_species_id': neighbor_species_id, 'center_species_name': species_names[icenter_spec], 'neighbor_species_name': species_names[ineighbor_spec], 'nvalid_centers_per_slice': [nvalid_centers_real], 'left_log_dens_pvals': [[np.log10(left_p_val)]], 'right_log_dens_pvals': [[np.log10(right_p_val)]]})

                    # Output flattening progress
                    progress_int = round((iroi + 1) / nrois_in_slide * (islide + 1) / num_slides * 100)
                    if progress_int != progress_int_old:
                        print('{:d}% done flattening calculated metrics'.format(progress_int), end="\r")
                        progress_int_old = progress_int

        # Return a Pandas dataframe from the filedata structure
        self.df_density_pvals = pd.DataFrame(data=filedata_list)

        # Define some items we used to pass using "metadata"
        self.all_species_names = species_names
        self.all_species_ids = plotting_map_species
        self.nall_species = num_all_species


    def flatten_roi_plotting_data(self):

        # Import relevant libraries
        import numpy as np
        import pandas as pd

        # Define variables already defined as attributes
        df_data_by_roi = self.df_data_by_roi
        data = self.data
        plotting_map = self.plotting_map

        print('Plotting Map:')
        print(self.plotting_map)

        # Initialize lists that will become columns of the main dataframe
        roi_index_holder = []
        roi_name_holder = []
        slide_name_holder = []
        num_rois_in_slide_holder = []
        x_roi_holder, y_roi_holder, x_range_holder, y_range_holder, species_roi_holder, spec2plot_roi_holder = [], [], [], [], [], []

        # For every unique slide...
        for slide_name in df_data_by_roi['unique_slide'].unique():

            # Determine the unique ROIs in the current slide
            unique_rois = df_data_by_roi[df_data_by_roi['unique_slide'] == slide_name]['unique_roi'].unique()

            # For every unique ROI in the current slide...
            for roi_name in unique_rois:

                # Note that we're performing flattening on the current ROI
                print('Flattening ROI {}...'.format(roi_name))

                # Get specific data for the current ROI
                # The following block is based on calculate_metrics_for_roi()
                data_roi = data.iloc[np.nonzero((data['tag'] == roi_name).to_numpy())[0]].reset_index(drop=True)
                x_roi = np.array(data_roi['Cell X Position'])
                y_roi = np.array(data_roi['Cell Y Position'])
                roi_x_range, roi_y_range, _ = roi_checks_and_output(x_roi, y_roi, do_printing=False, do_main_printing=False)
                x_range = np.array(roi_x_range)
                y_range = np.array(roi_y_range)
                species_roi = np.array(data_roi['Species int'], dtype='uint64')
                spec2plot_roi = [x[0] for x in plotting_map if x[0] in np.unique(species_roi)]

                # Use the "df_data_by_roi" dataframe to determine the df_data_by_roi index of the current ROI
                roi_index = df_data_by_roi.loc[df_data_by_roi['unique_roi'] == roi_name, :].index[0]

                # Check consistency of slide names
                assert slide_name == df_data_by_roi.loc[roi_index, 'unique_slide'], 'ERROR: slide name is inconsistent'

                # Append variables defined above into the running lists
                roi_index_holder.append(roi_index)
                roi_name_holder.append(roi_name)
                slide_name_holder.append(slide_name)
                num_rois_in_slide_holder.append(len(unique_rois))  # get the number of ROIs in the current slide
                x_roi_holder.append(x_roi)
                y_roi_holder.append(y_roi)
                x_range_holder.append(x_range)
                y_range_holder.append(y_range)
                species_roi_holder.append(species_roi)
                spec2plot_roi_holder.append(spec2plot_roi)

        # Create a Pandas dataframe holding all the per-ROI data so far
        # Columns already in the df: unique_roi, ncases, ncases_tot, unique_case, nslides, nslides_tot, unique_slide, nrois, nrois_tot, width, height, x_min_prior_to_decimation, y_min_prior_to_decimation, x_max_prior_to_decimation, y_max_prior_to_decimation
        df_data_by_roi = pd.concat(
            [
                df_data_by_roi.sort_index(),
                pd.DataFrame({'roi_name': roi_name_holder, 'slide_name': slide_name_holder, 'num_rois_in_slide': num_rois_in_slide_holder, 'x_roi': x_roi_holder, 'y_roi': y_roi_holder, 'x_range': x_range_holder, 'y_range': y_range_holder, 'species_roi': species_roi_holder, 'spec2plot_roi': spec2plot_roi_holder}, index=roi_index_holder).sort_index()  # new columns, but after checks just below roi_name and slide_name should be dropped
            ], axis='columns'
        )

        # Run some checks
        if df_data_by_roi['roi_name'].equals(df_data_by_roi['unique_roi']):
            print('roi_name and unique_roi are equal; dropping roi_name')
            df_data_by_roi = df_data_by_roi.drop('roi_name', axis='columns')
        else:
            print('NOTE: roi_name and unique_roi are not equal!')
        if df_data_by_roi['slide_name'].equals(df_data_by_roi['unique_slide']):
            print('slide_name and unique_slide are equal; dropping slide_name')
            df_data_by_roi = df_data_by_roi.drop('slide_name', axis='columns')
        else:
            print('NOTE: slide_name and unique_slide are not equal!')
        if df_data_by_roi['x_range'].apply(lambda x: x[1] - x[0]).equals(df_data_by_roi['width']):
            print('width is consistent')
        else:
            print('WARNING: width is INconsistent')
            assert (df_data_by_roi['x_max_prior_to_decimation'] - df_data_by_roi['x_min_prior_to_decimation']).equals(df_data_by_roi['width']), 'ERROR: it''s REALLY inconsistent!'
            print('...that''s likely due to decimation... at least the original range/width is consistent!')
        if df_data_by_roi['y_range'].apply(lambda x: x[1] - x[0]).equals(df_data_by_roi['height']):
            print('height is consistent')
        else:
            print('WARNING: height is INconsistent')
            assert (df_data_by_roi['y_max_prior_to_decimation'] - df_data_by_roi['y_min_prior_to_decimation']).equals(df_data_by_roi['height']), 'ERROR: it''s REALLY inconsistent!'
            print('...that''s likely due to decimation... at least the original range/height is consistent!')

        # Return the variable of interest
        return df_data_by_roi


    def average_dens_pvals_over_rois_for_each_slide(self, figsize=(10, 4), dpi=100, img_file_suffix='', plot_real_data=True, do_plotting=True, weight_rois_by_num_valid_centers=False, input_datafile='../data/dummy_txt_or_csv.csv'):
        """Average the P values over all ROIs for each slide. Note I have confirmed that this yields the same results as a non-weighted version of the original method.

        Args:
            figsize (tuple, optional): P value figure size in inches. Defaults to (14, 4).
            dpi (int, optional): DPI of the figure. Defaults to 150.
            img_file_suffix (str, optional): String you want to append to the image names. Defaults to ''.
            plot_real_data (bool, optional): Whether to plot the real or simulated data. This is likely meaningless now but was useful historically. Defaults to True.
            do_plotting (bool, optional): Whether to save the average P value plots in addition to calculating the average P values. Defaults to True.
        """

        # Import relevant libraries
        import numpy as np
        import pandas as pd
        import os
        import numpy.ma as ma

        # Define variables from object properties
        df_density_pvals_arrays = self.df_density_pvals_arrays
        df_data_by_roi = self.df_data_by_roi
        num_valid_centers_minimum = self.num_valid_centers_minimum
        log_pval_range = self.log_pval_range
        webpage_dir = self.webpage_dir
        all_species_names = self.all_species_names
        pickle_dir = self.pickle_dir

        # Constants
        entity = 'slide'
        pickle_file = 'data_for_input_into_correlation_analyzer.pkl'

        # Extract the basename of the input datafile
        input_datafile_basename = os.path.basename(input_datafile).split(os.path.extsep)[0]

        # Determine the slides in the full dataset
        unique_slides = df_data_by_roi['unique_slide'].unique()

        # Determine some properties of the final array from a sample log_dens_pvals_arr array
        sample_log_dens_pvals_arr_idx = df_density_pvals_arrays[df_density_pvals_arrays['log_dens_pvals_arr'].apply(lambda x: type(x) == np.ndarray)]['log_dens_pvals_arr'].index[0]
        sample_num_valid_centers_idx = df_density_pvals_arrays[df_density_pvals_arrays['num_valid_centers'].apply(lambda x: type(x) == np.ndarray)]['num_valid_centers'].index[0]
        sample_log_dens_pvals_arr = df_density_pvals_arrays.loc[sample_log_dens_pvals_arr_idx, 'log_dens_pvals_arr']
        sample_num_valid_centers = df_density_pvals_arrays.loc[sample_num_valid_centers_idx, 'num_valid_centers']
        constant_shape_log_dens_pvals_arr = sample_log_dens_pvals_arr.shape
        constant_shape_num_valid_centers = sample_num_valid_centers.shape
        pvals_dtype = sample_log_dens_pvals_arr.dtype
        num_centers_dtype = sample_num_valid_centers.dtype

        # Print what we're doing
        print('  * Minimum number of valid centers: {}'.format(num_valid_centers_minimum))
        print('  * Logarithm of the P value range: {}'.format(log_pval_range))
        print('  * Weighting the ROIs by the number of valid centers: {}'.format(weight_rois_by_num_valid_centers))

        # Generate an image title suffix with plotting/analysis settings
        title_suffix = ' - ROI weighting: {}, min # valid centers: {}'.format(weight_rois_by_num_valid_centers, num_valid_centers_minimum)

        # Get a list of the indexes of all ROIs having filtered (e.g., minimum number of valid centers) P value data
        roi_indexes_with_pval_data = df_density_pvals_arrays.index.to_list()

        # Initialize the array that will hold the numpy array averages for each slide
        log_dens_pvals_arr_per_slide = []

        # For every slide in the dataset, regardless of whether any valid P value data actually exists...
        for slide_name in unique_slides:

            # Determine the indexes of all possible ROIs in the slide
            roi_indexes_in_slide = df_data_by_roi[df_data_by_roi['unique_slide'] == slide_name].index.to_list()

            # Filter these out based on whether actual valid P value data exist and determine how many such valid ROIs there are
            roi_indexes_found = [x for x in roi_indexes_in_slide if x in roi_indexes_with_pval_data]  # could have used set intersection
            num_rois_found = len(roi_indexes_found)

            # Print out how many ROIs have been found for the current slide
            print('{} ROIs (indexes {}) have been found for slide {}'.format(num_rois_found, roi_indexes_found, slide_name))

            # If there's at least one ROI for the current slide...
            if num_rois_found >= 1:

                # Initialize the array holding the data for all the ROIs for the current slide
                log_dens_pvals_arr_all_rois_in_slide = np.zeros(constant_shape_log_dens_pvals_arr + (num_rois_found,), dtype=pvals_dtype)
                num_valid_centers_all_rois_in_slide = np.zeros(constant_shape_num_valid_centers + (num_rois_found,), dtype=num_centers_dtype)

                # For every valid ROI for the current slide, add that ROI's P value data to the overall array
                for iroi, (log_dens_pvals_arr, num_valid_centers) in enumerate(zip(df_density_pvals_arrays.loc[roi_indexes_found, 'log_dens_pvals_arr'], df_density_pvals_arrays.loc[roi_indexes_found, 'num_valid_centers'])):
                    log_dens_pvals_arr_all_rois_in_slide.reshape((-1, num_rois_found))[:, iroi] = log_dens_pvals_arr.flatten()
                    num_valid_centers_all_rois_in_slide.reshape((-1, num_rois_found))[:, iroi] = num_valid_centers.flatten()

                # Average this array holding all the ROI data over all the ROIs
                if weight_rois_by_num_valid_centers:
                    num_valid_centers_all_rois_in_slide = np.tile(num_valid_centers_all_rois_in_slide[:, :, np.newaxis, :, :], (1, 1, 2, 1, 1))
                else:
                    num_valid_centers_all_rois_in_slide = np.ones(log_dens_pvals_arr_all_rois_in_slide.shape)

                # Perform weighted averaging using the number of valid centers as weights while taking into account np.nan values in the log of the density P values
                nan_locations = np.isnan(log_dens_pvals_arr_all_rois_in_slide)
                log_dens_pvals_arr_all_rois_in_slide = ma.array(log_dens_pvals_arr_all_rois_in_slide, mask=nan_locations)
                num_valid_centers_all_rois_in_slide = ma.array(num_valid_centers_all_rois_in_slide, mask=nan_locations)
                weights = num_valid_centers_all_rois_in_slide / num_valid_centers_all_rois_in_slide.sum(axis=-1, keepdims=True)  # .data add this to the end if things don't work
                # log_dens_pvals_arr_per_slide.append(np.nanmean(log_dens_pvals_arr_all_rois_in_slide, axis=-1))  # good way to do it if not using weighting and masked array module

                # Need the following three lines (this involving "ans") or else all-nan slices seem to get filled in by zeros
                ans = (log_dens_pvals_arr_all_rois_in_slide * weights).sum(axis=-1)
                ans.set_fill_value(np.nan)
                log_dens_pvals_arr_per_slide.append(ans.filled())

            # If there is no data for the current slide, just return None
            else:
                log_dens_pvals_arr_per_slide.append(None)

        # Convert the slide P values to a Pandas dataframe and save this as a property of the object
        df_log_dens_pvals_arr_per_slide = pd.DataFrame(data={'log_dens_pvals_arr': log_dens_pvals_arr_per_slide}, index=unique_slides)
        self.df_log_dens_pvals_arr_per_slide = df_log_dens_pvals_arr_per_slide

        # Write the calculated data to disc for subsequent loading into correlation analyzer
        if not os.path.exists(os.path.join(pickle_dir, pickle_file)):
            make_pickle((df_log_dens_pvals_arr_per_slide, input_datafile_basename), pickle_dir, pickle_file)

        # Define the directory holding all the images of the averaged data
        savedir = os.path.join(webpage_dir, 'dens_pvals_per_{}'.format(entity))

        # If the directory holding the slide-based density P values does not exist, create and populate it
        if not os.path.exists(savedir):

            # Create the directory holding all the images of the averaged data
            os.makedirs(savedir)

            # Plot the averaged P values for every slide
            for slide_name, srs_slide_data in df_log_dens_pvals_arr_per_slide.iterrows():
                if srs_slide_data['log_dens_pvals_arr'] is not None:
                    print('Plotting P values for slide {} averaged over all its ROIs with valid P value data'.format(slide_name))
                    plot_density_pvals_simple(
                        log_dens_pvals_arr=srs_slide_data['log_dens_pvals_arr'],
                        log_pval_range=log_pval_range,
                        figsize=figsize,
                        dpi=dpi,
                        plots_dir=savedir,
                        plot_real_data=plot_real_data,
                        entity_name=slide_name,
                        img_file_suffix=img_file_suffix,
                        entity=entity,
                        entity_index=-1,
                        all_species_names=all_species_names,
                        title_suffix=title_suffix
                    )
                else:
                    print('Not plotting P values for slide {} averaged over all its ROIs with valid P value data, because there are no ROIs with valid P value data'.format(slide_name))

        # Print that we're not running the plotting part of this component
        else:
            print('Directory {} already exists; not plotting slide-based density heatmaps now'.format(savedir))


    def check_and_prepare_metrics_for_plotting(self, num_valid_centers_minimum=1, log_pval_range=(-50, 0), correct_flooring=True):

        # Import relevant libraries
        import os
        import subprocess
        import pandas as pd

        # Define variables from object properties
        pickle_dir = self.pickle_dir
        df_data_by_roi = self.df_data_by_roi
        nall_species = self.nall_species
        nslices = self.nslices
        df_density_pvals = self.df_density_pvals
        all_species_ids = self.all_species_ids

        # New pickle filename for the checked and put-into-array-form metrics
        pickle_file = 'density_pvals_arrays.pkl'

        # If the pickle file doesn't already exist...
        if not os.path.exists(os.path.join(pickle_dir, pickle_file)):

            # ---- Check and put into array form the metrics for all the ROIs that haven't already been processed, saving the results in individual pickle files

            # Determine the ROIs and indexes that are present in the metrics data and that have at least one center species with a minimum number of valid centers
            max_num_valid_centers_per_roi = df_density_pvals.groupby(by='roi_name')['nvalid_centers_per_slice'].agg(lambda x: x.max()[0])
            rois_with_at_least_one_valid_heatmap_cell = max_num_valid_centers_per_roi[max_num_valid_centers_per_roi >= num_valid_centers_minimum].index  # this prevents plotting of blank heatmaps per ROI by sending in to generate_dens_pvals_array_for_roi() only ROIs having at least one center species with a minimum number of valid centers. The averaging and subsequent plotting of heatmaps for slide averaged over ROI would be unaffected (since ROIs only contribute to the average if they contain center-neighbor data); this line only prevents blank per-ROI heatmaps. In other words, doing this prevents df_density_pvals_arrays from having all nan values for a ROI in the already reduced set of ROIs, which could otherwise occur if the number of valid centers for all species is less than the minimum
            index_holder = []
            for roi_name in rois_with_at_least_one_valid_heatmap_cell:
                df_one_row = df_data_by_roi[df_data_by_roi['unique_roi'] == roi_name]
                assert len(df_one_row) == 1, 'ERROR: Not one ({}) ROI with roi_name {}'.format(len(df_one_row), roi_name)
                index_holder.append(df_one_row.index[0])
            index_holder_set = set(index_holder)
            assert len(index_holder) == len(index_holder_set), 'ERROR: Not purely unique ROIs were determined while processing df_density_pvals'

            # Determine the ROIs that need to be processed, i.e., those whose corresponding pickle files are not present on the filesystem
            retval = subprocess.run(['ls {}/dens_pvals_array-roi_index_*.pkl'.format(pickle_dir)], shell=True, capture_output=True)
            roi_ids_not_present = index_holder_set - set([int(x.split('.pkl')[0].split('_')[-1]) for x in retval.stdout.decode().split('\n')[:-1]])

            # Generate a list of tuple arguments each of which is inputted into calculate_metrics_for_roi() to be run by a single worker
            constant_tuple = (pickle_dir, nall_species, all_species_ids, nslices, num_valid_centers_minimum, log_pval_range, False, df_data_by_roi, df_density_pvals, correct_flooring)
            list_of_tuple_arguments = [constant_tuple + (x,) for x in roi_ids_not_present]

            # Farm out the metrics calculations to the worker CPUs. This ensures that a pickle file gets created for each ROI
            # print('Running {} function calls using {} workers'.format(len(list_of_tuple_arguments), nworkers))
            print('Running {} function calls using {} workers (this is forced to be {} in check_and_prepare_metrics_for_plotting())'.format(len(list_of_tuple_arguments), 1, 1))
            # with mp.get_context('spawn').Pool(nworkers) as pool:  # This weirdly seems to use only 2 workers when 31 are requested and takes ~14 min total. Note I actually think 31 workers are actually being used (albeit not reported on the node or in the HPC dashboard) because when I tell the single function to sleep for five seconds at the end, I seem to get about 31 ROI plots generated every 5 seconds. This implies that the single function takes almost negligible time, which may explan why using nworkers=1 (~4 min total) seems to be faster than nworkers=31: it actually may be that each instance completes so quickly that the overhead of parallelism is detrimental. Note that when I sleep for five seconds in the single function the single calls seem to call the ROIs in a random order, but when I do not, the single calls seem to call the ROIs in order, which may corroborate the each-instance-completing-so-quickly-that-paralleism-messes-up theory here.
            utils.execute_data_parallelism_potentially(function=generate_dens_pvals_array_for_roi, list_of_tuple_arguments=list_of_tuple_arguments, nworkers=0, task_description='checking and preparing the metrics for plotting')

            # ---- Load the resulting individual pickle files into a new, single pickle file called density_pvals_arrays.pkl

            # For each slide...
            data_by_roi = []
            roi_pickle_files = []  # store filenames to delete later
            for roi_index in index_holder:
                print('Reading ROI {}...'.format(roi_index))

                # Load the appropriate pickle file
                roi_pickle_file = 'dens_pvals_array-roi_index_{:06}.pkl'.format(roi_index)
                roi_pickle_files.append(roi_pickle_file)
                roi_data_item = load_pickle(pickle_dir, roi_pickle_file)

                # Save the loaded data
                data_by_roi.append(roi_data_item)

            # Create a dataframe out of the list of dictionaries
            df_density_pvals_arrays = pd.DataFrame(data=data_by_roi, index=index_holder).rename_axis('roi_index', axis='index')

            # Create the single pickle file saving all the data
            make_pickle((df_density_pvals_arrays, num_valid_centers_minimum, log_pval_range), pickle_dir, pickle_file)

            # If the overall pickle file was successfully created, delete all intermediate pickle files for the ROIs
            if os.path.exists(os.path.join(pickle_dir, pickle_file)):
                for roi_pickle_file in roi_pickle_files:
                    os.remove(os.path.join(pickle_dir, roi_pickle_file))

        # If the pickle file already exists, load it
        else:
            df_density_pvals_arrays, num_valid_centers_minimum, log_pval_range = load_pickle(pickle_dir, pickle_file)

        # Save the calculated data as properties of the class object
        self.df_density_pvals_arrays = df_density_pvals_arrays
        self.num_valid_centers_minimum = num_valid_centers_minimum
        self.log_pval_range = log_pval_range


    def plot_whole_slide_patches(self, roi_figsize=(15, 10), depatch=True, yaxis_dir=-1):
        """Plot the outlines of every ROI on top of the slides

        I moved this from utils.py to here on 4/26/23.

        Sample call: slices.plot_whole_slide_patches(roi_figsize=(15, 10), depatch=False)

        Args:
            roi_figsize (tuple, optional): Size (width, height) in inches of the output figures. Defaults to (15, 10).
            depatch (bool, optional): Whether to perform depatching first. Defaults to True.
        """

        # Import relevant library
        import os

        # Define variables already defined as attributes
        df = self.data  # get a shortcut for the main data dataframe
        plotting_map = self.plotting_map
        num_colors = self.num_colors
        mapping_dict = self.mapping_dict
        coord_units_in_microns = self.coord_units_in_microns
        webpage_dir = self.webpage_dir

        # If the directory holding the whole slide patches does not exist, create and populate it
        if not os.path.exists(os.path.join(webpage_dir, 'whole_slide_patches')):

            # Get rid of any overlap created by the patches so that the same cells aren't plotted multiple times, showing a gridded pattern
            if depatch:
                depatched, _ = undo_patching_overlaps_and_decompounding(df)
            else:
                depatched = df

            # For kicks, first just plot the ROIs with no grid lines overplotted, which is the main point of this function
            # Note that roi_id_col='Slide ID' is set that way because we want to treat the entire slide as a single ROI, both here and when this function is also called below
            plot_just_rois(depatched, plotting_map, num_colors, webpage_dir, mapping_dict, coord_units_in_microns, slide_id_col='Slide ID', roi_id_col='Slide ID', x_coord_col='Cell X Position', y_coord_col='Cell Y Position', plot_real_data=True, roi_figsize=roi_figsize, marker_size_step=0.80, alpha=0.1, edgecolors=None, default_marker_size_fac=1, roi_dpi=200, yaxis_dir=yaxis_dir, boxes_to_plot0=None, webpage_dir_addon='whole_slide_patches')

            # Get the unique slides in the dataframe
            unique_slides = df['Slide ID'].unique()

            # Initialize the list that will hold one dataframe per slide whose rows will contain the coordinates of all the ROIs in that slide
            rectangle_data = []

            # For each slide...
            for unique_slide in unique_slides:

                # Get the data for the current slide
                curr_df = df[df['Slide ID'] == unique_slide]

                # Get the ROIs and the coordinates of the objects grouped by ROI
                grouped_rois = curr_df[['tag', 'Cell X Position', 'Cell Y Position']].groupby(by='tag')

                # Save the smallest coordinates of each ROI since this will be used twice below
                roi_min_coords = grouped_rois.min()

                # Save the width and the height of every ROI into a new dataframe called curr_rectangle_data
                curr_rectangle_data = (grouped_rois.max() - roi_min_coords).rename({'Cell X Position': 'width', 'Cell Y Position': 'height'}, axis='columns')

                # Add a column containing a tuple of the minimum x- and y-coordinates, just as patches.Rectangle() requires
                curr_rectangle_data['xy'] = roi_min_coords.apply(lambda x: tuple(x), axis='columns')

                # Save this ROI coordinate dataframe to the master list
                rectangle_data.append(curr_rectangle_data)

            # Again plot the whole slides, but this time overplot the rectangles of all the ROIs
            plot_just_rois(depatched, plotting_map, num_colors, webpage_dir, mapping_dict, coord_units_in_microns, slide_id_col='Slide ID', roi_id_col='Slide ID', x_coord_col='Cell X Position', y_coord_col='Cell Y Position', plot_real_data=True, roi_figsize=roi_figsize, marker_size_step=0.80, alpha=0.1, edgecolors=None, default_marker_size_fac=1, roi_dpi=200, yaxis_dir=yaxis_dir, boxes_to_plot0=rectangle_data, filename_suffix='-patched', webpage_dir_addon='whole_slide_patches')

        # Print that we're not running this component
        else:
            print('Directory {} already exists; not plotting whole slide patches now'.format(os.path.join(webpage_dir, 'whole_slide_patches')))



    def plot_outline_for_single_roi_on_whole_slide(self, roi_figsize=(10, 6.67), depatch=True, nworkers=1, use_multiprocessing=True, yaxis_dir=-1):
        """Plot the outline of a single ROI on top of the slides.

        This is based on the plot_whole_slide_patches() method above.

        Args:
            roi_figsize (tuple, optional): Size (width, height) in inches of the output figures. Defaults to (15, 10).
            depatch (bool, optional): Whether to perform depatching first. Defaults to True.
        """

        # Import relevant library
        import os

        # Define variables already defined as attributes
        df = self.data  # get a shortcut for the main data dataframe
        plotting_map = self.plotting_map
        num_colors = self.num_colors
        mapping_dict = self.mapping_dict
        coord_units_in_microns = self.coord_units_in_microns
        webpage_dir = self.webpage_dir

        # If the directory holding the ROI plots does not exist, create and populate it
        if not os.path.exists(os.path.join(webpage_dir, 'single_roi_outlines_on_whole_slides')):

            # Get rid of any overlap created by the patches so that the same cells aren't plotted multiple times, showing a gridded pattern
            if depatch:
                depatched, _ = undo_patching_overlaps_and_decompounding(df)
            else:
                depatched = df

            # Get the unique slides in the dataframe
            unique_slides = df['Slide ID'].unique()

            # Initialize the list that will hold one dataframe per slide whose rows will contain the coordinates of all the ROIs in that slide
            rectangle_data = []

            # For each slide...
            for unique_slide in unique_slides:

                # Get the data for the current slide
                curr_df = df[df['Slide ID'] == unique_slide]

                # Get the ROIs and the coordinates of the objects grouped by ROI
                grouped_rois = curr_df[['tag', 'Cell X Position', 'Cell Y Position']].groupby(by='tag')

                # Save the smallest coordinates of each ROI since this will be used twice below
                roi_min_coords = grouped_rois.min()

                # Save the width and the height of every ROI into a new dataframe called curr_rectangle_data
                curr_rectangle_data = (grouped_rois.max() - roi_min_coords).rename({'Cell X Position': 'width', 'Cell Y Position': 'height'}, axis='columns')

                # Add a column containing a tuple of the minimum x- and y-coordinates, just as patches.Rectangle() requires
                curr_rectangle_data['xy'] = roi_min_coords.apply(lambda x: tuple(x), axis='columns')

                # Save this ROI coordinate dataframe to the master list
                rectangle_data.append(curr_rectangle_data)

            # Plot the whole slides, but this time overplot the rectangles of the ROIs
            # Note that roi_id_col='Slide ID' is set that way because we want to treat the entire slide as a single ROI
            plot_just_rois(depatched, plotting_map, num_colors, webpage_dir, mapping_dict, coord_units_in_microns, slide_id_col='Slide ID', roi_id_col='Slide ID', x_coord_col='Cell X Position', y_coord_col='Cell Y Position', plot_real_data=True, roi_figsize=roi_figsize, marker_size_step=0.80, alpha=0.1, edgecolors=None, default_marker_size_fac=1, roi_dpi=200, yaxis_dir=yaxis_dir, boxes_to_plot0=rectangle_data, plot_just_single_roi_outline=True, webpage_dir_addon='single_roi_outlines_on_whole_slides', nworkers=nworkers, use_multiprocessing=use_multiprocessing)

        # Print that we're not running this component
        else:
            print('Directory {} already exists; not plotting outlines for each ROI on the whole slides now'.format(os.path.join(webpage_dir, 'single_roi_outlines_on_whole_slides')))


    def average_over_rois_per_annotation_region(self, annotations_csv_files, phenotyping_method, phenotype_identification_file, weight_rois_by_annotation=True, marker_column_names_str='[\'NOS2\', \'COX2\', \'CD8\']', marker_column_names_list=[], annotation_coord_units_in_microns=0.325, alpha=0.4, axis_buffer_frac=0.05, figsize=(10.0, 10.0), annotation_microns_per_integer_unit=0.325, settings__analysis__thickness=40, save_figures=True, also_save_pixel_plot=True, check_averaging_method=True, plot_real_data=True, settings__plotting__pval_figsize=(10, 4), settings__plotting__pval_dpi=100, min_log_pval_for_plotting=-50):

        # Import relevant libraries
        import numpy as np
        import pandas as pd
        import sys
        import annotations
        import os

        # Define variables already defined as attributes
        df_data_by_roi = self.df_data_by_roi
        df_density_pvals_arrays = self.df_density_pvals_arrays
        webpage_dir = self.webpage_dir
        pickle_dir = self.pickle_dir
        all_species_names = self.all_species_names
        log_pval_range = self.log_pval_range

        # Constants
        weight_types_cols = ['slide', 'annotation', 'weight_column_prefix', 'do_log_transform']
        weighted_average_col = ['log_dens_pvals_arr']
        zero = 1e-8
        entity = 'annotation'

        # New pickle filename for the checked and put-into-array-form metrics
        pickle_file = 'density_pvals_averaged_per_region.pkl'

        # If the pickle file doesn't already exist...
        if not os.path.exists(os.path.join(pickle_dir, pickle_file)):

            # Add raw weights values to df_data_by_roi and execute (and reformat) weighted averaging using a "new" calculation method, returning both the final and intermediate values in df_data_holder_new
            df_data_holder_new, df_data_by_roi = annotations.average_over_rois_per_annotation_region_for_all_slides_and_annotations(df_data_by_roi, df_density_pvals_arrays, annotations_csv_files, phenotyping_method, phenotype_identification_file, annotation_coord_units_in_microns=annotation_coord_units_in_microns, alpha=alpha, axis_buffer_frac=axis_buffer_frac, figsize=figsize, annotation_microns_per_integer_unit=annotation_microns_per_integer_unit, settings__analysis__thickness=settings__analysis__thickness, save_figures=save_figures, also_save_pixel_plot=also_save_pixel_plot, equal_weighting_for_all_rois=(not weight_rois_by_annotation), webpage_dir=webpage_dir, log_pval_range=log_pval_range)
            df_new = df_data_holder_new.rename({'weighted_average': weighted_average_col[0]}, axis='columns')[weight_types_cols + weighted_average_col].sort_values(weight_types_cols).reset_index(drop=True)

            # If we want to compare the above (newer) weighting method to one of the original weighting methods as a sanity check...
            if check_averaging_method:

                # Combine the two main dataframes into one
                df_joined = df_data_by_roi.set_index('unique_roi').join(df_density_pvals_arrays.set_index('roi_name'))

                # Determine the number of possible species in the dataset
                # nspecies = df_joined.iloc[0]['log_dens_pvals_arr'].shape[0]
                srs_tmp = df_joined['log_dens_pvals_arr']
                nspecies = srs_tmp[srs_tmp.apply(lambda x: type(x) == np.ndarray)].iloc[0].shape[0]

                # Structures holding iterables of all possible parameters
                annotation_types = [x.removeprefix('num_ann_objects_within_roi_') for x in df_joined.columns if x.startswith('num_ann_objects_within_roi_')]
                weight_column_prefixes = ['num_ann_objects_within_roi_', 'footprint_integer_area_']
                do_log_transforms = [False, True]
                unique_slides = df_joined['unique_slide'].unique()

                # Iterate over all possible parameters
                data_holder_old = []
                for selected_annotation in annotation_types:
                    for weight_column_prefix in weight_column_prefixes:
                        for do_log_transform in do_log_transforms:
                            current_weights_for_all_slides = []  # create a holder of the weights for all slides
                            for unique_slide in unique_slides:

                                # Obtain a fresh copy of the original joined dataframe for the current slide
                                df_joined_slide = df_joined[df_joined['unique_slide'] == unique_slide].copy()

                                # Determine the column name in df_data_by_roi to use as the weights, which were calculated from the annotation data
                                column_for_roi_colors = '{}{}'.format(weight_column_prefix, selected_annotation)

                                # Depending on whether the users selects to use log transformed data, set the transformation function
                                if do_log_transform:
                                    transformation_func = np.log10
                                else:
                                    transformation_func = lambda x: x

                                # For the selected slide of data and for the selected annotation and column to use for the weights, apply the transform to the weights (either the identity or log)
                                weights = df_joined_slide[column_for_roi_colors].apply(transformation_func)

                                # If there's at least 1 ROI with a valid weight...
                                if weights.notna().sum() >= 1:

                                    # Force the weights into the range [0, 1]
                                    weights, _, _ = annotations.standardize_weights_range(weights)

                                    # Rename the weights column to "weights"
                                    weights = weights.rename('weights').apply(lambda weight: np.ones((nspecies, nspecies, 1)) * weight)

                                # If there are no ROIs with a valid (scalar) weight, set the weights for all the ROIs to np.nan
                                else:
                                    print('NOTE: Parameter combination {}-{}-{}-{} has no valid weights; setting weights to np.nan'.format(unique_slide, selected_annotation, weight_column_prefix, do_log_transform))
                                    weights.loc[:] = np.nan

                                # Add the current slide's weights to the weights holder for all the slides
                                current_weights_for_all_slides.append(weights)

                            # Get a Pandas Series of the weights for all slides for the current annotation type, weight column prefix, and do log transform
                            srs_weights_all_slides = pd.concat(current_weights_for_all_slides)

                            # Change the index of df_density_pvals_arrays (1081 rows) from the integer index (but save it) to the ROI name
                            df_density_pvals_arrays_roi_index = df_density_pvals_arrays.reset_index().set_index('roi_name')

                            # From the series srs_weights_all_slides holding the new weights for all the slides (1299 rows), pick just the ROIs/rows present in df_density_pvals_arrays
                            srs_weights_all_slides_reduced_set = srs_weights_all_slides.loc[df_density_pvals_arrays_roi_index.index]

                            # Check that the indexes are the same
                            print('All indexes are the same:', df_density_pvals_arrays_roi_index.index.equals(srs_weights_all_slides_reduced_set.index))

                            # Change the num_valid_centers column to a "fake" one containing the new weights in a (3, 3, 1) np.ndarray
                            df_density_pvals_arrays_roi_index['num_valid_centers'] = srs_weights_all_slides_reduced_set

                            # Reset the index of this dataframe back to the original integer index so that it's identical to the original dataframe df_density_pvals_arrays aside from the num_valid_centers column
                            df_density_pvals_arrays_slide = df_density_pvals_arrays_roi_index.reset_index().set_index('roi_index')

                            # Run a function for weighted averaging based on the original averaging over all ROIs in a slide and obtain the final weighted averages and intermediate results
                            final_and_intermediate_vals_per_slide = average_dens_pvals_over_rois_for_each_slide_and_annotation_type(df_data_by_roi, df_density_pvals_arrays_slide, weight_rois_by_num_valid_centers=weight_rois_by_annotation)

                            # Separate the results back out into a per-slide basis and add the data for each slide to the data holder for the old way of calculating the averages
                            for unique_slide in unique_slides:
                                srs_curr_slide = final_and_intermediate_vals_per_slide.loc[unique_slide, :]
                                data_holder_old.append({'slide': unique_slide, 'annotation': selected_annotation, 'weight_column_prefix': weight_column_prefix, 'do_log_transform': do_log_transform, 'log_dens_pvals_arr': srs_curr_slide['log_dens_pvals_arr'], 'nan_locations': srs_curr_slide['nan_locations'], 'log_dens_pvals_arr_intermediate': srs_curr_slide['log_dens_pvals_arr_intermediate'], 'num_valid_centers': srs_curr_slide['num_valid_centers'], 'weights': srs_curr_slide['weights']})

                # Convert (and reformat) the averages and intermediate data calculated the old way into a Pandas dataframe
                df_data_holder_old = pd.DataFrame(data_holder_old)
                df_old = df_data_holder_old[weight_types_cols + weighted_average_col].sort_values(weight_types_cols).reset_index(drop=True)

                # For every weight type...
                for irow in range(len(df_new)):

                    # Store the new and old weighted averages
                    log_dens_pvals_arr_new = df_new.iloc[irow, -1]
                    log_dens_pvals_arr_old = df_old.iloc[irow, -1]

                    # If either of these is None, ensure they're *both* None
                    current_weighted_averages_are_equal = False
                    if (log_dens_pvals_arr_new is None) or (log_dens_pvals_arr_old is None):
                        current_weighted_averages_are_equal = (log_dens_pvals_arr_new is None) and (log_dens_pvals_arr_old is None)

                    # If they're both np.ndarrays, ensure all elements are basically equal
                    else:
                        arr_is_equal = True
                        for curr_old, curr_new in zip(log_dens_pvals_arr_old.flatten(), log_dens_pvals_arr_new.flatten()):
                            if np.isnan(curr_old) or np.isnan(curr_new):
                                if not (np.isnan(curr_old) and np.isnan(curr_new)):
                                    arr_is_equal = False
                                    break
                            else:
                                if np.abs(curr_new - curr_old) > zero:
                                    arr_is_equal = False
                                    break
                        current_weighted_averages_are_equal = arr_is_equal and np.allclose(log_dens_pvals_arr_old, log_dens_pvals_arr_new, equal_nan=True)  # test np.allclose() simultaneously now and just use that in the future if all is good for a while

                    # If the current old and new weighted averages are not equal, this is bad, so quit the program
                    if not current_weighted_averages_are_equal:
                        print('ERROR: The old and new methods for calculating the weighted averages of the logs of the P values do not agree')
                        sys.exit()

                # If we've gotten this far, then the weighted averages for all weight types are equal, so note this and return one of the results
                print('Great! The weighted averages of the logs of the density P values calculated using very different methods agree with each other!')

            # Create a pickle file containing the results
            make_pickle((df_new, df_data_by_roi), pickle_dir, pickle_file)

            # Define the directory holding all the images of the averaged data
            plots_dir = os.path.join(webpage_dir, 'dens_pvals_per_{}'.format(entity))
            if not os.path.exists(plots_dir):
                os.makedirs(plots_dir)

            # Plot the averaged P values for every slide and annotation/weight type
            for _, row_data in df_new.iterrows():
                slide_name = row_data['slide']
                selected_annotation = row_data['annotation']
                weight_column_prefix = row_data['weight_column_prefix']
                do_log_transform = row_data['do_log_transform']
                log_dens_pvals_arr = row_data[weighted_average_col[0]]
                entity_name = '{}-{}-{}-{}'.format(slide_name, selected_annotation, weight_column_prefix, do_log_transform)
                if log_dens_pvals_arr is not None:
                    print('Plotting density P values for weights {} averaged over all its ROIs with valid P value data'.format(entity_name))
                    plot_density_pvals_simple(
                        log_dens_pvals_arr=log_dens_pvals_arr,
                        # log_pval_range=log_pval_range,
                        log_pval_range=(min_log_pval_for_plotting, 0),
                        figsize=settings__plotting__pval_figsize,
                        dpi=settings__plotting__pval_dpi,
                        plots_dir=plots_dir,
                        plot_real_data=plot_real_data,
                        entity_name=entity_name,
                        img_file_suffix='',
                        entity=entity,
                        entity_index=-1,
                        all_species_names=all_species_names,
                        title_suffix=''
                    )
                else:
                    print('NOT plotting P values for annotation {} because there is no such data'.format(entity_name))

            # # Print that we're not running the plotting part of this component
            # else:
            #     print('Directory {} already exists; not plotting annotation-based density heatmaps now'.format(webpage_dir2))

        # If the pickle file already exists, load the data from it
        else:
            df_new, df_data_by_roi = load_pickle(pickle_dir, pickle_file)

        # Save the calculated data as properties of the class object
        self.df_pvals_averaged_per_annotation = df_new  # weighted averages (and intermediate results)
        self.df_data_by_roi = df_data_by_roi  # updated df_data_by_roi dataframe containing the raw weights values


def average_dens_pvals_over_rois_for_each_slide_and_annotation_type(df_data_by_roi, df_density_pvals_arrays, weight_rois_by_num_valid_centers=True):
    """Average the P values over all ROIs for each slide. Note I have confirmed that this yields the same results as a non-weighted version of the original method.

    Args:
        figsize (tuple, optional): P value figure size in inches. Defaults to (14, 4).
        dpi (int, optional): DPI of the figure. Defaults to 150.
        img_file_suffix (str, optional): String you want to append to the image names. Defaults to ''.
        plot_real_data (bool, optional): Whether to plot the real or simulated data. This is likely meaningless now but was useful historically. Defaults to True.
        do_plotting (bool, optional): Whether to save the average P value plots in addition to calculating the average P values. Defaults to True.
    """

    # Import relevant libraries
    import numpy as np
    import pandas as pd
    import numpy.ma as ma

    # Determine the slides in the full dataset
    unique_slides = df_data_by_roi['unique_slide'].unique()

    # Determine some properties of the final array from a sample log_dens_pvals_arr array
    sample_log_dens_pvals_arr_idx = df_density_pvals_arrays[df_density_pvals_arrays['log_dens_pvals_arr'].apply(lambda x: type(x) == np.ndarray)]['log_dens_pvals_arr'].index[0]
    sample_num_valid_centers_idx = df_density_pvals_arrays[df_density_pvals_arrays['num_valid_centers'].apply(lambda x: type(x) == np.ndarray)]['num_valid_centers'].index[0]
    sample_log_dens_pvals_arr = df_density_pvals_arrays.loc[sample_log_dens_pvals_arr_idx, 'log_dens_pvals_arr']
    sample_num_valid_centers = df_density_pvals_arrays.loc[sample_num_valid_centers_idx, 'num_valid_centers']
    constant_shape_log_dens_pvals_arr = sample_log_dens_pvals_arr.shape
    constant_shape_num_valid_centers = sample_num_valid_centers.shape
    pvals_dtype = sample_log_dens_pvals_arr.dtype
    num_centers_dtype = sample_num_valid_centers.dtype

    # Print what we're doing
    print('  * Weighting the ROIs by the number of valid centers: {}'.format(weight_rois_by_num_valid_centers))

    # Get a list of the indexes of all ROIs having filtered (e.g., minimum number of valid centers) P value data
    roi_indexes_with_pval_data = df_density_pvals_arrays.index.to_list()

    # Initialize the arrays that will hold the numpy array averages for each slide (plus all the intermediate value holders)
    log_dens_pvals_arr_per_slide, nan_locations_per_slide, log_dens_pvals_arr_per_slide2, num_valid_centers_per_slide, weights_per_slide = [], [], [], [], []

    # For every slide in the dataset, regardless of whether any valid P value data actually exists...
    for slide_name in unique_slides:

        # Determine the indexes of all possible ROIs in the slide
        roi_indexes_in_slide = df_data_by_roi[df_data_by_roi['unique_slide'] == slide_name].index.to_list()

        # Filter these out based on whether actual valid P value data exist and determine how many such valid ROIs there are
        roi_indexes_found = [x for x in roi_indexes_in_slide if x in roi_indexes_with_pval_data]  # could have used set intersection
        num_rois_found = len(roi_indexes_found)

        # Print out how many ROIs have been found for the current slide
        print('{} ROIs (indexes {}) have been found for slide {}'.format(num_rois_found, roi_indexes_found, slide_name))

        # print('num_valid_centers', df_density_pvals_arrays.loc[roi_indexes_found, 'num_valid_centers'])
        # print('all weights are nan:', df_density_pvals_arrays.loc[roi_indexes_found, 'num_valid_centers'].apply(lambda x: type(x) == float).sum() == len(roi_indexes_found))
        all_weights_are_nan = df_density_pvals_arrays.loc[roi_indexes_found, 'num_valid_centers'].apply(lambda x: type(x) == float).sum() == len(roi_indexes_found)

        # If there's at least one ROI for the current slide and there's at least one real weight...
        if (num_rois_found >= 1) and (not all_weights_are_nan):

            # Initialize the array holding the data for all the ROIs for the current slide
            log_dens_pvals_arr_all_rois_in_slide = np.zeros(constant_shape_log_dens_pvals_arr + (num_rois_found,), dtype=pvals_dtype)
            num_valid_centers_all_rois_in_slide = np.zeros(constant_shape_num_valid_centers + (num_rois_found,), dtype=num_centers_dtype)

            # For every valid ROI for the current slide, add that ROI's P value data to the overall array
            for iroi, (log_dens_pvals_arr, num_valid_centers) in enumerate(zip(df_density_pvals_arrays.loc[roi_indexes_found, 'log_dens_pvals_arr'], df_density_pvals_arrays.loc[roi_indexes_found, 'num_valid_centers'])):
                log_dens_pvals_arr_all_rois_in_slide.reshape((-1, num_rois_found))[:, iroi] = log_dens_pvals_arr.flatten()
                num_valid_centers_all_rois_in_slide.reshape((-1, num_rois_found))[:, iroi] = num_valid_centers.flatten()

            # Average this array holding all the ROI data over all the ROIs
            if weight_rois_by_num_valid_centers:
                num_valid_centers_all_rois_in_slide = np.tile(num_valid_centers_all_rois_in_slide[:, :, np.newaxis, :, :], (1, 1, 2, 1, 1))
            else:
                num_valid_centers_all_rois_in_slide = np.ones(log_dens_pvals_arr_all_rois_in_slide.shape)

            # Perform weighted averaging using the number of valid centers as weights while taking into account np.nan values in the log of the density P values
            nan_locations = np.isnan(log_dens_pvals_arr_all_rois_in_slide)
            log_dens_pvals_arr_all_rois_in_slide = ma.array(log_dens_pvals_arr_all_rois_in_slide, mask=nan_locations)
            num_valid_centers_all_rois_in_slide = ma.array(num_valid_centers_all_rois_in_slide, mask=nan_locations)
            weights = num_valid_centers_all_rois_in_slide / num_valid_centers_all_rois_in_slide.sum(axis=-1, keepdims=True)  # .data add this to the end if things don't work
            # log_dens_pvals_arr_per_slide.append(np.nanmean(log_dens_pvals_arr_all_rois_in_slide, axis=-1))  # good way to do it if not using weighting and masked array module

            # Need the following three lines (this involving "ans") or else all-nan slices seem to get filled in by zeros
            ans = (log_dens_pvals_arr_all_rois_in_slide * weights).sum(axis=-1)
            ans.set_fill_value(np.nan)

            # Store the final and intermediate values
            log_dens_pvals_arr_per_slide.append(ans.filled())
            nan_locations_per_slide.append(nan_locations)
            log_dens_pvals_arr_per_slide2.append(log_dens_pvals_arr_all_rois_in_slide)
            num_valid_centers_per_slide.append(num_valid_centers_all_rois_in_slide)
            weights_per_slide.append(weights)

        # If there is no data for the current slide, just return None
        else:
            log_dens_pvals_arr_per_slide.append(None)
            nan_locations_per_slide.append(None)
            log_dens_pvals_arr_per_slide2.append(None)
            num_valid_centers_per_slide.append(None)
            weights_per_slide.append(None)

    # Convert the slide P values to a Pandas dataframe and save this as a property of the object
    final_and_intermediate_vals_per_slide = pd.DataFrame(data={'log_dens_pvals_arr': log_dens_pvals_arr_per_slide, 'nan_locations': nan_locations_per_slide, 'log_dens_pvals_arr_intermediate': log_dens_pvals_arr_per_slide2, 'num_valid_centers': num_valid_centers_per_slide, 'weights': weights_per_slide}, index=unique_slides)

    # Return all the final and intermediate data
    return final_and_intermediate_vals_per_slide


def incorporate_species_equivalents(species_equivalents, plotting_map, data, mapping_dict=None):
    '''
    Redefine some species as other species since the thresholding that Houssein sets can sometimes be ambiguous, e.g., large CD163 membranes sometimes result in cells positive for CD163 even though they're not.

    Call like: plotting_map, data, num_colors, unique_species = incorporate_species_equivalents(species_equivalents, plotting_map, data)
    '''

    # Import relevant libraries
    import pandas as pd
    import numpy as np

    # # Define some arrays from already-defined attributes
    # plotting_map = self.plotting_map
    # data = self.data

    # Get a dataframe version of the plotting map
    df_plotting_map = pd.DataFrame(data=plotting_map, columns=['species_id', 'positive_markers', 'species_count', 'marker_ids', 'circle_sizes']).set_index('species_id', drop=False)

    # Sort by species prevalence
    df_plotting_map = df_plotting_map.sort_values('species_count', ascending=False)

    # Ensure the attributes of the object are overwritten
    plotting_map = df_plotting_map.to_numpy()
    num_colors = (df_plotting_map['circle_sizes'] == 1).sum()
    unique_species = df_plotting_map['species_id'].to_numpy(dtype=np.int64)

    return(plotting_map, data, num_colors, unique_species)


# Get the average spacing on either side of each datapoint in an array
def get_avg_spacing(arr):
    if len(arr) >= 2:
        import numpy as np
        arr2 = np.concatenate(([2 * arr[0] - arr[1]], arr, [2 * arr[-1] - arr[-2]]))
        return((arr2[2:] - arr2[0:-2]) / 2)
    else:
        print('Not actually getting average spacing in arr because len(arr) < 2; returning 1')
        return([1])


# Read in the Consolidated_data.txt TSV file into a Pandas dataframe
def get_consolidated_data(csv_file):
    import pandas as pd
    return(pd.read_csv(csv_file, sep='\t'))  # read in the data


# Given a list of phenotypes in a species, return the A+/B+ etc. string version
def phenotypes_to_string(phenotype_list):
    return '/'.join(sorted(phenotype_list))


# Given a list of phenotypes in a species, return the nicely formatted version, if there's a known cell type corresponding to the species
# Note these are purely labels; the species themselves are determined by allow_compound_species as usual
def get_descriptive_cell_label(phenotype_list, mapping_dict):
    # Note: CD163+/CD4+ REALLY ISN'T ANYTHING COMPOUND --> MAKE IT OVERLAPPING SPECIES (i.e., it shouldn't be in the dictionary below)!!!!
    return phenotypes_to_string(phenotype_list), False


# Obtain the plotting map, total number of unique colors needed for plotting, the list of unique species (in the same order as in plotting_map), and a correctly sorted list of slides (e.g., 1,2,15 instead of 1,15,2)
# Note that individual unique species are specified by the allow_compound_species keyword, which in turn affects which of the 'Species int' columns of the Pandas dataframe are actually unique
# Don't get confused by the mapping_dict variable, which only affects plotting of the species... it doesn't affect what is actually considered a unique species or not!
def get_dataframe_info(data, phenotypes, species_int_to_pheno_name=None):

    # Import relevant modules
    import numpy as np
    from operator import itemgetter
    import pandas as pd

    # Create an ndarray containing all the unique species in the dataset in descending order of frequency with columns: integer label, string list, frequency, color(s), circle size(s)
    plotting_map = [[-(list(data['Species int']).count(x)), list(int2list(phenotypes, x)), x] for x in np.unique(data['Species int'])]  # create a list of the unique species in the dataset with columns: -frequency, string list, integer label
    plotting_map.sort(key=itemgetter(0, 1))  # sort by decreasing frequency (since the frequency column is negative) and then by the string list
    plotting_map = [[-x[0], x[1], x[2]] for x in plotting_map]  # make the frequency column positive

    print(plotting_map)

    # Get the colors of the species that are already known to us; use a -1 if the species isn't known
    colors = []
    known_phenotypes = []
    known_colors = []
    icolor = 0
    for item in plotting_map:
        phenotype_list = item[1]
        # If the species (each row of plotting_map) is known to us (i.e., in the inputted mapping_dict variable, which simply assigns a cell label to any single or compound species)...
        # ...give that species its own color, and make a note if the species is also a single, non-compound species (i.e., a single phenotype)
        colors.append(-1)

    # Get the colors of the rest of the species using the colors of the already-known single-phenotype species
    # I.e., if the species is not known to us (i.e., not in mapping_dict), do not give the species its own color (unless it contains a phenotype that's not in known_phenotypes)
    # Instead, give each phenotype in the species either the corresponding color in known_phenotypes (if it's in there) or a new color (if it's not in known_phenotypes)
    # Assign the corresponding circle sizes as well
    colors2 = []
    circle_sizes = []
    for item, color in zip(plotting_map, colors):
        phenotype_list = item[1]
        if color == -1:
            curr_colors = []
            for single_phenotype in phenotype_list:
                if single_phenotype in known_phenotypes:
                    curr_colors.append(known_colors[known_phenotypes.index(single_phenotype)])
                else:
                    curr_colors.append(icolor)
                    known_phenotypes.append(single_phenotype)
                    known_colors.append(icolor)
                    icolor = icolor + 1
        else:
            curr_colors = [color]

        # Always have the most prevalent single species (if a lower color number really implies higher prevalence, it should generally at least) first in the color list, and make the corresponding circle size the largest (but in the background of course)
        curr_colors.sort()
        curr_sizes = list(np.arange(start=len(curr_colors), stop=0, step=-1))

        colors2.append(curr_colors)
        circle_sizes.append(curr_sizes)
    colors = colors2

    # Store the total number of unique colors to plot
    num_colors = icolor

    # Finalize the plotting map
    plotting_map = np.array([[item[2], item[1], item[0], (color if (len(color) != 1) else color[0]), (circle_size if (len(circle_size) != 1) else circle_size[0])] for item, color, circle_size in zip(plotting_map, colors, circle_sizes)], dtype=object)

    # Use the plotting map to extract just the unique species in the data
    unique_species = np.array([x[0] for x in plotting_map])  # get a list of all the unique species in the dataset in the correct order

    # Get the unique slides sorted correctly
    tmp = [[int(x.split('-')[0][0:len(x.split('-')[0]) - 1]), x] for x in np.unique(data['Slide ID'])]
    tmp.sort(key=(lambda x: x[0]))
    unique_slides = [x[1] for x in tmp]

    # Get a dataframe version of the plotting map
    df_plotting_map = pd.DataFrame(data=plotting_map, columns=['species_id', 'positive_markers', 'species_count', 'marker_ids', 'circle_sizes']).set_index('species_id', drop=False)

    # Sort by species prevalence
    df_plotting_map = df_plotting_map.sort_values('species_count', ascending=False)

    # Overwrite the positive marker list with the assigned phenotype name from new_phenotyping_lib.apply_phenotyping()
    if species_int_to_pheno_name is not None:
        df_plotting_map['positive_markers'] = df_plotting_map['species_id'].apply(lambda x: [species_int_to_pheno_name[x]])
        df_plotting_map['marker_ids'] = range(len(df_plotting_map))
        df_plotting_map['circle_sizes'] = 1

    # Ensure the attributes of the object are overwritten
    plotting_map = df_plotting_map.to_numpy()
    num_colors = (df_plotting_map['circle_sizes'] == 1).sum()
    unique_species = df_plotting_map['species_id'].to_numpy(dtype=np.int64)

    return(plotting_map, num_colors, unique_species, unique_slides, data)


# Given an array of densities, for each of which we will generate a different ROI of the corresponding density, return a Pandas dataframe of the simulated data
# Only a single slide will be returned (but in general with multiple ROIs)
def get_simulated_data(doubling_type, densities, max_real_area, min_coord_spacing, mult):

    # Import relevant modules
    import numpy as np
    import pandas as pd

    # Warn if we're doubling up the data
    if doubling_type != 0:
        print('NOTE: DOUBLING UP THE SIMULATED DATA!!!!')

    # Specify the columns that are needed based on what's used from consolidated_data.txt
    columns = ['tag', 'Cell X Position', 'Cell Y Position', 'Slide ID', 'Phenotype A', 'Phenotype B']

    # Determine a constant number of cells in a simulated ROI using the maximum-sized ROI of the real data and the largest simulated density
    N = int(np.int(densities[-1] * max_real_area) * mult)

    # Initialize the Pandas dataframe
    data = pd.DataFrame(columns=columns)

    # For each ROI density and average either-side density spacing...
    perc_error = []
    for pop_dens, avg_spacing in zip(densities, get_avg_spacing(densities)):
        print('Current population density:', pop_dens)

        # Get the Cartesian coordinates of all the cells from random values populated until the desired density is reached
        tot_area = N / pop_dens
        side_length = np.sqrt(tot_area)
        tmp = np.round(side_length / min_coord_spacing)  # now we want N random integers in [0,tmp]
        coords_A = np.random.randint(tmp, size=(int(2 * N / 3), 2)) * min_coord_spacing
        coords_B = np.random.randint(tmp, size=(int(1 * N / 3), 2)) * min_coord_spacing
        x_A = coords_A[:, 0]
        y_A = coords_A[:, 1]
        x_B = coords_B[:, 0]
        y_B = coords_B[:, 1]

        # Set the ROI name from the current density
        tag = 'pop_dens_{:09.7f}'.format(pop_dens)  # IN GENERAL THIS SHOULD PROBABLY INCLUDE THE SLIDE ID JUST LIKE THE CONSOLIDATED_DATA FILES DO! (i.e., something like '1A-only_slide_pop_dens_{:09.7f}'). This is because in general the ROI name should be slide-specific.

        # Add the simulated data to a nested list and then convert to a Pandas dataframe and add it to the master dataframe
        # Number in bracket is number of actually unique dataset (three actually unique datasets)
        # (1) [1] doubling_type=0, allow_compound_species=True:  coords=(A,B), labels=(A,B) # two species with different coordinates
        # (2) [1] doubling_type=0, allow_compound_species=False: coords=(A,B), labels=(A,B)
        # (3) [2] doubling_type=1, allow_compound_species=True:  coords=(A,A), labels=(A,B) # two overlapping species
        # (4) [2] doubling_type=1, allow_compound_species=False: coords=(A,A), labels=(A,B)
        # (5) [3] doubling_type=2, allow_compound_species=True:  coords=(A),   labels=(AB)  # one compound species (AB = compound species)
        # (6) [2] doubling_type=2, allow_compound_species=False: coords=(A,A), labels=(A,B)
        if doubling_type != 2:
            list_set = [[tag, curr_x_A, curr_y_A, '1A-only_slide', 'A+', 'B-'] for curr_x_A, curr_y_A in zip(x_A, y_A)]
            if doubling_type == 0:
                list_set = list_set + [[tag, curr_x_B, curr_y_B, '1A-only_slide', 'A-', 'B+'] for curr_x_B, curr_y_B in zip(x_B, y_B)]
            elif doubling_type == 1:
                list_set = list_set + [[tag, curr_x_A, curr_y_A, '1A-only_slide', 'A-', 'B+'] for curr_x_A, curr_y_A in zip(x_A, y_A)]
        else:
            list_set = [[tag, curr_x_A, curr_y_A, '1A-only_slide', 'A+', 'B+'] for curr_x_A, curr_y_A in zip(x_A, y_A)]
        tmp = pd.DataFrame(list_set, columns=columns)
        data = data.append(tmp, ignore_index=True)  # may need to fix due to deprecation of .append(); use pd.concat() instead

        # Calculate the percent error in the actual density from the desired density
        if doubling_type == 0:
            x_tmp = np.r_[x_A, x_B]
            y_tmp = np.r_[y_A, y_B]
        else:
            x_tmp = x_A
            y_tmp = y_A
        perc_error.append((N / ((x_tmp.max() - x_tmp.min()) * (y_tmp.max() - y_tmp.min())) - pop_dens) / avg_spacing * 100)

    print('Percent error:', perc_error)
    print('Maximum percent error:', np.max(perc_error))

    return(data)


# Convert integer numbers defined in species bit-wise to a string list based on phenotypes
def int2list(phenotypes, species):
    return(phenotypes[[bool(int(char)) for char in ('{:0' + str(len(phenotypes)) + 'b}').format(species)]])


# Load some data from a pickle file
def load_pickle(pickle_dir, pickle_file):
    import pickle
    import os
    filename = os.path.join(pickle_dir, pickle_file)
    print('Reading pickle file ' + filename + '...')
    with open(filename, 'rb') as f:
        try:
            data_to_load = pickle.load(f)
        except ModuleNotFoundError:
            import pandas as pd
            try:
                data_to_load = pd.read_pickle(f)  # solution per https://stackoverflow.com/questions/75953279/modulenotfounderror-no-module-named-pandas-core-indexes-numeric-using-metaflo
            except EOFError:
                data_to_load = pd.compat.pickle_compat.load(f)  # another solution per https://stackoverflow.com/questions/75953279/modulenotfounderror-no-module-named-pandas-core-indexes-numeric-using-metaflo
    return(data_to_load)


# Write a pickle file from some data
def make_pickle(data_to_save, pickle_dir, pickle_file):
    import pickle
    import os
    filename = os.path.join(pickle_dir, pickle_file)
    print('Creating pickle file ' + filename + '...')
    if not os.path.exists(pickle_dir):
        os.makedirs(pickle_dir)
    with open(filename, 'wb') as f:
        pickle.dump(data_to_save, f)


def plot_roi(fig, spec2plot, species, x, y, plotting_map, colors, x_range, y_range, title, marker_size_step, default_marker_size, dpi, mapping_dict, coord_units_in_microns, filepath=None, do_plot=True, alpha=1, edgecolors='k', yaxis_dir=1, boxes_to_plot=None, pval_params=None, roi_pval_alpha=0.5, num_colors=2**16, title_suffix=None):
    '''
    For the raw data (coordinates) for a given ROI, plot a circle (scatter plot) representing each species, whether known (in which case get_descriptive_cell_label() is used) or unknown; plot a legend too

    Save the figure to a file as well

    Note that in plotting_map (which is a list made up of five-element lists, one per species in the experiment, in order of decreasing frequency over the entire experiment), each row is what ends up as "plotting_data" below. The fourth and fifth elements (the colors and sizes, respectively) are ordered correctly relative to each other, and (for compound species) are in the order of increasing color, which, due to how plotting_map is defined, means of decreasing frequency. (In general; the first criterion in plotting_map is actually whether the species is "known", i.e., in mapping_dict.) This generally means, for instance, that the more prevalent species is plotting larger and in the back of the less prevalent species in the case of the plotting of a compound species.

    The point is, however, the second element of the plotting_data lists are NOT ordered correctly relative to those two potential-lists, so we take care below to ensure the subspecies are LABELED correctly. After this is done, we do not need the variable icircle, which we are now removing.
    '''

    if do_plot:

        # Import relevant library
        import numpy as np
        import matplotlib.patches as patches
        import seaborn as sns  # we need to do this in order to get the rocket colormap below even though sns is never used

        # Axis
        ax = fig.axes[0]
        ax.cla()

        # Get the lookup table for the choosing of the correct species for labeling below (where label_lookup_for_color is used). This is what we are referring to in the docstring for this function.
        primary_species = [not isinstance(x[3], list) for x in plotting_map]
        plotting_map_primary_noncompound_species = [x for (x, y) in zip(plotting_map, primary_species) if (y and (len(x[1]) == 1))]
        label_lookup_for_color = [[x[1][0], x[3]] for x in plotting_map_primary_noncompound_species]

        # For each unique species in the current ROI (in the correct order)...
        plotted_colors = []
        plots_for_legend = []
        labels_for_legend = []
        for spec in spec2plot:

            # Get the data for that species
            spec_ind = np.nonzero(species == spec)
            x_pos = x[spec_ind]
            y_pos = y[spec_ind]
            plotting_data = plotting_map[[x[0] for x in plotting_map].index(spec)]

            # Ensure the colors and marker sizes are lists and determine whether a single circle will be plotted for the potentially compound species
            if not isinstance(plotting_data[3], list):
                all_colors = [plotting_data[3]]
                all_sizes = [plotting_data[4]]
                is_primary = True
            else:
                all_colors = plotting_data[3]
                all_sizes = plotting_data[4]
                is_primary = False

            # For each circle to plot within the current species...
            for curr_color, curr_size in zip(all_colors, all_sizes):

                # Obtain the actual color and marker size to plot
                curr_color2 = colors[curr_color]
                curr_size2 = (((curr_size - 1) * marker_size_step + 1) * default_marker_size) ** 2

                # Plot the current species
                # Note: coord_units_in_microns is # of microns per unit (same as coord_units_in_microns, the new input parameter to the entire workflow)
                curr_plt = ax.scatter(x_pos * coord_units_in_microns, y_pos * coord_units_in_microns, s=curr_size2, c=curr_color2, edgecolors=edgecolors, alpha=alpha)

                # If we're on a primary species (in which a single circle is plotted for a potentially compound species), add the current plot to the legend
                if is_primary:
                    curr_label = get_descriptive_cell_label(plotting_data[1], mapping_dict)[0]
                    plotted_colors.append(curr_color)  # keep a record of the colors we've plotted so far in order to add a minimal number of non-primary species to the legend
                    plots_for_legend.append(curr_plt)
                    labels_for_legend.append(curr_label)

                # If we're on a non-primary species, only add it to the legend if the color hasn't yet been plotted
                else:

                    # If the color hasn't yet been plotted...
                    if curr_color not in plotted_colors:

                        # Get the correct label for the current phenotype within the non-primary species
                        # curr_label = get_descriptive_cell_label([plotting_data[1][icircle]], mapping_dict)[0] # this assumes the order in plotting_data[1] is consistent with that in plotting_data[3] and plotting_data[4], which is not necessarily true!
                        curr_label = get_descriptive_cell_label([label_lookup_for_color[[x[1] for x in label_lookup_for_color].index(curr_color)][0]], mapping_dict)[0]  # this performs the correct lookup

                        # If the current color to add to the legend was NOT a minimal size, first make a dummy plot of one of the minimal size
                        if not curr_size == 1:
                            curr_plt = ax.scatter(x_range[0] * 2 * coord_units_in_microns, y_range[0] * 2 * coord_units_in_microns, s=(default_marker_size**2), c=curr_color2, edgecolors='k', alpha=alpha)

                        # Add the plot to the legend
                        plotted_colors.append(curr_color)  # keep a record of the colors we've plotted so far in order to add a minimal number of non-primary species to the legend
                        plots_for_legend.append(curr_plt)
                        labels_for_legend.append(curr_label)

        # Complete the plot
        ax.set_aspect('equal')
        ax.set_xlim(tuple(x_range * coord_units_in_microns))
        ax.set_ylim(tuple(y_range * coord_units_in_microns))
        ax.set_xlabel('X coordinate (microns)')
        ax.set_ylabel('Y coordinate (microns)')
        # ax.set_title('ROI ' + title + ('' if title_suffix is None else title_suffix))
        ax.set_title(title + ('' if title_suffix is None else title_suffix))
        ax.legend(plots_for_legend, labels_for_legend, loc='upper left', bbox_to_anchor=(1, 0, 1, 1))
        ax.set_ylim(ax.get_ylim()[::yaxis_dir])

        # Optionally overplot boxes, as in, e.g., the outlines of the ROIs when an entire slide is treated as the ROI referred to by plot_roi()
        # If ROI coordinates are present...
        if boxes_to_plot is not None:

            # For each ROI...
            for row_index, row_data in boxes_to_plot.iterrows():

                # Initialize the facecolor of the ROI to transparent
                facecolor = 'none'

                # If we are inputting into the "ROI" plot a dataframe P value column to plot as well as the P value colorbar extremes...
                if pval_params is not None:

                    # Get the column in boxes_to_plot to plot (e.g., either left_log_dens_pval or right_log_dens_pval), the data minimum and maximum values for setting the colormap extremes, and whether to print the P value for each ROI for debugging purposes
                    pval_colname, vmin, vmax, print_pval_for_each_roi = pval_params

                    # Get the scalar data value to plot
                    log_dens_pval = row_data[pval_colname]

                    # log_dens_pval_is_valid = not np.isnan(log_dens_pval)  # will throw an error if log_dens_pval is ever None, which it should never be, so we're good
                    log_dens_pval_is_valid = not ((log_dens_pval is None) or np.isnan(log_dens_pval))

                    # If this data point exists...
                    if log_dens_pval_is_valid:

                        # Overwrite the facecolor of the ROI to that corresponding to the datapoint, the minimum/maximum values of the colormap extremes, and the colormap name, appending the desired transparency to the final tuple
                        # facecolor = get_rgba_for_data_on_colormap(log_pval=log_dens_pval, vmin=vmin, vmax=vmax, colormap_name='rocket', print_byte_value=False)[:3] + (roi_pval_alpha,)
                        facecolor = get_properly_interpolated_rgba_from_data(val_to_map=log_dens_pval, vmin=vmin, vmax=vmax, colormap_name='rocket', print_byte_value=False, N=num_colors)[:3] + (roi_pval_alpha,)

                        # Print the log density P value of the current ROI for testing purposes
                        if print_pval_for_each_roi:
                            print('ROI {} has a log density P value of {} and is color {}'.format(row_index, log_dens_pval, facecolor))

                # Plot the ROI outlines and the corresponding P value colors if they are valid
                # do this if (1) we're not overplotting pvals (pval_params is None) or (2) if we're overplotting pvals (pval_params is not None), if the data point exists (log_dens_pval_is_valid)
                if (pval_params is None) or ((pval_params is not None) and log_dens_pval_is_valid):
                    ax.add_patch(patches.Rectangle(row_data['xy'], row_data['width'], row_data['height'], linewidth=1, edgecolor='k', facecolor=facecolor))

        # Save the figure
        if filepath is not None:
            fig.savefig(filepath, dpi=dpi, bbox_inches='tight')


# For the inputted ROI coordinates, print the minimum coordinate spacing, coordinate range, number of cells, ROI area, and density
def roi_checks_and_output(x_roi, y_roi, do_printing=True, do_main_printing=True):
    import numpy as np
    # coord_spacing_check = 0.5
    x_roi2 = x_roi.copy()
    y_roi2 = y_roi.copy()
    x_roi2.sort()
    y_roi2.sort()
    unique_spacing_x = np.unique(x_roi2[1:] - x_roi2[0:-1])[0:2]
    unique_spacing_y = np.unique(y_roi2[1:] - y_roi2[0:-1])[0:2]

    if (len(unique_spacing_x) < 2) or (len(unique_spacing_y) < 2):
        min_coordinate_spacing = -7  # even though this block happens, min_coordinate_spacing is never actually used so it shouldln't matter if it's -7
        print('NOTE: The length of the unique spacing array in x or y is less than 2; printing x and y coords:')
        print('  x coords:', x_roi)
        print('  y coords:', y_roi)
    else:
        if unique_spacing_x[1] != unique_spacing_y[1]:
            if do_main_printing:
                print('NOTE: Minimum spacing in x coordinates ({}) is different from that in y coordinates ({})'.format(unique_spacing_x[1], unique_spacing_y[1]))
        min_coordinate_spacing = unique_spacing_x[1]
    if do_main_printing:
        print('Calculated minimum coordinate spacing: {}'.format(min_coordinate_spacing))

    # expected_unique_spacing = [0,coord_spacing_check]
    # if (not (unique_spacing_x==expected_unique_spacing).all()) or (not (unique_spacing_y==expected_unique_spacing).all()): # checks that coord_spacing is coord_spacing_check
    #     print('ERROR: Coordinate spacing is not', coord_spacing_check)
    #     exit()
    x_range = [x_roi2.min(), x_roi2.max()]
    y_range = [y_roi2.min(), y_roi2.max()]
    if do_printing:
        print('x range:', x_range)
        print('y range:', y_range)
    ncells = len(x_roi2)
    area = (x_roi2.max() - x_roi2.min()) * (y_roi2.max() - y_roi2.min())

    if area < 1e-6:
        print('ROI has zero area!')

    if do_printing:
        print('[ncells, area, density]:', [ncells, area, ncells / area])
    # return(area, unique_spacing_x[1])
    return(x_range, y_range, min_coordinate_spacing)


def calculate_metrics_from_coords(min_coord_spacing, input_coords=None, neighbors_eq_centers=False, ncenters_roi=1300, nneighbors_roi=220, nbootstrap_resamplings=0, rad_range=(2.2, 5.1), use_theoretical_counts=False, roi_edge_buffer_mult=1, roi_x_range=(1.0, 100.0), roi_y_range=(0.5, 50.0), silent=False, log_file_data=None, keep_unnecessary_calculations=False, old_method=False):
    '''
    Given a set of coordinates (whether actual coordinates or ones to be simulated), calculate the P values and Z scores.

    See the calculate_metrics() method of the TIMECellInteraction class for further documentation.
    '''

    # Import relevant libraries
    import numpy as np
    import scipy.spatial
    import scipy.stats
    import utils

    # Constant
    tol = 1e-8

    # Create a numpy array version of rad_range
    radii = np.array(rad_range)

    # Determine whether logging is to be done
    if log_file_data is not None:
        do_logging = True
        log_file_handle, roi_index, uroi, center_species, neighbor_species = log_file_data
    else:
        do_logging = False

    # if do_logging:
    #     log_file_handle.write('ROI {:06d} ({}): center {} and neighbor {} both exist in the ROI to any extent\n'.format(roi_index, uroi, center_species, neighbor_species))

    # If input_coords is None, then simulate the coordinates
    if input_coords is not None:
        coords_centers, coords_neighbors = input_coords
        ncenters_roi = coords_centers.shape[0]
        if coords_neighbors is not None:
            nneighbors_roi = coords_neighbors.shape[0]
        else:
            nneighbors_roi = None
        simulate_coords = False
    else:
        simulate_coords = True

    # Calculate some properties of the ROI itself
    roi_area_adj = (roi_x_range[1] - roi_x_range[0] + min_coord_spacing) * (roi_y_range[1] - roi_y_range[0] + min_coord_spacing) - min_coord_spacing**2
    ngridpoints_x = int((roi_x_range[1] - roi_x_range[0]) / min_coord_spacing + 1)
    ngridpoints_y = int((roi_y_range[1] - roi_y_range[0]) / min_coord_spacing + 1)
    grid_indices = np.indices((ngridpoints_x, ngridpoints_y)).reshape((2, -1)).T

    # Properties of the slice
    slice_area = np.pi * (rad_range[1]**2 - rad_range[0]**2)
    slice_area_adj = slice_area - min_coord_spacing**2
    if slice_area_adj < 0:
        print('ERROR: Slice area is too small')
        exit()

    # Calculate the expected number of neighbors in the slice; see physical notebook notes on 1/7/21 for details
    roi_area_used = roi_area_adj
    if (not neighbors_eq_centers) and (rad_range[0] > 0):
        npossible_neighbors = nneighbors_roi
        slice_area_used = slice_area
    elif (not neighbors_eq_centers) and (rad_range[0] < tol):
        npossible_neighbors = nneighbors_roi
        slice_area_used = slice_area_adj
    elif (neighbors_eq_centers) and (rad_range[0] > 0):
        npossible_neighbors = ncenters_roi - 1
        slice_area_used = slice_area
    elif (neighbors_eq_centers) and (rad_range[0] < tol):
        # NOTE THAT IN THIS CASE WE MUST SUBTRACT 1 FROM NNEIGHBORS (if not use_simulated_counts)!!
        npossible_neighbors = ncenters_roi - 1
        slice_area_used = slice_area_adj
    else:
        print('ERROR: Impossible situation encountered')
        exit()
    nexpected = npossible_neighbors / roi_area_used * slice_area_used

    # Initialize the random number generator
    if simulate_coords or (nbootstrap_resamplings != 0):
        rng = np.random.default_rng()

    # Generate random coordinates for the centers and neighbors; by sampling without replacement, we're simulating our real data in that no two species occupy the same point on the grid
    if neighbors_eq_centers:
        if simulate_coords:
            coords_centers = (rng.choice(grid_indices, size=ncenters_roi, replace=False, shuffle=False) * min_coord_spacing) + np.array([roi_x_range[0], roi_y_range[0]])[np.newaxis, :]
        coords_neighbors = coords_centers
    else:
        if simulate_coords:
            coords_tmp = (rng.choice(grid_indices, size=(ncenters_roi + nneighbors_roi), replace=False, shuffle=False) * min_coord_spacing) + np.array([roi_x_range[0], roi_y_range[0]])[np.newaxis, :]
            coords_centers = coords_tmp[:ncenters_roi, :]
            coords_neighbors = coords_tmp[ncenters_roi:, :]

    # Determine and count the centers that are within the larger radius of the ROI boundaries
    if not silent:
        if np.abs(roi_edge_buffer_mult - 1) > tol:
            print('WARNING: roi_edge_buffer_mult should be set to at least the default value of 1 (1 keeps the most data) in order to correctly account for edge effects; it is currently set to {}'.format(roi_edge_buffer_mult))
    valid_centers = (coords_centers[:, 0] <  (roi_x_range[1] + min_coord_spacing / 2 - roi_edge_buffer_mult * rad_range[1])) & \
                    (coords_centers[:, 0] >= (roi_x_range[0] - min_coord_spacing / 2 + roi_edge_buffer_mult * rad_range[1])) & \
                    (coords_centers[:, 1] <  (roi_y_range[1] + min_coord_spacing / 2 - roi_edge_buffer_mult * rad_range[1])) & \
                    (coords_centers[:, 1] >= (roi_y_range[0] - min_coord_spacing / 2 + roi_edge_buffer_mult * rad_range[1]))
    nvalid_centers = valid_centers.sum()

    # As long as there's at least one valid center (i.e., one sample)...
    if nvalid_centers >= 1:

        if do_logging:
            log_file_handle.write('ROI {:06d} (split 03, {}): there are {} valid centers for center {} and neighbor {}\n'.format(roi_index, uroi, nvalid_centers, center_species, neighbor_species))

        # Calculate the number of neighbors between the slice radii around all valid centers
        if use_theoretical_counts:
            if not silent:
                print('NOTE: Using artificial distribution')
            nneighbors = scipy.stats.poisson.rvs(nexpected, size=(nvalid_centers,))
        else:
            if old_method:
                dist_mat = scipy.spatial.distance.cdist(coords_centers[valid_centers, :], coords_neighbors, 'euclidean')  # calculate the distances between the valid centers and all the neighbors
                nneighbors = ((dist_mat >= rad_range[0]) & (dist_mat < rad_range[1])).sum(axis=1)  # count the number of neighbors in the slice around every valid center
            else:
                nneighbors = utils.calculate_neighbor_counts_with_possible_chunking(center_coords=coords_centers[valid_centers, :], neighbor_coords=coords_neighbors, radii=radii, single_dist_mat_cutoff_in_mb=200, verbose=False)[:, 0]  # (num_centers,)
            if (neighbors_eq_centers) and (rad_range[0] < tol):
                nneighbors = nneighbors - 1  # we're always going to count the center as a neighbor of itself in this case, so account for this; see also physical notebook notes on 1/7/21

        # From all possible numbers of neighbors calculate the bin edges
        if keep_unnecessary_calculations:
            nbins = nneighbors.max() + 1
            edges = np.arange(nbins + 1) - 0.5
        else:
            edges = None

        # Redefine the counts of the neighbors around the centers (i.e., the densities), either by adding a singleton dimension or performing bootstrap resampling
        if nbootstrap_resamplings == 0:
            nneighbors = nneighbors[:, np.newaxis]  # first add a new single-dimension axis so that there's only a single "sample"
        else:
            nneighbors = rng.choice(nneighbors, size=(nvalid_centers, nbootstrap_resamplings), replace=True)  # first perform bootstrap resampling

        # Generate histograms (i.e., the non-normalized PMFs) of the numbers of neighbors
        if keep_unnecessary_calculations:
            nsamples = nneighbors.shape[1]
            pmf = np.zeros((nbins, nsamples))
            for isample in range(nsamples):
                pmf[:, isample], _ = np.histogram(nneighbors[:, isample], bins=edges, density=False)  # this quantity (the number of centers having j neighbors) is binomial-distributed and is not technically the PMF since we are not normalizing it here

        # Calculate the Z scores and P values for the densities and PMFs
        density_metrics = calculate_density_metrics(nneighbors, nexpected, keep_unnecessary_calculations=keep_unnecessary_calculations)
        if keep_unnecessary_calculations:
            pmf_metrics = calculate_PMF_metrics(pmf, nexpected, nvalid_centers)
        else:
            pmf_metrics = None

    # If there are no valid centers, don't return anything for the density and PMF metrics
    else:

        if do_logging:
            log_file_handle.write('ROI {:06d} (split 03, {}): there are no valid centers for center {} and neighbor {}. Setting density_metrics, pmf_metrics, and edges to None.\n'.format(roi_index, uroi, center_species, neighbor_species))

        density_metrics = None
        pmf_metrics = None
        edges = None

    return(density_metrics, pmf_metrics, nexpected, nvalid_centers, coords_centers, coords_neighbors, valid_centers, edges, npossible_neighbors, roi_area_used, slice_area_used)


def calculate_density_metrics(nneighbors, nexpected, keep_unnecessary_calculations=False):
    '''
    Calculate the Z score (not generally applicable) and the left and right P values for the "density".
    '''

    # Import relevant libraries
    import numpy as np
    import scipy.stats

    # Calculate the number of samples (and number of valid centers) from the input array, which may be 1 particularly in the case where we've decided not to bootstrap prior to calling this function
    nvalid_centers, nsamples = nneighbors.shape

    # Define the array holding the values of interest: the Z score, the "less" or "left" P value, and the "greater" or "right" P value
    metrics = np.zeros((nsamples, 3))

    # Define the variable that is Poisson-distributed
    k = nneighbors.sum(axis=0)  # (nsamples,)

    # Parameter of the Poisson distribution
    mu = nvalid_centers * nexpected  # scalar

    # Calculate a z-score, which is meaningful only if the distribution is approximately Gaussian
    if keep_unnecessary_calculations:
        mean = mu  # scalar, property of the Poisson distribution
        std = np.sqrt(mu)  # scalar, property of the Poisson distribution
        metrics[:, 0] = (k - mean) / std  # z-score; this is meaningless if the distribution is not approximately Gaussian

    # Calculate the P values
    metrics[:, 1] = scipy.stats.poisson.cdf(k, mu)  # "less" or "left" P value
    metrics[:, 2] = scipy.stats.poisson.sf(k - 1, mu)  # "greater" or "right" P value

    # Average over all the samples to reduce sensitivity in favor of specificity (if bootstrapping)
    return(metrics.mean(axis=0))  # (3,)


def calculate_PMF_metrics(pmf, nexpected, nvalid_centers):
    '''
    Calculate the Z score (not generally applicable) and the left and right P values for the "PMF".

    Note I found for the random points that the return value never had any nans.
    '''

    # Import relevant libraries
    import numpy as np
    import scipy.stats

    # Calculate some array sizes
    nbins, nsamples = pmf.shape

    # Define the metrics array containing the Z score and the left and right P values
    metrics = np.zeros((nbins, nsamples, 3))

    # Define the variable that is binomial-distributed
    k = pmf  # (nbins, nsamples)

    # Parameters of the binomial distribution
    n = nvalid_centers  # scalar
    p = scipy.stats.poisson.pmf(np.arange(nbins), nexpected)[:, np.newaxis]  # (nbins,1)

    # Properties of the binomial distribution
    mean = n * p  # (nbins,1)
    std = np.sqrt(n * p * (1 - p))  # (nbins,1)

    # Calculate the metrics
    metrics[:, :, 0] = (k - mean) / std  # z-score; this is meaningless if the distribution is not approximately Gaussian
    metrics[:, :, 1] = scipy.stats.binom.cdf(k, n, p)  # "less" or "left" P value
    metrics[:, :, 2] = scipy.stats.binom.sf(k - 1, n, p)  # "greater" or "right" P value

    # Average over all the samples to reduce sensitivity in favor of specificity (if bootstrapping)
    return(metrics.mean(axis=1))  # (nbins,3)


def plot_pvals(fig, data_by_slice, log_pval_range, name=None, calculate_empty_bin_pvals=False, max_nbins_over_slices=None, square=True, yticklabels=2):
    '''
    Plot the four-axis P value heatmaps.

    The top row has the left and right density P values and the bottom row has the left and right PMF P values.

    Here we carefully account for when the number of valid centers is 0, which can be true when this function is called; this function is only not called in the first place when there are no centers of a current species in the ROI... it is thus still called even if there are no VALID centers in the ROI.

    Further, if there are no valid centers for the smallest slice size, then there can't be any valid centers for any larger slice sizes, and therefore there are no valid centers for the current center species period. In this case, both pvals_dens and pvals_pmf below will be None throughout. Furthermore, as the validity of the centers is unaffected by the neighbor species, then there would be no data for this center species as a whole (for any neighbor species), and therefore the entire species as a center could be dropped from any analysis plots. In this case, max_nbins_over_slices would remain 0.

    Of course, if there were some valid centers for the smallest slice size but not for larger slice sizes, then there would indeed be data in these analysis plots, and we should plot the data for this center species as usual. In this case, max_nbins_over_slices would not end up as 0.
    '''

    # Import relevant libraries
    import numpy as np
    import seaborn as sns

    # Calculate some variables
    nslices = len(data_by_slice)  # number of slices from the input data
    (vmin, vmax) = log_pval_range

    # Initialize some variables of interest
    any_valid_centers_over_slices = False  # whether there are any valid centers of all the slices
    nvalid_centers_holder = np.zeros((nslices,), dtype=np.uint32)  # number of valid centers in each slice

    # Define and populate the array holding the density P values for all the slices for the current center/neighbor combination
    pvals_dens = np.zeros((1, nslices, 2))
    if max_nbins_over_slices is None:
        max_nbins_over_slices = 0  # determine the maximum number of bins over all nslices slices, initializing this value to zero
    for islice in range(nslices):  # for every slice...
        if data_by_slice[islice][3] >= 1:  # as long as there is at least 1 valid center...
            pvals_dens[:, islice, :] = data_by_slice[islice][0][1:]  # get the left and right pvals of the islice-th slice; (2,)
            curr_nbins = data_by_slice[islice][1].shape[0]  # get the number of PMF bins for the current slice
            if curr_nbins > max_nbins_over_slices:  # update the maximum number of bins over all the slices
                max_nbins_over_slices = curr_nbins
            any_valid_centers_over_slices = True  # if we've gotten here, then there are SOME valid centers over all the slices
            nvalid_centers_holder[islice] = data_by_slice[islice][3]  # store this number of valid centers for the current slice
        else:
            pvals_dens[:, islice, :] = None

    # Define and populate the array holding the PMF P values for all the slices for the current center/neighbor combination
    if any_valid_centers_over_slices:  # as long as there are some valid centers over all the slices...
        pvals_pmf = np.zeros((max_nbins_over_slices, nslices, 2))
        for islice in range(nslices):  # for every slice...
            if data_by_slice[islice][3] >= 1:  # as long as there is at least 1 valid center...
                curr_nbins = data_by_slice[islice][1].shape[0]  # get the number of PMF bins for the current slice
                pvals_pmf[:curr_nbins, islice, :] = data_by_slice[islice][1][:, 1:]  # get the left and right pvals of the islice-th slice for all the current number of bins; (curr_nbins,2)
                if not calculate_empty_bin_pvals:  # optionally fill in the P values for higher j-values (higher-bins) whose PMFs are zero since these bins aren't even returned in the first place
                    pvals_pmf[curr_nbins:, islice, :] = None  # set the bins that are not set for the current slice (but exist here because over all slices there are bins this large) to None
                else:
                    nexpected, nvalid_centers = data_by_slice[islice][2:4]
                    pvals_pmf[curr_nbins:, islice, :] = calculate_PMF_metrics(np.zeros((max_nbins_over_slices, 1)), nexpected, nvalid_centers)[curr_nbins:, 1:]  # (nbins_remaining,2)
            else:
                pvals_pmf[:, islice, :] = None

    # Initialize the current figure by clearing it and all its axes
    fig.clf()
    ax = fig.subplots(nrows=2, ncols=2)

    # Create the four-axis figure, where the top row has the left and right density P values and the bottom row has the left and right PMF P values
    if any_valid_centers_over_slices:

        # Plot the log10 of the left density P values
        left_log_dens_pvals = np.log10(pvals_dens[:, :, 0])
        sns.heatmap(left_log_dens_pvals, vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[0, 0], cbar=True, yticklabels=True, square=True)
        ax[0, 0].set_title('log10(\"less\" density pvals)')
        ax[0, 0].set_xlabel('Slice')

        # Plot the log10 of the right density P values
        right_log_dens_pvals = np.log10(pvals_dens[:, :, 1])
        sns.heatmap(right_log_dens_pvals, vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[0, 1], cbar=True, yticklabels=True, square=True)
        ax[0, 1].set_title('log10(\"greater\" density pvals)')
        ax[0, 1].set_xlabel('Slice')

        # Plot the log10 of the left PMF P values
        left_log_pmf_pvals = np.log10(pvals_pmf[:, :, 0])
        sns.heatmap(left_log_pmf_pvals, vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[1, 0], cbar=True, yticklabels=yticklabels, square=square)
        ax[1, 0].set_title('log10(\"less\" PMF pvals)')
        ax[1, 0].set_xlabel('Slice')
        ax[1, 0].set_ylabel('Number of neighbors')

        # Plot the log10 of the right PMF P values
        right_log_pmf_pvals = np.log10(pvals_pmf[:, :, 1])
        sns.heatmap(right_log_pmf_pvals, vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[1, 1], cbar=True, yticklabels=yticklabels, square=square)
        ax[1, 1].set_title('log10(\"greater\" PMF pvals)')
        ax[1, 1].set_xlabel('Slice')
        ax[1, 1].set_ylabel('Number of neighbors')

        # Place a descriptive title on the figure
        if name is not None:
            fig.suptitle(name + '\nnum valid centers per slice: {}'.format(nvalid_centers_holder))

    else:
        # Initialize the return variables
        left_log_dens_pvals = None
        right_log_dens_pvals = None
        left_log_pmf_pvals = None
        right_log_pmf_pvals = None

    return(nvalid_centers_holder, left_log_dens_pvals, right_log_dens_pvals, left_log_pmf_pvals, right_log_pmf_pvals)


def get_max_nbins_for_center_neighbor_pair(data_by_slice, max_nbins_over_slices):
    '''
    Get the largest possible number of bins so that we can plot all the experiment's data in the same exact way (having the same number of rows for the PMF P value plots).

    These are calculated from the sizes of the PMF P value arrays and therefore the same sizes (nbins) as the calculated PMFs, where nbins = nneighbors.max() + 1.
    '''

    # Define and populate the array holding the density P values for all the slices for the current center/neighbor combination
    for islice in range(len(data_by_slice)):  # for every slice...
        if data_by_slice[islice][3] >= 1:  # as long as there is at least 1 valid center...
            curr_nbins = data_by_slice[islice][1].shape[0]  # get the number of PMF bins for the current slice
            if curr_nbins > max_nbins_over_slices:  # update the maximum number of bins over all the slices
                max_nbins_over_slices = curr_nbins

    return(max_nbins_over_slices)


def preprocess_dataset(format, input_datafile, coord_units_in_microns, images_to_analyze, min_coord_spacing=None, species_equivalents={}, mapping_dict={}, sep='\t', roi_width=None, overlap=0, phenotype_identification_tsv_file=None):
    """Given dataset-specific attributes, adhere to the format of the dataset required for processing in the SIP library time_cell_interaction_lib.py

    I should probably do away with species_equivalents and mapping_dict in favor of phenotype_identification_tsv_file.

    As of 3/13/24 this function should finally be completely unnecessary due to doing the equivalent more efficiently in 03a_Tool_parameter_selection.py.

    Args:
        format (str): Class name in dataset_formats.py describing the dataset, e.g., `Native`, `GMBSecondGeneration`, `OMAL`.
        input_datafile (str): Full pathname to the authors' original datafile.
        coord_units_in_microns (float): Factor by which to multiply the coordinates in the input datafile in order to obtain units of microns.
        min_coord_spacing (float, optional): Minimum coordinate spacing after conversion to microns. If set to None, this value is automatically calculated. Defaults to None.
        species_equivalents (dict, optional): Define equivalent species, which are identified by integer numbers, but sometimes differently-marked species should actually be the same species; specify that here. Unlike `mapping_dict`, this DOES affect the actual determination of the species and therefore the calculated data such as the P values. Species can be accessed via slices.plotting_map. Defaults to {}.
        mapping_dict (dict, optional): Dictionary for mapping markers to known species names in order to make the plots more clear. This does not affect the actual determination of the species; this is just labeling. Defaults to {}.
        sep (str, optional): Text separator in the input datafile. Defaults to '\t'.
        roi_width (float, optional): Desired ROI width in microns for patching up the input slide. If None, no patching is performed. Defaults to None.
        overlap (float, optional): If patching is requested (roi_width!=None), amount in microns you want consecutive ROIs to overlap. Note that for the SIP code, in order to utilize all data most efficiently, this should equal twice the density radius (e.g., the "thickness" value in main.ipynb). Defaults to 0.
        phenotype_identification_tsv_file (str, optional): Name of the TSV file generated by pasting into a text editor columns from Excel in the order specified in the read_csv() method below. Defaults to None.

    Returns:
        Instance of the class specified by the `format` parameter: Object describing the dataset.
    """

    # Import relevant libraries
    import dataset_formats
    import os

    # Determine the path of the pickle file holding the dataset object using the input data file name
    # pickle_path = os.path.splitext(input_datafile)[0] + '.pkl'
    # pickle_dir = os.path.dirname(pickle_path)
    pickle_dir = os.path.join('.', 'output', 'checkpoints')
    pickle_file = os.path.basename(os.path.splitext(input_datafile)[0] + '.pkl')

    # If the pickle file doesn't already exist, create the dataset object and store it in a pickle file
    if not os.path.exists(os.path.join(pickle_dir, pickle_file)):

        if not os.path.exists(pickle_dir):
            os.makedirs(pickle_dir)

        # If some apps still use the OMAL format, support it, even though it really means HALO
        if format == 'OMAL':
            format = 'HALO'

        # Obtain the class describing the dataset according to the `format` parameter
        dataset_class = getattr(dataset_formats, format)

        # Get an instance of the appropriate class
        dataset_obj = dataset_class(input_datafile, coord_units_in_microns, images_to_analyze, min_coord_spacing=min_coord_spacing, species_equivalents=species_equivalents, mapping_dict=mapping_dict, sep=sep, roi_width=roi_width, overlap=overlap, phenotype_identification_tsv_file=phenotype_identification_tsv_file)

        # Perform the actual format conversion/adherence
        dataset_obj.process_dataset()

        # Create a pickle file containing the dataset object
        make_pickle(dataset_obj, pickle_dir, pickle_file)

    # If the pickle file already exists, load it
    else:
        dataset_obj = load_pickle(pickle_dir, pickle_file)

    # Return the dataset object
    return dataset_obj


def refine_plotting_map_from_mapping_dict(df_plotting_map, mapping_dict):
    """Use the mapping dictionary containing labels/phenotypes with hyphen-separated entities to calculate the markers and circle sizes accordingly for the plotting map

    Args:
        df_plotting_map (Pandas dataframe): Dataframe version of the plotting map, which is traditionally a list
        mapping_dict (dict): Dictionary that takes positive markers to more descriptive labels/phenotypes

    Returns:
        Pandas dataframe: Updated plotting map
    """
    # Import relevant library
    import numpy as np

    # Use the mapping dictionary for more clearly labeling species to create a column of the plotting map dataframe containing these labels
    df_plotting_map['label'] = df_plotting_map['positive_markers'].apply(lambda x: replace_extra_characters(str(x)).replace(',', '/')).map(mapping_dict)

    # Get a dataframe of just the "simple" phenotypes, which only contain a single entity that could exist between hyphens in the labels
    simple_phenotype_markers = df_plotting_map.loc[df_plotting_map['label'].apply(lambda x: len(x.split('-'))) == 1, 'label'].reset_index().reset_index().set_index('label')

    # Add a column to the plotting map that contains a list of the simple phenotypes in decreasing frequency order
    df_plotting_map['label_sorted'] = df_plotting_map['label'].apply(lambda x: [y for y in simple_phenotype_markers.index if y in x.split('-')])

    # Use the simple phenotype dataframe to create a corresponding list of markers
    df_plotting_map['marker_list'] = df_plotting_map['label_sorted'].apply(lambda x: [simple_phenotype_markers.loc[y, 'index'] for y in x])

    # Extract just the first element if it's a list of length 1
    df_plotting_map['marker_ids'] = df_plotting_map['marker_list'].apply(lambda x: (x if len(x) > 1 else x[0]))

    # Create a list of circle sizes in decreasing frequency order
    df_plotting_map['circle_size_list'] = df_plotting_map['marker_list'].apply(lambda x: list(np.arange(start=len(x), stop=0, step=-1)))

    # Extract just the first element if it's a list of length 1
    df_plotting_map['circle_sizes'] = df_plotting_map['circle_size_list'].apply(lambda x: (x if len(x) > 1 else x[0]))

    # Drop the intermediate columns
    df_plotting_map = df_plotting_map.drop(['label_sorted', 'marker_list', 'circle_size_list'], axis='columns')

    return df_plotting_map


def replace_extra_characters(input_str):
    """Delete single quotes, square brackets, and spaces from a string

    Args:
        input_str (str): Arbitrary string

    Returns:
        str: String with extra characters removed
    """
    output_str = input_str.replace('\'', '').replace('[', '').replace(']', '').replace(' ', '')
    return output_str


def generate_case_slide_roi_contents(data, roi_colname='tag', num_last_chars_in_case_id_to_drop=1):
    """Print (and save to a dataframe) a summary of the cases, slides, ROIs, and ROI sizes in a dataset

    Call like, e.g.:

      contents = time_cell_interaction_lib.generate_case_slide_roi_contents(slices.data)

    Args:
        data (Pandas dataframe): Dataframe stored in slices.data

    Returns:
        Pandas dataframe: "Contents" dataframe containing all unique cases, slides, and ROIs and their widths/heights
    """

    # Import relevant library
    import pandas as pd

    # Rename the data dataframe
    df = data

    # Create a column containing the case ID
    if num_last_chars_in_case_id_to_drop != 0:
        df['Case ID'] = df['Slide ID'].apply(lambda x: x.split('-')[0][:-num_last_chars_in_case_id_to_drop])
    else:
        df['Case ID'] = df['Slide ID'].apply(lambda x: x.split('-')[0][:])

    # Initialize the variables holding the total numbers of cases, slides, and ROIs
    ncases_tot, nslides_tot, nrois_tot = 0, 0, 0

    # Initialize a list holding all the data
    data_holder = []

    # Determine the unique cases in the dataset and loop over them
    unique_cases = df['Case ID'].unique()
    ncases = 0
    for unique_case in unique_cases:

        # Update the total and current case count
        ncases_tot = ncases_tot + 1
        ncases = ncases + 1

        # Get the rows in the dataframe corresponding to the current case
        case_rows = df['Case ID'] == unique_case

        # Determine the unique slides in the current case and loop over them
        unique_slides = df.loc[case_rows, 'Slide ID'].unique()
        nslides = 0
        for unique_slide in unique_slides:

            # Update the total and current slide count
            nslides_tot = nslides_tot + 1
            nslides = nslides + 1

            # Get the rows in the dataframe corresponding to the current slide in the current case
            slide_rows = case_rows & (df['Slide ID'] == unique_slide)

            # Determine the unique ROIs in the current slide in the current case and loop over them
            unique_rois = df.loc[slide_rows, roi_colname].unique()
            nrois = 0
            for unique_roi in unique_rois:

                # Update the total and current ROI count
                nrois_tot = nrois_tot + 1
                nrois = nrois + 1

                # Print out the case and slide information for the current ROI
                print('Case (curr={}, tot={}): {}; slide (curr={}, tot={}): {}; ROI (curr={}, tot={}): {}'.format(ncases, ncases_tot, unique_case, nslides, nslides_tot, unique_slide, nrois, nrois_tot, unique_roi))

                # Save the current set of data
                data_holder.append([ncases, ncases_tot, unique_case, nslides, nslides_tot, unique_slide, nrois, nrois_tot, unique_roi])

    # At the end, print out the total number of cases, slides, and ROIs
    print('Totals: cases={}, slides={}, rois={}'.format(ncases_tot, nslides_tot, nrois_tot))

    # Create a Pandas dataframe to save all the data so far
    contents = pd.DataFrame(data=data_holder, columns=['ncases', 'ncases_tot', 'unique_case', 'nslides', 'nslides_tot', 'unique_slide', 'nrois', 'nrois_tot', 'unique_roi'])

    # Group the coordinates by ROI name
    grouped = df[[roi_colname, 'Cell X Position', 'Cell Y Position']].groupby(by=roi_colname)

    # Set the original contents index name to "orig_index" so that we can later get that index back
    contents.index.name = 'orig_index'

    # Push orig_index of contents into the columns and set the ROI names as the new index
    contents_by_roi = contents.reset_index().set_index('unique_roi')

    # Get the width and height of each ROI
    roi_sizes_by_roi = (grouped.max() - grouped.min()).rename({'Cell X Position': 'width', 'Cell Y Position': 'height'}, axis='columns')

    # To address decimated case having slightly different P values, include the original x and y ranges of each ROI (prior to decimation)
    roi_coord_mins_by_roi = (grouped.min()).rename({'Cell X Position': 'x_min_prior_to_decimation', 'Cell Y Position': 'y_min_prior_to_decimation'}, axis='columns')
    roi_coord_maxs_by_roi = (grouped.max()).rename({'Cell X Position': 'x_max_prior_to_decimation', 'Cell Y Position': 'y_max_prior_to_decimation'}, axis='columns')

    # Concatenate the two dataframes by ROI name
    contents_new = pd.concat([contents_by_roi, roi_sizes_by_roi, roi_coord_mins_by_roi, roi_coord_maxs_by_roi], axis='columns')

    # Reset the name of the index to the corresponding column name in the original contents dataframe: unique_roi
    contents_new.index.name = contents_by_roi.index.name

    # Push the index into the columns and set the original index as the new index, finally sorting by this original index
    contents_new = contents_new.reset_index().set_index('orig_index').sort_index()

    # Re-define and return the contents dataframe with the ROI sizes appended in columns
    return contents_new


def calculate_metrics_for_roi(args_as_single_tuple):
    """Calculate the metrics for a single ROI

    Args:
        args_as_single_tuple (tuple): Tuple of arguments to be unpacked below, in this format so that the metrics can be calculated using the multiprocessing library in the traditional way. See the calculate_metrics() method of the TIMECellInteraction class above for more details. Arguments are: pickle_dir, nslices, thickness, min_coord_spacing, data, all_species_list, nall_species, df_data_by_roi, keep_unnecessary_calculations, roi_index
    """

    # Import relevant modules
    import os
    import numpy as np
    import time
    import contextlib
    from datetime import datetime
    import tci_squidpy_supp_lib

    # Unpack the arguments
    pickle_dir, nslices, thickness, n_neighs, radius_instead_of_knn, min_coord_spacing, all_species_list, nall_species, do_logging, use_analytical_significance, df_data_by_roi, keep_unnecessary_calculations, n_jobs, roi_index = args_as_single_tuple

    # Constants which I can later turn into a parameter if desired
    my_seed = 42
    z_hardcode = 0  # this shouldn't matter anyway since I don't do anything with Z scores

    # Determine the pickle filename using the ROI index
    pickle_file = 'calculated_metrics-roi_index_{:06}.pkl'.format(roi_index)
    log_file = 'calculated_metrics-roi_index_{:06}.log'.format(roi_index)

    with (open(file=os.path.join(pickle_dir, log_file), mode='wt') if do_logging else contextlib.nullcontext()) as log_file_handle:

        # Determine the ROI name from the ROI index
        uroi = df_data_by_roi.loc[roi_index, 'unique_roi']

        if do_logging:
            log_file_handle.write('ROI {:06d} (split 00, {}): ROI processing started at {}\n'.format(roi_index, uroi, datetime.fromtimestamp(datetime.timestamp(datetime.now())).strftime("%Y-%m-%d, %H:%M:%S")))

        # If the pickle file doesn't already exist...
        if not os.path.exists(os.path.join(pickle_dir, pickle_file)):

            # Save the starting time
            start_time = time.time()

            # Note what we're doing
            print('Calculating metrics for ROI {} (ROI index {})'.format(uroi, roi_index))

            if do_logging:
                log_file_handle.write('ROI {:06d} (split 01, {}): pickle file {}/{} does not already exist for the ROI\n'.format(roi_index, uroi, pickle_dir, pickle_file))

            # Get the needed ROI data
            x_roi = df_data_by_roi.loc[roi_index, 'x_roi']
            y_roi = df_data_by_roi.loc[roi_index, 'y_roi']
            species_roi = df_data_by_roi.loc[roi_index, 'species_roi']
            roi_x_range_prior_to_decimation = np.array([df_data_by_roi.loc[roi_index, 'x_min_prior_to_decimation'], df_data_by_roi.loc[roi_index, 'x_max_prior_to_decimation']])
            roi_y_range_prior_to_decimation = np.array([df_data_by_roi.loc[roi_index, 'y_min_prior_to_decimation'], df_data_by_roi.loc[roi_index, 'y_max_prior_to_decimation']])
            unique_species_in_roi = np.unique(species_roi)
            num_unique_species_in_roi = len(unique_species_in_roi)
            coords_roi = np.c_[x_roi, y_roi]

            # Run some checks
            roi_x_range, roi_y_range, roi_min_coord_spacing = roi_checks_and_output(x_roi, y_roi, do_printing=False)  # note I don't believe roi_min_coord_spacing is ever actually used for anything other than saving the value to a text file
            assert (roi_x_range == df_data_by_roi.loc[roi_index, 'x_range']).all(), 'ERROR: ROI x range is not consistent ({} != {})'.format(roi_x_range, df_data_by_roi.loc[roi_index, 'x_range'])
            assert (roi_y_range == df_data_by_roi.loc[roi_index, 'y_range']).all(), 'ERROR: ROI y range is not consistent ({} != {})'.format(roi_y_range, df_data_by_roi.loc[roi_index, 'y_range'])
            # if (not ([x_roi.min(),x_roi.max()]==roi_x_range)) or (not ([y_roi.min(),y_roi.max()]==roi_y_range)) or (not (0.5==roi_min_coord_spacing)):
            if (not ([x_roi.min(), x_roi.max()] == roi_x_range)) or (not ([y_roi.min(), y_roi.max()] == roi_y_range)):
                print('ERROR: A basic check failed')
                # print(roi_x_range, roi_y_range, roi_min_coord_spacing)
                print(roi_x_range, roi_y_range)
                exit()

            # Define the arrays of interest holding the metrics and other data
            real_data = np.empty((nall_species, nall_species, nslices), dtype=object)
            sim_data = np.empty((nall_species, nall_species, nslices), dtype=object)

            if use_analytical_significance:

                # For every center and neighbor species in the entire experiment...
                for icenter_spec, center_species in enumerate(all_species_list):
                    for ineighbor_spec, neighbor_species in enumerate(all_species_list):

                        # If the current center and neighbor species exist in the current ROI...
                        if (center_species in unique_species_in_roi) and (neighbor_species in unique_species_in_roi):

                            if do_logging:
                                log_file_handle.write('ROI {:06d} (split 02, {}): center {} and neighbor {} both exist in the ROI to any extent\n'.format(roi_index, uroi, center_species, neighbor_species))

                            # Determine the coordinates of the current center and neighbor species
                            coords_centers = coords_roi[species_roi == center_species, :]
                            coords_neighbors = coords_roi[species_roi == neighbor_species, :]

                            # For every radius/slice...
                            for islice in range(nslices):

                                # Define the inner and outer radii of the current slice
                                small_rad = islice * thickness
                                large_rad = (islice + 1) * thickness

                                # Count the neighbors, calculate the PMFs, and from these determine the P values of interest for the single set of real data
                                density_metrics_real, pmf_metrics_real, nexpected_real, nvalid_centers_real, coords_centers_real, coords_neighbors_real, valid_centers_real, edges_real, npossible_neighbors_real, roi_area_used_real, slice_area_used_real = \
                                    calculate_metrics_from_coords(min_coord_spacing, input_coords=(coords_centers, coords_neighbors), neighbors_eq_centers=(neighbor_species == center_species), nbootstrap_resamplings=0, rad_range=(small_rad, large_rad), use_theoretical_counts=False, roi_edge_buffer_mult=1, roi_x_range=roi_x_range_prior_to_decimation, roi_y_range=roi_y_range_prior_to_decimation, silent=False, log_file_data=(log_file_handle, roi_index, uroi, center_species, neighbor_species), keep_unnecessary_calculations=keep_unnecessary_calculations)

                                # Save the results, plus some other data, into a primary array of interest
                                real_data[icenter_spec, ineighbor_spec, islice] = (density_metrics_real, pmf_metrics_real, nexpected_real, nvalid_centers_real, coords_centers_real, coords_neighbors_real, valid_centers_real, edges_real, npossible_neighbors_real, roi_area_used_real, slice_area_used_real, center_species, neighbor_species, small_rad, large_rad, islice)

                                if keep_unnecessary_calculations:

                                    # Count the neighbors, calculate the PMFs, and from these determine the P values of interest for a single set of simulated data with the same properties as the real data
                                    density_metrics_sim, pmf_metrics_sim, nexpected_sim, nvalid_centers_sim, coords_centers_sim, coords_neighbors_sim, valid_centers_sim, edges_sim, npossible_neighbors_sim, roi_area_used_sim, slice_area_used_sim = \
                                        calculate_metrics_from_coords(min_coord_spacing, input_coords=None, neighbors_eq_centers=(neighbor_species == center_species), ncenters_roi=coords_centers.shape[0], nneighbors_roi=coords_neighbors.shape[0], nbootstrap_resamplings=0, rad_range=(small_rad, large_rad), use_theoretical_counts=False, roi_edge_buffer_mult=1, roi_x_range=roi_x_range_prior_to_decimation, roi_y_range=roi_y_range_prior_to_decimation, silent=False, log_file_data=(log_file_handle, roi_index, uroi, center_species, neighbor_species), keep_unnecessary_calculations=keep_unnecessary_calculations)

                                    # Save the results, plus some other data, into a primary array of interest
                                    sim_data[icenter_spec, ineighbor_spec, islice] = (density_metrics_sim, pmf_metrics_sim, nexpected_sim, nvalid_centers_sim, coords_centers_sim, coords_neighbors_sim, valid_centers_sim, edges_sim, npossible_neighbors_sim, roi_area_used_sim, slice_area_used_sim, center_species, neighbor_species, small_rad, large_rad, islice)

                                else:

                                    sim_data[icenter_spec, ineighbor_spec, islice] = (None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None)

                        else:

                            if do_logging:
                                log_file_handle.write('ROI {:06d} (split 02, {}): center {} ({}) or neighbor {} ({}) does not exist in the ROI to any extent\n'.format(roi_index, uroi, center_species, ('exists' if (center_species in unique_species_in_roi) else 'does not exist'), neighbor_species, ('exists' if (neighbor_species in unique_species_in_roi) else 'does not exist')))

            else:  # use Squidpy's numerical method

                # The numerical method exits (gracefully) if there's only one species (i.e., type of cell) present in the ROI
                if (num_unique_species_in_roi == 1):
                    print('NOTE: There is only a single species ({}) in the current ROI, which is undefined for the numerical method; no P values have been calculated'.format(unique_species_in_roi))
                    if do_logging:
                        log_file_handle.write('ROI {:06d} ({}): only a single species ({}) exists in the current ROI, which is undefined for the numerical method; no P values have been calculated\n'.format(roi_index, uroi, unique_species_in_roi))

                # If the number of coordinates is 2, Squidpy functionality sometimes dies (NOT gracefully), e.g., sq.pl.spatial_scatter() dies depending on the species labels due to its sometimes incorrect setting of color_vector, which is usually manifested when there are 2-3 coordinates. The error in spatial_scatter() is generally something like: "ValueError: 'c' argument must be a color, a sequence of colors, or a sequence of numbers, not ['aa' 'bb']". Note the neighborhood enrichment also dies with something like: "TypeError: Expected `adata.obs['species']` to be `categorical`, found `string`"
                elif len(coords_roi) <= 2:
                    print('NOTE: There is more than a single unique species in the ROI but the total number of cells ({}) is 2 or fewer, which can cause non-graceful issues in Squidpy; no P values have been calculated'.format(len(coords_roi)))
                    if do_logging:
                        log_file_handle.write('ROI {:06d} ({}): there is more than a single unique species in the ROI but the total number of cells ({}) is 2 or fewer, which can cause non-graceful issues in Squidpy; no P values have been calculated\n'.format(roi_index, uroi, unique_species_in_roi))

                else:

                    # Create and store the name of the directory in which to store Squidpy's scatter plots and heatmaps
                    squidpy_dir = os.path.join('.', 'output', 'images', 'squidpy')
                    if not os.path.exists(squidpy_dir):
                        os.makedirs(squidpy_dir)
                    image_path_prefix = os.path.join(squidpy_dir, 'roi_index-{}__roi_name-{}-'.format(roi_index, uroi))
                    
                    if do_logging:
                        log_file_handle.write('ROI {:06d} (split 02, {}): calculating P values using SquidPy for all centers and neighbors\n'.format(roi_index, uroi))

                    # Calculate the left and right P values using Squidpy's neighborhood enrichment method
                    pvals_left, pvals_right, heatmap_rowcol_labels = tci_squidpy_supp_lib.calculate_squidpy_pvals(coordinates=coords_roi, labels=species_roi, seed=my_seed, label_name='species', circle_size=50, image_path_prefix=image_path_prefix, dpi=150, radius=thickness, n_neighs=n_neighs, radius_instead_of_knn=radius_instead_of_knn, annotate_heatmap=True, print_heatmap_data=False, close_figs=True, n_jobs=n_jobs)
                
                    # Conform to the data storage format already set up in the codebase: for every center and neighbor species in the entire experiment...
                    islice = 0
                    for icenter_spec, center_species in enumerate(all_species_list):
                        for ineighbor_spec, neighbor_species in enumerate(all_species_list):

                            # If the current center and neighbor species exist in the current ROI...
                            if (center_species in unique_species_in_roi) and (neighbor_species in unique_species_in_roi) and (pvals_left is not None):

                                if do_logging:
                                    log_file_handle.write('ROI {:06d} (split 02, {}): center {} and neighbor {} both exist in the ROI to any extent\n'.format(roi_index, uroi, center_species, neighbor_species))

                                # Determine the indexes of the Squidpy output corresponding to the current center and neighbor species
                                sp_heatmap_idx_center = heatmap_rowcol_labels.index(center_species)
                                sp_heatmap_idx_neighbor = heatmap_rowcol_labels.index(neighbor_species)

                                # Determine the corresponding P values
                                pval_left = pvals_left[sp_heatmap_idx_center, sp_heatmap_idx_neighbor]
                                pval_right = pvals_right[sp_heatmap_idx_center, sp_heatmap_idx_neighbor]

                                # Determine the coordinates of the current center and neighbor species
                                nvalid_centers = (species_roi == center_species).sum()  # in this case, set the number of valid centers to the total number of centers

                                # For the current center/neighbor species pair, appropriately store the left and right P values
                                real_data[icenter_spec, ineighbor_spec, islice] = (np.array([z_hardcode, pval_left, pval_right]), None, None, nvalid_centers, None, None, None, None, None, None, None, center_species, neighbor_species, 0, thickness, islice)
                                
                                # Duplicate the real P values on the simulated data, which is meaningless in this case because Squidpy doesn't do a simulation like I do, which I don't need anyway but am doing it for posterity and to aid debugging
                                # sim_data[icenter_spec, ineighbor_spec, islice] = real_data[icenter_spec, ineighbor_spec, islice]
                                sim_data[icenter_spec, ineighbor_spec, islice] = (None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None)

                            else:

                                if do_logging:
                                    log_file_handle.write('ROI {:06d} (split 02, {}): center {} ({}) or neighbor {} ({}) does not exist in the ROI to any extent, or Squidpy failed for the ROI\n'.format(roi_index, uroi, center_species, ('exists' if (center_species in unique_species_in_roi) else 'does not exist'), neighbor_species, ('exists' if (neighbor_species in unique_species_in_roi) else 'does not exist')))

            # Save the data to a data structure
            roi_data_item = {'roi_min_coord_spacing': roi_min_coord_spacing, 'unique_species_in_roi': unique_species_in_roi, 'real_data': real_data, 'sim_data': sim_data, 'all_species_list': all_species_list, 'roi_x_range': roi_x_range_prior_to_decimation, 'roi_y_range': roi_y_range_prior_to_decimation}

            # Run a check on the unique species in the ROI
            if not (set(df_data_by_roi.loc[roi_index, 'spec2plot_roi']) == set(unique_species_in_roi)):
                print('spec2plot_roi:', df_data_by_roi.loc[roi_index, 'spec2plot_roi'])
                print('unique_species_in_roi:', unique_species_in_roi)
            assert set(df_data_by_roi.loc[roi_index, 'spec2plot_roi']) == set(unique_species_in_roi), 'ERROR: spec2plot_roi does not equal unique_species_in_roi!'

            # Create a pickle file saving the data that we just calculated
            make_pickle(roi_data_item, pickle_dir, pickle_file)

            # Output the metrics calculation time for the current ROI
            duration = time.time() - start_time
            print('Metrics calculation for ROI {} (ROI ID {}) took {} seconds'.format(uroi, roi_index, duration))

        else:

            # The pickle file already exists
            print('The pickle file {} in directory {} already exists'.format(pickle_file, pickle_dir))

            if do_logging:
                log_file_handle.write('ROI {:06d} (split 01, {}): pickle file {}/{} already exists for the ROI\n'.format(roi_index, uroi, pickle_dir, pickle_file))


def save_figs_and_corresp_data_for_roi(args_as_single_tuple):
    """Save the ROI and P value figures and the corresponding data for a single ROI; this way we can parallelize over all ROIs.

    Args:
        args_as_single_tuple (tuple): tuple containing all necessary data (see "Unpack the arguments" block below)
    """

    # Import relevant modules
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import time

    # Unpack the arguments
    pickle_dir, plotting_map, roi_figsize, pval_figsize, colors, marker_size_step, default_marker_size, roi_dpi, mapping_dict, coord_units_in_microns, alpha, edgecolors, yaxis_dir, webpage_dir, plot_real_data, nall_species, species_names, log_pval_range, calculate_empty_bin_pvals, max_nbins, square, yticklabels, all_species_list, save_individual_pval_plots, pval_dpi, df, roi_index = args_as_single_tuple

    # Determine the pickle filename using the ROI index
    pickle_file = 'figure_data-{}-roi_index_{:06}.pkl'.format(('real' if plot_real_data else 'simulated'), roi_index)

    # If the pickle file doesn't already exist...
    if not os.path.exists(os.path.join(pickle_dir, pickle_file)):

        # Save the starting time
        start_time = time.time()

        # Assign the data for the current ROI to useful variables
        roi_name = df.loc[roi_index, 'roi_name']
        roi_data = df.loc[roi_index, 'roi_data']
        num_rois_in_slide = df.loc[roi_index, 'num_rois_in_slide']

        # Note what we're doing
        print('Saving figures and corresponding data for ROI {} (ROI index {})'.format(roi_name, roi_index))

        # Collect the data needed to plot the ROIs, plot the ROI, and save it to disk
        x_range = np.array(roi_data[0][0])
        y_range = np.array(roi_data[0][1])
        roi_min_coord_spacing = roi_data[0][2]
        x_roi = roi_data[0][5]
        y_roi = roi_data[0][6]
        species_roi = roi_data[0][7]
        uroi = roi_name
        spec2plot_roi = [x[0] for x in plotting_map if x[0] in roi_data[0][3]]

        # Define the figures to use for plotting the ROIs and the P values
        fig_roi = plt.subplots(figsize=roi_figsize)[0]
        fig_pvals = plt.subplots(nrows=2, ncols=2, figsize=pval_figsize)[0]

        # Plot and save the current ROI
        plot_roi(fig_roi, spec2plot_roi, species_roi, x_roi, y_roi, plotting_map, colors, x_range, y_range, uroi, marker_size_step, default_marker_size, roi_dpi, mapping_dict, coord_units_in_microns, filepath=None, do_plot=True, alpha=alpha, edgecolors=edgecolors, yaxis_dir=yaxis_dir)
        roi_fig_filename = 'roi_{}.png'.format(roi_name)
        roi_fig_pathname = os.path.join(webpage_dir, roi_fig_filename)
        fig_roi.savefig(roi_fig_pathname, dpi=roi_dpi, bbox_inches='tight')

        # Close the figure for the ROI
        plt.close(fig_roi)

        # Determine which dataset to plot, each of which is (nall_species, nall_species, nslices)
        # The elements of the dataset ndarray are either a tuple (if both the center and neighbor species are in the ROI) or None
        if plot_real_data:
            data_to_plot = roi_data[1]  # this is the single set of real data
        else:
            data_to_plot = roi_data[2]  # this is the corresponding single set of simulated data

        # Determine whether both the center and neighbor species are in the dataset
        center_and_neighbor_species_in_dataset = ((data_to_plot != None).sum(axis=2)) > 0  # (nall_species, nall_species)... since we want to do element-wise comparisons here, don't listen to linting when it says the right way to do the comparison is "data_to_plot is not None"

        # Define the data we want to save for the ROI
        filedata_for_roi = []

        # For all combinations of centers and neighbors...
        for icenter_spec in range(nall_species):
            center_name = species_names[icenter_spec]
            for ineighbor_spec in range(nall_species):
                neighbor_name = species_names[ineighbor_spec]

                # If data exist for at least one of the slices for the current center/neighbor combination...
                if center_and_neighbor_species_in_dataset[icenter_spec, ineighbor_spec]:

                    # Create and save the P value figure containing four heatmaps as well as their corresponding data
                    print(icenter_spec, ineighbor_spec)
                    figure_title = '{} data\nSlide/ROI: {}\ncenter={}, neighbor={}'.format(('Real' if plot_real_data else 'Simulated'), roi_name, center_name, neighbor_name)
                    nvalid_centers_per_slice, left_log_dens_pvals, right_log_dens_pvals, left_log_pmf_pvals, right_log_pmf_pvals = plot_pvals(fig_pvals, data_to_plot[icenter_spec, ineighbor_spec, :], log_pval_range, name=figure_title, calculate_empty_bin_pvals=calculate_empty_bin_pvals, max_nbins_over_slices=max_nbins, square=square, yticklabels=yticklabels)
                    if nvalid_centers_per_slice.sum(axis=0) == 0:
                        print('NOTE: Even though there are centers of species {} present in the ROI, none of them are valid; we can eliminate this species as a center (for any neighbor) from any analyses!'.format(center_name))
                    else:
                        pvals_fig_filename = 'pvals_{}_center-{}_neighbor-{}.png'.format(roi_name, all_species_list[icenter_spec], all_species_list[ineighbor_spec])
                        pvals_fig_pathname = os.path.join(webpage_dir, pvals_fig_filename)
                        if save_individual_pval_plots:  # make this an option since it takes significant time in aggregate
                            fig_pvals.savefig(pvals_fig_pathname, dpi=pval_dpi, bbox_inches='tight')
                        filedata_for_roi.append({
                            'roi_fig_pathname': roi_fig_pathname, 'nrois_in_slide': num_rois_in_slide, 'roi_name': roi_name, 'unique_species_in_roi': spec2plot_roi, 'roi_x_range': x_range, 'roi_y_range': y_range, 'roi_spacing': roi_min_coord_spacing,
                            'pvals_fig_pathname': pvals_fig_pathname, 'center_species_id': all_species_list[icenter_spec], 'neighbor_species_id': all_species_list[ineighbor_spec], 'center_species_name': species_names[icenter_spec], 'neighbor_species_name': species_names[ineighbor_spec], 'nvalid_centers_per_slice': nvalid_centers_per_slice, 'left_log_dens_pvals': left_log_dens_pvals, 'right_log_dens_pvals': right_log_dens_pvals, 'left_log_pmf_pvals': left_log_pmf_pvals, 'right_log_pmf_pvals': right_log_pmf_pvals
                        })

        # Close the figure for the P values
        plt.close(fig_pvals)

        # Create a pickle file saving the data that we just calculated
        make_pickle(filedata_for_roi, pickle_dir, pickle_file)

        # Output the metrics calculation time for the current ROI
        duration = time.time() - start_time
        print('Figures and corresponding data saving for ROI {} (ROI ID {}) took {} seconds'.format(uroi, roi_index, duration))

    # If the pickle file DOES already exist...
    else:

        # The pickle file already exists
        print('The pickle file {} in directory {} already exists'.format(pickle_file, pickle_dir))


def check_roi_order(data, location_string='', roi_col_name='tag'):
    """Check whether the ROIs i.e. tags in a dataframe are sorted.

    Args:
        data (Pandas dataframe): Dataframe containing a roi_col_name column.
        location_string (str, optional): String describing the location in the code in order to help debug. Defaults to ''.
    """
    unique_rois = data[roi_col_name].unique()  # note that this unique() method does not sort the results
    unique_rois_sorted = unique_rois.copy()
    unique_rois_sorted.sort()  # this sorts the ROI names
    print('----')
    print('Location: {}'.format(location_string))
    print('  Dataframe length: {}'.format(len(data)))
    print('  Number of unique ROIs: {}'.format(len(unique_rois)))
    print('  ROIs in sorted order? {}'.format((unique_rois_sorted == unique_rois).sum() == len(unique_rois)))


def create_colorbar(colormap_name='rocket', figsize=(2, 12), vmin=-50, vmax=0, label='log(P)', save_filename=None):
    """Plot a colorbar and return the object.

    In order to save the colorbar to a file, you can call like, e.g.:

      import time_cell_interaction_lib as tci
      tci.create_colorbar(colormap_name='rocket', figsize=(2, 12), vmin=-50, vmax=0, label='log(P)', save_filename='colorbar.png')

    Args:
        colormap_name (str, optional): Colormap name. Defaults to 'rocket', the default for Seaborn's heatmap.
        figsize (tuple, optional): Figure size for the colorbar. Defaults to (1, 6).
        vmin (int, optional): Minimum data value for the colorbar. Defaults to -50.
        vmax (int, optional): Maximum data value for the colorbar. Defaults to 0.
        label (str, optional): Colobar label. Defaults to 'log(P)'.
        save_filename (str, optional): Name of the file to which you want to save the image. Defaults to None, which means not to save the colorbar as an image on the filesystem.

    Returns:
        matplotlib.colorbar.Colorbar: Colorbar object.
    """

    # Import relevant libraries
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import seaborn  # needed to get e.g. the "rocket" colormap which is Seaborn's default for heatmaps

    # Sample colormap to a specified amount
    x_bit_color = 16  # default appears to be 8
    lut = 2 ** x_bit_color  # same name as argument to get_cmap() below

    # Get the colormap (rocket is the default for Seaborn's heatmap)
    # Note this can be obtained from an existing heatmap using cmap = hm.collections[-1].colorbar.cmap
    cmap = mpl.cm.get_cmap(colormap_name, lut=lut)

    # Get a figure/axis for the colorbar
    fig, ax = plt.subplots(figsize=figsize)

    # Get the appropriate normalization object
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    # Plot the colorbar
    colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax, orientation='vertical', label=label)

    # Save the colorbar to a file if so desired
    if save_filename is not None:
        ax.set_ylabel(ax.get_ylabel(), fontsize=24)
        plt.yticks(fontsize=16)
        fig.savefig(save_filename, facecolor='white', bbox_inches='tight')

    # Return the colorbar object
    return colorbar


def get_rgba_for_data_on_colormap(log_pval=-140, vmin=-200, vmax=0, colormap_name='rocket', print_byte_value=False):
    """Get the RGBA data for a datapoint on a particular colormap. The following block is modified from https://stackoverflow.com/questions/42125433/seaborn-heatmap-get-array-of-color-codes-values.

    **NOTE**: This is probably deprecated and should be replaced by, e.g., get_properly_interpolated_rgba_from_data(-0.33, vmin=-200, vmax=0, N=2**16, colormap_name='rocket', print_byte_value=False).

    Args:
        log_pval (int, optional): The data. Defaults to -140.
        vmin (int, optional): Minimum data value for the colormap. Defaults to -200.
        vmax (int, optional): Maximum data value for the colormap. Defaults to 0.
        colormap_name (str, optional): Colormap name. Defaults to 'rocket', the default for Seaborn's heatmap.
        print_byte_value (bool, optional): Whether to print a 0-255-ranged RGBA value for testing. Defaults to False.

    Returns:
        tuple: Desired RGBA data as a four-element tuple
    """

    # Import relevant libraries
    import matplotlib as mpl
    import numpy as np

    # Sample colormap to a specified amount
    x_bit_color = 16  # default appears to be 8
    lut = 2 ** x_bit_color  # same name as argument to get_cmap() below

    # Get the colormap (rocket is the default for Seaborn's heatmap)
    # Note this can be obtained from an existing heatmap using cmap = hm.collections[-1].colorbar.cmap
    cmap = mpl.cm.get_cmap(colormap_name, lut=lut)

    print(len(cmap.colors))

    # Get the appropriate normalization object
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    # Get the needed rgba values
    rgba_values = cmap(norm(log_pval))

    # If we want to print out the RGBA values for testing, do so (0-255)
    if print_byte_value:
        print((np.array(rgba_values) * 255).round())

    # Return the 0-1-ranged RGBA values that we wanted
    return rgba_values


def plot_simple_heatmap(log_pval=-140, vmin=-200, vmax=0):
    """Plot a simple, one-element heatmap using Seaborn.

    Args:
        log_pval (int, optional): The data. Defaults to -140.
        vmin (int, optional): Minimum data value for the colormap. Defaults to -200.
        vmax (int, optional): Maximum data value for the colormap. Defaults to 0.
    """

    # Import relevant libraries
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np

    # Get an axis on which to plot the heatmap
    _, ax = plt.subplots()

    # Format the data as a 2D numpy array
    data = np.array([log_pval]).reshape((1, 1))

    # Create the one-element heatmap using Seaborn, including the colorbar
    sns.heatmap(data, vmin=vmin, vmax=vmax, ax=ax, cbar=True)


def undo_patching_overlaps_and_decompounding(df):
    """From the data dataframe in the time_cell_interaction.TIMECellInteraction class, first remove duplicates that were created due to overlapping patching, and then remove duplicates that were created due to de-compounding the species.

    Call like:
        depatched, depatched_and_compounded = undo_patching_overlaps_and_decompounding(slices.data)

    Args:
        df (Pandas dataframe): Data dataframe created by instantiation of an object of the TIMECellInteraction class

    Returns:
        Pandas dataframe: Dataframe with patching removed
        Pandas dataframe: Dataframe with both patching and decompounding removed
    """

    # Import relevant libraries
    import pandas as pd

    # Everything will be done by slide so get the unique slides in the entire dataset
    unique_slides = df['Slide ID'].unique()

    # Initialize arrays that will hold new dataframes to be concatenated later
    depatched = []
    depatched_and_compounded = []

    # For each slide...
    for unique_slide in unique_slides:

        # Get the data for just the current slide
        curr_df = df[df['Slide ID'] == unique_slide]

        # For the objects in the slide with the exact same coordinates AND species ID, keep only the first, thereby eliminating duplicated objects due to overlapping patches
        curr_depatched = curr_df.groupby(by=['Cell X Position', 'Cell Y Position', 'Species int']).first().reset_index()

        # Save this for later
        depatched.append(curr_depatched)

        # Then, group by just the coordinates, which if compounding was originally not allowed, their should be duplicates due to multiple positive surface markers per object
        grouped = curr_depatched.groupby(by=['Cell X Position', 'Cell Y Position'])

        # Get the first row of data within each of these groups
        curr_compounded = grouped.first()

        # Then, using the grouped data, simply add up the species IDs (which must be pure powers of two) in order to get the compounded spedcies ID, for which the "Species string" and "Phenotype ..." columns still apply since we never modified these as it was unnecessary. Save these species IDs as the species ID for each of the groups
        curr_compounded['Species int'] = grouped['Species int'].sum()

        # Push the indexes into columns
        curr_compounded = curr_compounded.reset_index()

        # Save this for later
        depatched_and_compounded.append(curr_compounded)

    # Combine the data for all the slides
    depatched = pd.concat(depatched)
    depatched_and_compounded = pd.concat(depatched_and_compounded)

    # Print out the species counts in all three of these views of the dataset
    print('Species counts for patched and decompounded data:')
    print(df['Species int'].value_counts())
    print('Total: {}'.format(len(df)))
    print('')
    print('Species counts for de-patched but still decompounded data:')
    print(depatched['Species int'].value_counts())
    print('Total: {}'.format(len(depatched)))
    print('')
    print('Species counts for de-patched and compounded data:')
    print(depatched_and_compounded['Species int'].value_counts())
    print('Total: {}'.format(len(depatched_and_compounded)))

    # Return the "reduced" dataframes
    return depatched, depatched_and_compounded


def plot_just_rois(df, plotting_map, num_colors, webpage_dir, mapping_dict, coord_units_in_microns, slide_id_col='Slide ID', roi_id_col='tag.x', x_coord_col='Cell X Position', y_coord_col='Cell Y Position', plot_real_data=True, roi_figsize=(6, 4), marker_size_step=0.80, alpha=1, roi_dpi=200, edgecolors='k', default_marker_size_fac=1, yaxis_dir=1, boxes_to_plot0=None, filename_suffix='', pval_params=None, title_suffix=None, plot_just_single_roi_outline=False, webpage_dir_addon=None, nworkers=1, use_multiprocessing=True):
    '''
    This is based on save_figs_and_corresp_data() and includes just the pieces of that function that plot the ROI, as opposed to any metrics such as P values.

    * alpha is the transparency for the circles in the ROI plots (0 is fully transparent, 1 is fully opaque)
    * marker_size_step=0.80 means the radius should be 80% larger for cells plotted behind other cells

    Sample calls:

      Sample #1:
      tci.plot_just_rois(slices.data, slices.plotting_map, slices.num_colors, slices.webpage_dir, slices.mapping_dict, slices.coord_units_in_microns, slide_id_col='Slide ID', roi_id_col='tag', x_coord_col='Cell X Position', y_coord_col='Cell Y Position', plot_real_data=True, roi_figsize=(6, 4), marker_size_step=0.80, alpha=1, roi_dpi=200)

      Sample #2 --> should reproduce what's generated in utils.plot_whole_slide_patches(slices, roi_figsize=(15, 10), depatch=True) (not with ROIs outlined):
      tci.plot_just_rois(tci.undo_patching_overlaps_and_decompounding(slices.data)[0], slices.plotting_map, slices.num_colors, slices.webpage_dir, slices.mapping_dict, slices.coord_units_in_microns, slide_id_col='Slide ID', roi_id_col='Slide ID', x_coord_col='Cell X Position', y_coord_col='Cell Y Position', plot_real_data=True, roi_figsize=(15, 10), marker_size_step=0.80, alpha=0.1, edgecolors=None, default_marker_size_fac=1, roi_dpi=200, yaxis_dir=-1, boxes_to_plot0=None)

      Sample #3 --> should be the same thing except with darker and outlined circles:
      tci.plot_just_rois(tci.undo_patching_overlaps_and_decompounding(slices.data)[0], slices.plotting_map, slices.num_colors, slices.webpage_dir, slices.mapping_dict, slices.coord_units_in_microns, slide_id_col='Slide ID', roi_id_col
      ='Slide ID', x_coord_col='Cell X Position', y_coord_col='Cell Y Position', plot_real_data=True, roi_figsize=(15, 10), marker_size_step=0.80, alpha=0.5, edgecolors='k', default_marker_size_fac=1, roi_dpi=200, yaxis_dir=-1, boxes_to_plot0=None)
    '''

    # Antonio's fix to enable plot generation in SLURM's batch mode
    import matplotlib
    matplotlib.use('Agg')

    # Import relevant libraries
    import numpy as np
    import os

    # Define the directory holding all the images for the webpage and the filename of the file holding all the corresponding figure data
    if webpage_dir_addon is None:
        savedir = os.path.join(webpage_dir, ('real' if plot_real_data else 'simulated'))
    else:
        savedir = os.path.join(webpage_dir, webpage_dir_addon)
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    # Extract the correct number of colors from the default color palette
    ielem = 0
    colors = []
    for elem in matplotlib.rcParams['axes.prop_cycle']():
        color = elem['color']
        colors.append(color)
        ielem = ielem + 1
        if ielem == num_colors:
            break
    default_marker_size = matplotlib.rcParams['lines.markersize'] * default_marker_size_fac

    # Get the unique slides in the dataset
    unique_slides = df[slide_id_col].unique()

    # For each slide...
    for islide, slide_name in enumerate(unique_slides):

        print('On slide {} of {}...'.format(islide + 1, len(unique_slides)))

        # Get the boolean indexes of the current slide
        slide_loc = df[slide_id_col] == slide_name

        # Get the unique ROIs in the current slide
        unique_rois = df.loc[slide_loc, roi_id_col].unique()

        # For each ROI in the slide...
        for iroi, roi_name in enumerate(unique_rois):

            print('  On ROI {} of {}...'.format(iroi + 1, len(unique_rois)))

            # Get the boolean indexes of the current ROI of the current slide
            roi_loc = slide_loc & (df[roi_id_col] == roi_name)

            # Collect the data needed to plot the ROIs, plot the ROI, and save it to disk
            x_roi = df.loc[roi_loc, x_coord_col].to_numpy()
            y_roi = df.loc[roi_loc, y_coord_col].to_numpy()
            x_range = np.array([x_roi.min(), x_roi.max()])
            y_range = np.array([y_roi.min(), y_roi.max()])
            species_roi = np.array(df.loc[roi_loc, 'Species int'], dtype='uint64')
            uroi = roi_name
            unique_species_in_roi = np.unique(species_roi)
            spec2plot_roi = [x[0] for x in plotting_map if x[0] in unique_species_in_roi]

            # Plot the ROI and potentially the ROI outlines, which are potentially colored
            if not plot_just_single_roi_outline:
                boxes_to_plot = (boxes_to_plot0[islide] if boxes_to_plot0 is not None else None)
                tag = roi_name
                plot_and_save_roi((roi_figsize, spec2plot_roi, species_roi, x_roi, y_roi, plotting_map, colors, x_range, y_range, uroi, marker_size_step, default_marker_size, roi_dpi, mapping_dict, coord_units_in_microns, alpha, edgecolors, yaxis_dir, boxes_to_plot, pval_params, title_suffix, tag, filename_suffix, savedir))
            else:
                boxes_for_slide = boxes_to_plot0[islide]
                list_of_tuple_arguments = []
                for iroi_box in range(len(boxes_for_slide)):
                    boxes_to_plot = boxes_for_slide.iloc[iroi_box].to_frame().T  # this is a one-row dataframe
                    tag = boxes_to_plot.index[0]
                    list_of_tuple_arguments.append((roi_figsize, spec2plot_roi, species_roi, x_roi, y_roi, plotting_map, colors, x_range, y_range, uroi, marker_size_step, default_marker_size, roi_dpi, mapping_dict, coord_units_in_microns, alpha, edgecolors, yaxis_dir, boxes_to_plot, pval_params, title_suffix, tag, filename_suffix, savedir))
                utils.execute_data_parallelism_potentially(function=plot_and_save_roi, list_of_tuple_arguments=list_of_tuple_arguments, nworkers=(0 if not use_multiprocessing else nworkers), task_description='plotting of the single ROI outlines on whole slide images')


def get_integer_index(val_to_map, vmin=-200, vmax=0, N=2**8):
    """Map a data value to an integer index, useful for indexing custom colormaps.

    Args:
        val_to_map (float): Data value we want to map to an integer index
        vmin (float, optional): Minimum allowed value of the data. Defaults to -200.
        vmax (float, optional): Maximum allowed value of the data. Defaults to 0.
        N (int, optional): Number of possible indexes to which to map val_to_map. Defaults to 2**8.

    Returns:
        int: Integer index in the range 0, 1, 2, ..., N-1
    """

    # Import relevant library
    import numpy as np

    # Clip to within the allowed data range
    val_to_map = max([vmin, val_to_map])
    val_to_map = min([val_to_map, vmax])

    # Return an integer that is closest to the allowed values: 0, 1, ..., N-1
    return int(np.round((val_to_map - vmin) / (vmax - vmin) * (N - 1)))


def get_properly_interpolated_colormap(colormap_name='rocket', N=2**8):
    """Generate a colormap based on a known colormap that can be properly interpolated beyond 256 unique colors. See https://stackoverflow.com/questions/56469836/generate-16-bit-color-map-in-matplotlib. The point of this is that lut doesn't work as expected in e.g. "mpl.cm.get_cmap('viridis', lut=1000)", which yields a colormap of length lut BUT there are only 256 unique values. E.g., see output of this as is basically outputted from get_rgba_for_data_on_colormap() above:

      cmap = matplotlib.cm.get_cmap('viridis', lut=2**16)
      np.unique(cmap(norm(np.linspace(-200, 0, 1000)))[:, 1]).shape

    Args:
        colormap_name (str, optional): Known colormap name. Defaults to 'rocket', which is the default for Seaborn's heatmap() function.
        N (int, optional): Desired number of unique colors generated from the known colormap. Defaults to 2**8.

    Returns:
        cmap: Colormap
    """

    # Import relevant libraries
    import matplotlib as mpl
    import seaborn  # needed to import the 'rocket' colormap

    # Return a cmap, without a lookup table, but which you can use with an integer index, as e.g. outputted from get_integer_index() above
    return mpl.colors.LinearSegmentedColormap.from_list('', mpl.cm.get_cmap(colormap_name).colors, N=N)


def get_properly_interpolated_rgba_from_data(val_to_map, vmin=-200, vmax=0, N=2**8, colormap_name='rocket', print_byte_value=False):
    """Get an RGBA tuple corresponding to an input data value, possible min/max values, a total colormap size, and the name of a known colormap.

    **NOTE**: This replaces get_rgba_for_data_on_colormap() above!!

    Args:
        val_to_map (float): Data value for which we want an RGBA tuple
        vmin (float, optional): Minimum allowed value of the data. Defaults to -200.
        vmax (float, optional): Maximum allowed value of the data. Defaults to 0.
        N (int, optional): Desired number of unique colors generated from a known colormap. Defaults to 2**8.
        colormap_name (str, optional): Known colormap name. Defaults to 'rocket', which is the default for Seaborn's heatmap() function.
        print_byte_value (bool, optional): Whether to print a 0-255-ranged RGBA value for testing. Defaults to False.

    Returns:
        tuple: RGBA tuple corresponding to the input data
    """

    # Import relevant library
    import numpy as np

    # Get the appropriately normalized integer index corresponding to a datapoint, min/max values, and total number of indexes
    custom_colormap_index = get_integer_index(val_to_map, vmin=vmin, vmax=vmax, N=N)

    # Get a custom colormap that is interpolated possibly beyond 256 unique colors based on a known colormap
    cmap = get_properly_interpolated_colormap(colormap_name=colormap_name, N=N)

    # Get the RGBA values
    rgba_tuple = cmap(custom_colormap_index)

    # If we want to print out the RGBA values for testing, do so (0-255)
    if print_byte_value:
        print((np.array(rgba_tuple) * 255).round())

    # Return the corresponding RGBA tuple
    return rgba_tuple


def plot_single_density_pvals(argument_tuple):
    """Plot a single set of left and right log density P values in a parameter format convenient for the original use of the multiprocessing module

    Args:
        argument_tuple (tuple): tuple containing the function parameters
    """

    # Get the arguments as a tuple, for original paramater format for the multiprocessing module
    log_pval_range, figsize, webpage_dir_for_dens_pvals_per_entity, plot_real_data, img_file_suffix, all_species_names, dpi, df_density_pvals_arrays, title_suffix, entity_index = argument_tuple

    # Define the ROI-specific variables from the main df_density_pvals_arrays parameter
    log_dens_pvals_arr = df_density_pvals_arrays.loc[entity_index, 'log_dens_pvals_arr']
    entity_name = df_density_pvals_arrays.loc[entity_index, 'roi_name']
    entity = 'roi'

    # Plot the heatmaps for the current set of data
    plot_density_pvals_simple(log_dens_pvals_arr, log_pval_range, figsize, dpi, webpage_dir_for_dens_pvals_per_entity, plot_real_data, entity_name, img_file_suffix, entity, entity_index, all_species_names, title_suffix=title_suffix)

    # Save the figure pathname; this may not work, in which case delete this line as it doesn't matter much
    # Didn't work, commenting out for now
    # df_data_by_roi.loc[entity_index, 'density_pval_fig_pathname'] = filename


def calculate_recommended_radius(contents, recommended_radius_frac_of_side_length=0.1):
    """Print a reasonable radius around each center (for the counting of the neighbors) based on the average side length of each ROI

    Args:
        contents (Pandas dataframe): "Contents" dataframe containing the width and height of every ROI (these columns must be present)
        recommended_radius_frac_of_side_length (float, optional): Recommended fraction of the side length for the radius. Defaults to 0.1.
    """
    import numpy as np
    # recommended_radius = (contents[['width', 'height']].mean().min() * recommended_radius_frac_of_side_length).round().astype(int)
    recommended_radius = np.round(contents[['width', 'height']].mean().min() * recommended_radius_frac_of_side_length).astype(int)
    print('Recommended "thickness" parameter in microns, e.g., radius in microns around each center for the counting of neighbors: {}'.format(recommended_radius))
    return(recommended_radius)


def plot_single_roi(args_as_single_tuple):
    """Plot and save the ROI plot for a single ROI; this way we can parallelize over all ROIs.

    Args:
        args_as_single_tuple (tuple): tuple containing all necessary data (see "Unpack the arguments" block below)
    """

    # Import relevant modules
    import matplotlib.pyplot as plt
    import os

    # Unpack the arguments
    plotting_map, roi_figsize, colors, marker_size_step, default_marker_size, roi_dpi, mapping_dict, coord_units_in_microns, alpha, edgecolors, yaxis_dir, savedir, df_data_by_roi, roi_index = args_as_single_tuple

    # Assign the data for the current ROI to useful variables
    roi_name = df_data_by_roi.loc[roi_index, 'unique_roi']  # df_data_by_roi.loc[roi_index, 'roi_name']
    roi_fig_filename = 'roi_plot_{}_{}.png'.format(roi_name, roi_index)
    roi_fig_pathname = os.path.join(savedir, roi_fig_filename)

    # If the image file doesn't already exist...
    # Well, actually I believe existence is basically checked for in the calling function plot_rois(), so this check should not be necessary!
    if not os.path.exists(roi_fig_pathname):

        # Note what we're doing
        print('Plotting and saving ROI {} (ROI index {})'.format(roi_name, roi_index))

        # Collect the data needed to plot the ROIs, plot the ROI, and save it to disk
        x_range = df_data_by_roi.loc[roi_index, 'x_range']
        y_range = df_data_by_roi.loc[roi_index, 'y_range']
        x_roi = df_data_by_roi.loc[roi_index, 'x_roi']
        y_roi = df_data_by_roi.loc[roi_index, 'y_roi']
        species_roi = df_data_by_roi.loc[roi_index, 'species_roi']
        spec2plot_roi = df_data_by_roi.loc[roi_index, 'spec2plot_roi']

        # Define the figure to use for plotting the ROI
        fig_roi = plt.subplots(figsize=roi_figsize)[0]

        # Plot and save the current ROI
        plot_roi(fig_roi, spec2plot_roi, species_roi, x_roi, y_roi, plotting_map, colors, x_range, y_range, roi_name, marker_size_step, default_marker_size, roi_dpi, mapping_dict, coord_units_in_microns, filepath=None, do_plot=True, alpha=alpha, edgecolors=edgecolors, yaxis_dir=yaxis_dir)
        fig_roi.savefig(roi_fig_pathname, dpi=roi_dpi, bbox_inches='tight')

        # Close the figure for the ROI
        plt.close(fig_roi)

        # Save the path to the ROI plot
        # I don't usually write to a common object that's being accessed in parallel, so if things get weird or die, try removing this line, or try initializing it in plot_rois()
        # It actually seems to not work at all (7/5/22); no such column is added; trying to add column instantiation for analogous metrics plotting functionality, but if that doesn't work either then can probably scrap both new column generations
        # 10-4-23: Just creating the column first in plot_rois()
        df_data_by_roi.loc[roi_index, 'roi_plot_fig_pathname'] = roi_fig_pathname

    # If the image file DOES already exist...
    else:

        # The image file already exists
        print('The image file {} already exists'.format(roi_fig_pathname))


def generate_dens_pvals_array_for_roi(args_as_single_tuple):
    """Perform some validity checks on the flattened metrics data.

    This function keeps in df_density_pvals_arrays all the same exact ROIs as in df_density_pvals except puts the data into numpy arrays, except when for a given center-neighbor pair there is fewer than num_valid_centers_minimum centers.

    Args:
        args_as_single_tuple (tuple): All arguments in tuple form
    """

    # Import relevant libraries
    import numpy as np

    # Obtain the parameters for the current ROI
    pickle_dir, nall_species, all_species_ids, nslices, num_valid_centers_minimum, log_pval_range, debug, df_data_by_roi, df_density_pvals, correct_flooring, iroi = args_as_single_tuple

    # Determine the pickle file to be created for the current ROI that will hold the ROI name and the final numpy arrays to use
    pickle_file = 'dens_pvals_array-roi_index_{:06}.pkl'.format(iroi)

    # Initialize the arrays of interest
    log_dens_pvals_arr = np.zeros((nall_species, nall_species, 2, nslices))
    log_dens_pvals_arr[:] = None
    num_valid_centers = np.zeros((nall_species, nall_species, nslices))  # this doesn't need a neighbor index but we're doing it for parallelism with log_dens_pvals_arr
    num_valid_centers[:] = None
    roi_center_neighbor_holder = np.ones((nall_species, nall_species, 2), dtype=np.int64) * -1

    # Use the "contents" dataframe to determine the current ROI name
    roi_name = df_data_by_roi.loc[iroi, 'unique_roi']

    # Get the slice of the density P values dataframe that corresponds to just the current ROI
    df_density_pvals_roi = df_density_pvals[df_density_pvals['roi_name'] == roi_name]

    # Keep only those entries in the data for which there is a minimum number of valid centers
    df_valid = df_density_pvals_roi[df_density_pvals_roi['nvalid_centers_per_slice'].apply(lambda x: x[0]) >= num_valid_centers_minimum]

    # For every possible center and neighbor species...
    for icenter_spec, center_spec in enumerate(all_species_ids[:]):
        for ineighbor_spec, neighbor_spec in enumerate(all_species_ids[:]):

            # Get a dataframe of, if present, the P value data for the current ROI/center/neighbor. There should never be more than one match, and there could be zero matches due to having too few center species present
            one_or_zero_rows = df_valid[(df_valid['center_species_id'] == center_spec) & (df_valid['neighbor_species_id'] == neighbor_spec)]
            nrows = len(one_or_zero_rows)

            # Ensure there is exactly zero or one matches in the data for the current ROI/center/neighbor
            if nrows not in {0, 1}:
                print('ERROR: There is more than one match for the current ROI/center/neighbor: {}/{}/{}'.format(roi_name, center_spec, neighbor_spec))
                exit()

            # For debugging purposes, if there are no matches for the current ROI/center/neighbor, print that information
            if nrows == 0:
                if debug:
                    print('NOTE: There is no data for ROI/center/neighbor: {}/{}/{}'.format(roi_name, center_spec, neighbor_spec))

            # If there's a single match for the current ROI/center/neighbor (nrows=1)
            else:

                # Store the current index of the single row since we use it twice below
                index = one_or_zero_rows.index[0]

                # Store the log of the left and right P values in the numpy array
                log_dens_pvals_arr[icenter_spec, ineighbor_spec, :, 0] = np.array([one_or_zero_rows.loc[index, 'left_log_dens_pvals'][0][0],
                                                                                   one_or_zero_rows.loc[index, 'right_log_dens_pvals'][0][0]])

                # I SHOULD PROBABLY MOVE THIS OUT OF THE ELSE STATEMENT!!!!, I doubt that would hurt anything and would rather be generally helpful!
                roi_center_neighbor_holder[icenter_spec, ineighbor_spec, :] = np.array([center_spec, neighbor_spec])

                # Store the current number of valid centers
                num_valid_centers[icenter_spec, ineighbor_spec, 0] = one_or_zero_rows.loc[index, 'nvalid_centers_per_slice'][0]

    # Set the (negative) infinite values to the darkest color (or else they won't be plotted, as inf values are not plotted)
    if not correct_flooring:
        log_dens_pvals_arr[np.isneginf(log_dens_pvals_arr)] = log_pval_range[0]  # this old way produces -50 a log but also -100, -200, etc.. The per-ROI heatmaps are still plotted correctly since their color gets floored to the minimum regardless, but e.g. the averaging is not done fairly when e.g. averaging over the ROIs for each slide.
    else:  # new, correct way
        # Since those values converted-to--50 etc. values are then greater than some of values that aren't negative infinity but are small (e.g., -100), also convert those values to the minimum possible value
        log_dens_pvals_arr[log_dens_pvals_arr < log_pval_range[0]] = log_pval_range[0]  # should produce -50 at minimum

    # Check there are no log P values that are greater than zero (check that the P value is never greater than 1)
    assert (log_dens_pvals_arr > log_pval_range[1]).sum() == 0

    # Create the pickle file holding the final numpy arrays
    make_pickle({'roi_name': roi_name, 'log_dens_pvals_arr': log_dens_pvals_arr, 'num_valid_centers': num_valid_centers, 'centers_neighbors_arr': roi_center_neighbor_holder}, pickle_dir, pickle_file)


def plot_density_pvals_simple(log_dens_pvals_arr, log_pval_range, figsize, dpi, plots_dir, plot_real_data, entity_name, img_file_suffix, entity, entity_index, all_species_names, title_suffix=''):

    # Import relevant libraries
    import matplotlib.pyplot as plt
    import seaborn as sns
    import os

    # Determine the number of slices from the main array to plot
    nslices = log_dens_pvals_arr.shape[3]

    # Set the log of the P value range for the color plotting
    vmin = log_pval_range[0]
    vmax = log_pval_range[1]

    # Initialize the figure and axes
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=figsize)

    # Plot the density P values for every slice
    for islice in range(nslices):

        # Determine the filename of the figure
        filename = os.path.join(plots_dir, 'density_pvals-{}-{}-slice_{:02d}_of_{:02d}{}-{}_index_{}.png'.format(('real' if plot_real_data else 'simulated'), entity_name, islice + 1, nslices, img_file_suffix, entity, entity_index))

        # Reset the figure/axes
        fig.clf()
        ax = fig.subplots(nrows=1, ncols=2)

        # Plot the log of the left/less P values
        sns.heatmap(log_dens_pvals_arr[:, :, 0, islice], vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[0], cbar=True, xticklabels=all_species_names, yticklabels=all_species_names, square=True)
        ax[0].set_xticklabels(all_species_names, rotation=45)
        ax[0].set_yticklabels(all_species_names, rotation=30)
        ax[0].set_title('log10(\"less\" density pvals)')
        ax[0].set_xlabel('Neighbor species')
        ax[0].set_ylabel('Center species')

        # Plot the log of the right/greater P values
        sns.heatmap(log_dens_pvals_arr[:, :, 1, islice], vmin=vmin, vmax=vmax, linewidths=.5, ax=ax[1], cbar=True, xticklabels=all_species_names, yticklabels=all_species_names, square=True, cbar_kws={'alpha': 0})  # note the colorbar transparency setting here does not seem to have any effect
        ax[1].set_xticklabels(all_species_names, rotation=45)
        ax[1].set_yticklabels(all_species_names, alpha=0, rotation=30)  # alpha=0 is a hack to make the y-axis label not appear yet still have the heatmap be the size it would otherwise be
        ax[1].set_title('log10(\"greater\" density pvals)')
        ax[1].set_xlabel('Neighbor species')
        ax[1].set_ylabel('Center species', alpha=0)  # alpha=0 is a hack to make the y-axis label not appear yet still have the heatmap be the size it would otherwise be

        # Set the figure title to the slide title and ensure the facecolor is white
        # fig.suptitle('Density P values - {} - {} data - slice {} of {}'.format(entity_name, ('real' if plot_real_data else 'simulated'), islice + 1, nslices))
        fig.suptitle('Density P values - {}{}'.format(entity_name, title_suffix))
        fig.patch.set_facecolor('white')

        # Save the figure to disk
        fig.savefig(filename, dpi=dpi, bbox_inches='tight')

    # Close the figure
    plt.close(fig)

def plot_pvals_over_slides_for_center_neighbor_pair(args_as_single_tuple):

    import numpy as np

    rectangle_data_orig, unique_slides, df_density_pvals_arrays, filename_suffix2, depatched, plotting_map, num_colors, webpage_dir, mapping_dict, coord_units_in_microns, roi_figsize, marker_size_step, alpha, edgecolors, default_marker_size_fac, roi_dpi, yaxis_dir, vmin, vmax, icenter, ineighbor, center_name, neighbor_name, webpage_dir_addon = args_as_single_tuple

    # Print what we're plotting
    print('Plotting P values of {} around {} over all slides'.format(neighbor_name, center_name))

    # Get a fresh copy of the rectangle data for the current center/neighbor pair
    rectangle_data = rectangle_data_orig.copy()

    # For each slide of rectangle data...
    for curr_rectangle_data, unique_slide in zip(rectangle_data, unique_slides):
        print('  On slide {}'.format(unique_slide))

        # For each ROI of non-zero area...
        for roi_name in curr_rectangle_data.index:

            # Get the row (or lack thereof) of the plotting P value dataframe corresponding to the current ROI and get the number of matches
            one_or_zero_rows = df_density_pvals_arrays[df_density_pvals_arrays['roi_name'] == roi_name]
            num_roi_matches = len(one_or_zero_rows)

            # Initialize the P values to be stored in the current rectangle data dataframe
            left_right_for_rectangle_data = [None, None]

            # If a ROI corresponding to the current ROI name exists and has P value data for *some* center-neighbor pair...
            if num_roi_matches == 1:

                # Print what we're doing
                print('ROI {} does indeed have P value data for *some* center-neighbor pair'.format(roi_name))

                # Get a two-element array of the left and right log density P value
                left_right = df_density_pvals_arrays.loc[one_or_zero_rows.index[0], 'log_dens_pvals_arr'][icenter, ineighbor, :, 0].flatten()  # (nall_species, nall_species, 2, nslices)

                # If there is heatmap data for the current ROI but not for the current center-neighbor pair, don't assign P values to the rectangles dataframe
                if np.isnan(left_right[0]):
                    print('No data are available for ROI {} for center {} and neighbor {}'.format(roi_name, center_name, neighbor_name))

                # Save the P values for the current ROI into the dataframes in the rectangles list
                else:
                    print('Data found for ROI {}, center {}, and neighbor {}'.format(roi_name, center_name, neighbor_name))
                    left_right_for_rectangle_data = left_right

            # If there is no match, set the P values for the current ROI to None so they will not be plotted
            else:

                # Print informational messages if there's not a single match for the current ROI
                if num_roi_matches == 0:
                    print('No ROIs with any heatmap data at all have been found for ROI {}'.format(roi_name))
                elif num_roi_matches > 1:
                    print('More than one ROI appears to match for a single ROI name; examine what\'s going on for ROI {}'.format(roi_name))

            # Assign the left and right P values to the rectangles dataframe
            curr_rectangle_data.loc[roi_name, 'left_log_dens_pval'] = left_right_for_rectangle_data[0]
            curr_rectangle_data.loc[roi_name, 'right_log_dens_pval'] = left_right_for_rectangle_data[1]

    # Do plotting (twice, one for left and one for right) for the current center and neighbor for all slides
    # Left log density P values
    left_or_right_string = 'left'
    title_suffix = ' - dispersion of {} around {}'.format(neighbor_name, center_name)
    filename_suffix = '-with_log_dens_pvals_per_roi__center_{}__neighbor_{}__{}_pvals{}'.format(center_name.strip(), neighbor_name.strip(), left_or_right_string, filename_suffix2)
    pval_colname = '{}_log_dens_pval'.format(left_or_right_string)
    plot_just_rois(depatched, plotting_map, num_colors, webpage_dir, mapping_dict, coord_units_in_microns, slide_id_col='Slide ID', roi_id_col='Slide ID', x_coord_col='Cell X Position', y_coord_col='Cell Y Position', plot_real_data=True, roi_figsize=roi_figsize, marker_size_step=marker_size_step, alpha=alpha, edgecolors=edgecolors, default_marker_size_fac=default_marker_size_fac, roi_dpi=roi_dpi, yaxis_dir=yaxis_dir, boxes_to_plot0=rectangle_data, filename_suffix=filename_suffix, pval_params=(pval_colname, vmin, vmax, False), title_suffix=title_suffix, webpage_dir_addon=webpage_dir_addon)

    # Right log density P values
    left_or_right_string = 'right'
    title_suffix = ' - aggregation of {} around {}'.format(neighbor_name, center_name)
    filename_suffix = '-with_log_dens_pvals_per_roi__center_{}__neighbor_{}__{}_pvals{}'.format(center_name.strip(), neighbor_name.strip(), left_or_right_string, filename_suffix2)
    pval_colname = '{}_log_dens_pval'.format(left_or_right_string)
    plot_just_rois(depatched, plotting_map, num_colors, webpage_dir, mapping_dict, coord_units_in_microns, slide_id_col='Slide ID', roi_id_col='Slide ID', x_coord_col='Cell X Position', y_coord_col='Cell Y Position', plot_real_data=True, roi_figsize=roi_figsize, marker_size_step=marker_size_step, alpha=alpha, edgecolors=edgecolors, default_marker_size_fac=default_marker_size_fac, roi_dpi=roi_dpi, yaxis_dir=yaxis_dir, boxes_to_plot0=rectangle_data, filename_suffix=filename_suffix, pval_params=(pval_colname, vmin, vmax, False), title_suffix=title_suffix, webpage_dir_addon=webpage_dir_addon)

def plot_and_save_roi(args_as_single_tuple):
    import matplotlib.pyplot as plt
    import os
    roi_figsize, spec2plot_roi, species_roi, x_roi, y_roi, plotting_map, colors, x_range, y_range, uroi, marker_size_step, default_marker_size, roi_dpi, mapping_dict, coord_units_in_microns, alpha, edgecolors, yaxis_dir, boxes_to_plot, pval_params, title_suffix, tag, filename_suffix, savedir = args_as_single_tuple
    fig_roi = plt.subplots(figsize=roi_figsize)[0]  # define the figures to use for plotting the ROIs and the P values
    plot_roi(fig_roi, spec2plot_roi, species_roi, x_roi, y_roi, plotting_map, colors, x_range, y_range, uroi, marker_size_step, default_marker_size, roi_dpi, mapping_dict, coord_units_in_microns, alpha=alpha, edgecolors=edgecolors, yaxis_dir=yaxis_dir, boxes_to_plot=boxes_to_plot, pval_params=pval_params, title_suffix=title_suffix)
    roi_fig_filename = '{}{}.png'.format(tag, filename_suffix)
    roi_fig_pathname = os.path.join(savedir, roi_fig_filename)
    fig_roi.savefig(roi_fig_pathname.replace(' ', '_'), dpi=roi_dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig_roi)  # close the figure
