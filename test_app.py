import streamlit as st
import os
from streamlit_extras.app_logo import add_logo
import streamlit_session_state_management
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import streamlit_utils
import install_missing_packages
import sandbox
import subprocess

install_missing_packages.live_package_installation()

# Note if any of the following imports having "  # slow" are not commented out, there is a delay in running the forking test
from pages2 import data_import_and_export
from pages2 import datafile_format_unifier
from pages2 import open_file
from pages2 import robust_scatter_plotter
from pages2 import multiaxial_gating
from pages2 import thresholded_phenotyping  # slow due to things ultimately importing umap
from pages2 import adaptive_phenotyping
from pages2 import Pheno_Cluster_a  # "slow" for forking test initialization
from pages2 import Pheno_Cluster_b  # "slow" for forking test initialization
from pages2 import Tool_parameter_selection
from pages2 import Run_workflow
from pages2 import Display_individual_ROI_heatmaps
from pages2 import Display_average_heatmaps
from pages2 import Display_average_heatmaps_per_annotation
from pages2 import Display_ROI_P_values_overlaid_on_slides
from pages2 import Neighborhood_Profiles  # slow due to things ultimately importing umap
from pages2 import UMAP_Analyzer  # slow due to things ultimately importing umap
from pages2 import Clusters_Analyzer  # slow due to things ultimately importing umap
from pages2 import memory_analyzer
from pages2 import radial_bins_plots
from pages2 import radial_profiles_analysis
from pages2 import preprocessing
from pages2 import results_transfer
# from pages2 import forking_test


def welcome_page():
    # Markdown text
    intro_markdown = ndl.read_markdown_file('markdown/MAWA_WelcomePage.md')
    st.markdown(intro_markdown, unsafe_allow_html=True)


def main():

    st.set_page_config(layout="wide")

    # Use the new st.naviation()/st.Page() API to create a multi-page app
    pg = st.navigation({
        'Home 🏠':
            [
                st.Page(welcome_page, title="Welcome", url_path='home')
            ],
        'File Handling 🗄️':
            [
                st.Page(data_import_and_export.main, title="Data Import and Export", url_path='data_import_and_export'),
                st.Page(datafile_format_unifier.main, title="Datafile Unification", url_path='datafile_unification'),
                st.Page(open_file.main, title="Open File", url_path='open_file')
            ],
        'Coordinate Scatter Plotter 🌟':
            [
                st.Page(robust_scatter_plotter.main, title="Coordinate Scatter Plotter", url_path='coordinate_scatter_plotter')
            ],
        'Phenotyping 🧬':
            [
                st.Page(multiaxial_gating.main, title="Using Raw Intensities", url_path='using_raw_intensities'),
                st.Page(thresholded_phenotyping.main, title="Using Thresholded Intensities", url_path='using_thresholded_intensities'),
                st.Page(adaptive_phenotyping.main, title="Adaptive Phenotyping", url_path='adaptive_phenotyping')
            ],
         'Phenotype Clustering Workflow ✨':
             [
                 st.Page(Pheno_Cluster_a.main, title="Unsupervised Phenotype Clustering", url_path='unsupervised_phenotype_clustering'),
                 st.Page(Pheno_Cluster_b.main, title="Differential Intensity", url_path='differential_intensity')
             ],
        'Spatial Interaction Tool 🗺️':
            [
                st.Page(Tool_parameter_selection.main, title="Tool Parameter Selection", url_path='tool_parameter_selection'),
                st.Page(Run_workflow.main, title="Run SIT Workflow", url_path='run_sit_workflow'),
                st.Page(Display_individual_ROI_heatmaps.main, title="Display Individual ROI Heatmaps", url_path='display_individual_roi_heatmaps'),
                st.Page(Display_average_heatmaps.main, title="Display Average Heatmaps", url_path='display_average_heatmaps'),
                st.Page(Display_average_heatmaps_per_annotation.main, title="Display Average Heatmaps per Annotation", url_path='display_average_heatmaps_per_annotation'),
                st.Page(Display_ROI_P_values_overlaid_on_slides.main, title="Display ROI P Values Overlaid on Slides", url_path='display_roi_p_values_overlaid_on_slides')
            ],
        'Neighborhood Profiles Workflow 🌳':
            [
                st.Page(Neighborhood_Profiles.main, title="Neighborhood Profiles", url_path='neighborhood_profiles'),
                st.Page(UMAP_Analyzer.main, title="UMAP Differences Analyzer", url_path='umap_differences_analyzer'),
                st.Page(Clusters_Analyzer.main, title="Clusters Analyzer", url_path='clusters_analyzer')
            ],
        'Radial Profiles 🌀':
            [
                st.Page(radial_bins_plots.main, title="Radial Bins Plots", url_path='radial_bins_plots'),
                st.Page(radial_profiles_analysis.main, title="Radial Profiles Analysis", url_path='radial_profiles_analysis')
            ],
        'Utilities 🛠️':
            [
                st.Page(preprocessing.main, title="Preprocessing", url_path='preprocessing'),
                st.Page(memory_analyzer.main, title="Memory Analyzer", url_path='memory_analyzer'),
                st.Page(results_transfer.main, title="Results Transfer", url_path='results_transfer'),
                # st.Page(forking_test.main, title="Forking Test", url_path='forking_test')
            ]
        })

    # Set the top directories, which may or may not be the actual input, output, and saved states directories for the session (instead, they could be subdirectories in those top directories)
    input_top_dir = os.path.join('.', 'input')
    output_top_dir = os.path.join('.', 'output')
    saved_states_top_dir = os.path.join('.', 'saved_streamlit_session_states')

    # Creates session state keys input_user_dir, output_user_dir, saved_states_user_dir. These should be the directories used in MAWA for filesystem interaction
    if 'input_user_dir' not in st.session_state:
        sandbox.set_up_user_directories(input_top_dir, output_top_dir, saved_states_top_dir)

    # Creates session state key num_cpus_for_sit
    if 'num_cpus_for_sit' not in st.session_state:
        st.session_state['num_cpus_for_sit'] = sandbox.get_num_cpus_for_sit()

    print(f"Input user directory: {st.session_state['input_user_dir']}")
    print(f"Output user directory: {st.session_state['output_user_dir']}")
    print(f"Saved states user directory: {st.session_state['saved_states_user_dir']}")
    print(f"Number of CPUs for SIT: {st.session_state['num_cpus_for_sit']}")
    print(subprocess.run(['uname', '-a'], capture_output=True, text=True, check=True).stdout.strip())

    # # Allow the user to delete their user directories, just uncomment if so
    # if st.button('Delete user directories'):
    #     delete_user_directories()

    # For widget persistence, we need always copy the session state to itself, being careful with widgets that cannot be persisted, like st.data_editor() (where we use the "__do_not_persist" suffix to avoid persisting it)
    for key in st.session_state.keys():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')):
            st.session_state[key] = st.session_state[key]

    # This is needed for the st.dataframe_editor() class (https://github.com/andrew-weisman/streamlit-dataframe-editor) but is also useful for seeing where we are and where we've been
    st.session_state['current_page_name'] = pg.url_path if pg.url_path != '' else 'Home'
    if 'previous_page_name' not in st.session_state:
        st.session_state['previous_page_name'] = st.session_state['current_page_name']

    # Add logo to sidebar
    add_logo('app_images/mawa_logo-width315.png', height=250)

    # Determine whether this is the first time the app has been run
    if 'app_has_been_run_at_least_once' not in st.session_state:
        st.session_state['app_has_been_run_at_least_once'] = True
        first_app_run = True
    else:
        first_app_run = False

    # Run session state management in the sidebar
    streamlit_session_state_management.execute(first_app_run)

    # Initalize session_state values for streamlit processing
    if 'init' not in st.session_state:
        st.session_state = ndl.init_session_state(st.session_state)

    # Sidebar organization
    with st.sidebar:
        st.write('**:book: [Documentation](https://ncats.github.io/multiplex-analysis-web-apps/)**')
        with st.expander('Advanced:'):
            benchmark_button = True
            if benchmark_button:
                st.button('Record Benchmarking', on_click = st.session_state.bc.save_run_to_csv)
            if st.button('Calculate memory used by Python session'):
                streamlit_utils.write_python_session_memory_usage()

    # Check the platform
    st.session_state = sandbox.check_for_platform(st.session_state)

    # Format tooltips
    tooltip_style = """
        <style>
        div[data-baseweb="tooltip"] {
        width: 250px;
        }
        </style>
    """
    st.markdown(tooltip_style,unsafe_allow_html=True)

    # On every page, display its title
    st.title(pg.title)

    # Render the select page
    pg.run()

    # Update the previous page location
    st.session_state['previous_page_name'] = st.session_state['current_page_name']


# Needed for rendering pages which use multiprocessing (https://docs.python.org/3/library/multiprocessing.html#the-spawn-and-forkserver-start-methods)
if __name__ == '__main__':
    main()
