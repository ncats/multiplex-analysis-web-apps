# Import relevant libraries.
#### Framework imports.
import streamlit as st
import generate_results
import framework.startup as startup
import framework.manage_sessions as manage_sessions
import framework.inspect_database_tables as inspect_database_tables
import framework.monitor_jobs as monitor_jobs
import framework.platform_abstraction as pa
#### App imports.
import re
import streamlit_utils
from pages2 import data_import_and_export
from pages2 import datafile_format_unifier
from pages2 import open_file
from pages2 import feature_creation
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
from streamlit_extras.app_logo import add_logo
import streamlit_session_state_management
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
from fast_neighborhood_profiles import load_unified_input_file
from fast_neighborhood_profiles import phenotype
from fast_neighborhood_profiles import run_spatial_umap
from fast_neighborhood_profiles import assign_neighborhood_types
from fast_neighborhood_profiles import sample_analysis_1
from fast_neighborhood_profiles import sample_analysis_2
from pathlib import Path
from fast_neighborhood_profiles import plot_neighborhood_types

ST_KEY_PREFIX = "app.py__"
ST_KEY_PREFIX_STARTUP = "startup.py__"


def welcome_page():
    '''
    First page displayed when the app opens
    '''
    # Markdown text
    with open("markdown/MAWA_WelcomePage.md", "r", encoding="utf-8") as f:
        md_content = f.read()

    parts = re.split(r"!\[(.*?)\]\((.*?)\)", md_content)
    for i, part in enumerate(parts):
        if i % 3 == 0:
            st.markdown(part, unsafe_allow_html=True)
        elif i % 3 == 1:
            title = part
        else:
            st.image(part)


def suggestions_page():
    md_path = Path("./fast_neighborhood_profiles/suggestions.md")
    md_text = md_path.read_text(encoding="utf-8")
    st.markdown(md_text)


# Define the main function.
def main():

    # Run one-time initialization.
    key = ST_KEY_PREFIX + "app_initialized"
    if key not in st.session_state:
        startup.initialize()
        st.session_state[key] = True
        first_app_run = True
    else:
        first_app_run = False

    # Define the pages for the navigation bar.
    pg = st.navigation(
        {
            "High-performance workflow": [
                st.Page(suggestions_page, title="Start Here!", url_path='suggestions', default=True),
                st.Page(manage_sessions.main, title="Manage sessions", url_path='manage_sessions'),
                st.Page(load_unified_input_file.main, title="Load unified input file", url_path='load_unified_input_file'),
                st.Page(phenotype.main, title="Phenotype", url_path='phenotype'),
                st.Page(run_spatial_umap.main, title="Run spatial UMAP", url_path='run_spatial_umap'),
                st.Page(assign_neighborhood_types.main, title="Assign neighborhood types", url_path='assign_neighborhood_types'),
                st.Page(plot_neighborhood_types.main, title="Plot neighborhood types", url_path='plot_neighborhood_types'),
                st.Page(monitor_jobs.main, title="Monitor jobs", url_path='monitor_jobs'),
                ],
            'Home 🏠':
                [
                    st.Page(sample_analysis_1.main, title="Sample analysis 1", url_path='sample_analysis_1'),
                    st.Page(sample_analysis_2.main, title="Sample analysis 2", url_path='sample_analysis_2'),
                    st.Page(generate_results.main, title="Generate results", url_path='generate_results'),
                    st.Page(welcome_page, title="Welcome", url_path='home')
                ],
            'File Handling 🗄️':
                [
                    st.Page(data_import_and_export.main, title="Data Import and Export", url_path='data_import_and_export'),
                    st.Page(datafile_format_unifier.main, title="Datafile Unification", url_path='datafile_unification'),
                    st.Page(open_file.main, title="Open File", url_path='open_file')
                ],
            'Dataset Investigation 🌟':
                [
                    st.Page(feature_creation.main, title="Feature Creation", url_path='feature_creation'),
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
                    st.Page(inspect_database_tables.main, title="Inspect database tables", url_path='inspect_database_tables'),
                    # st.Page(forking_test.main, title="Forking Test", url_path='forking_test')
                ],
        }
    )

    # For widget persistence between pages, we need always copy the session state to itself.
    for key in st.session_state:
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')):  # Could add things like: "(not key.endswith('_button'))".
            st.session_state[key] = st.session_state[key]

    # This is needed for the st.dataframe_editor() class (https://github.com/andrew-weisman/streamlit-dataframe-editor) but is also useful for seeing where we are and where we've been.
    st.session_state['current_page_name'] = pg.url_path if pg.url_path != '' else 'Home'
    if 'previous_page_name' not in st.session_state:
        st.session_state['previous_page_name'] = st.session_state['current_page_name']

    # Add logo to sidebar
    add_logo('app_images/mawa_logo-width315.png', height=250)

    # Run session state management in the sidebar
    streamlit_session_state_management.execute(first_app_run)

    # Dante's session state initialization to initalize session_state values for streamlit processing. Not putting in startup.py so it gets rerun properly when the user hits the "Reset app" button without my having to trace through everything it does and add manually to manage_sessions.reset_session_state().
    if 'init' not in st.session_state:
        st.session_state = ndl.init_session_state(st.session_state)

    # Sidebar organization
    with st.sidebar:

        # App-specific things.
        st.write('**:open_book: [Documentation](https://ncats.github.io/multiplex-analysis-web-apps/)**')
        with st.expander('Advanced:'):
            benchmark_button = True
            if benchmark_button:
                st.button('Record Benchmarking', on_click = st.session_state.bc.save_run_to_csv)
            if st.button('Calculate memory used by Python session'):
                streamlit_utils.write_python_session_memory_usage()

        # Allow user to shut down entire app cleanly.
        with st.container(horizontal=True):
            st.button("🔄 Refresh page", help="If you want to refresh the page, press this button, *not* your browser's refresh button.")
            st.button("🧹 Reset app", help="Reset the app to its initial state.", on_click=manage_sessions.reset_session_state)
            if st.button("🛑 Shut down app", help="Always save the app session prior to shutdown (unless you don't want to resume your work). Even if you have a running job, you can still shut down the app; just make sure you've saved the app session first so you can pick back up where you left off and load the completed job results as usual."):
                pa.record_explicit_shutdown_time(st.session_state[ST_KEY_PREFIX_STARTUP + "app_session_id"])
                pa.shut_down_app()

    # Write a banner if this app session is waiting on a job to complete.
    if "JOB_PENDING" in st.session_state:
        st.warning(f"This app session is awaiting the results of job {st.session_state['JOB_PENDING']['job_id']}. Please monitor this on the job monitor page.")

    # Display the title of the page.
    st.title(pg.title)

    # Display the page.
    pg.run()

    # Update the previous page location.
    st.session_state['previous_page_name'] = st.session_state['current_page_name']


# Run the main function.
if __name__ == "__main__":
    main()
