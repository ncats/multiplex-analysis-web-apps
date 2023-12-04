# Import relevant libraries
import pandas as pd
import seaborn as sns
import streamlit as st
import plotly.graph_objects as go

# Update the dependencies of the selectbox for the current analysis column
def update_dependencies_of_column_to_plot():
    df = st.session_state['mg__df']
    column_to_plot = st.session_state['mg__column_to_plot']
    st.session_state['mg__curr_column_range'] = [df[column_to_plot].min(), df[column_to_plot].max()]
    st.session_state['mg__selected_min_val'] = st.session_state['mg__curr_column_range'][0]
    st.session_state['mg__selected_max_val'] = st.session_state['mg__curr_column_range'][1]

# Set the input datafile path
input_datafile_path = './input/measurementsEpCAMLy51MHCII-exported.csv'

# Initialize some things in the session state
if 'mg__run_before' not in st.session_state:
    st.session_state['mg__run_before'] = False
if 'mg__dfs_to_plot' not in st.session_state:
    st.session_state['mg__dfs_to_plot'] = dict()

# Load the data only once, otherwise keep it in memory
if 'mg__df' not in st.session_state:
    st.session_state['mg__df'] = pd.read_csv(input_datafile_path).select_dtypes(include='number')
df = st.session_state['mg__df']

# Have a dropdown for the column on which to perform a kernel density estimate
st.selectbox(label='Column to plot:', options=df.columns, key='mg__column_to_plot', on_change=update_dependencies_of_column_to_plot)
column_to_plot = st.session_state['mg__column_to_plot']

# Set some session state keys required below
if not st.session_state['mg__run_before']:
    update_dependencies_of_column_to_plot()

# Output information on the column range
column_range = st.session_state['mg__curr_column_range']
st.write('Column range: {}'.format(column_range))

# Determine the x-y data to plot for the selected column, calculating the KDE for each column only once ever
if column_to_plot not in list(st.session_state['mg__dfs_to_plot']):
    line2d = sns.kdeplot(df, x=column_to_plot).get_lines()[0]
    curr_df = pd.DataFrame({'Value': line2d.get_xdata(), 'Density': line2d.get_ydata()})
    st.session_state['mg__dfs_to_plot'][column_to_plot] = curr_df[(curr_df['Value'] >= column_range[0]) & (curr_df['Value'] <= column_range[1])]  # needed because the KDE can extend outside the possible value range
df_to_plot_full = st.session_state['mg__dfs_to_plot'][column_to_plot]

# Write some text boxes for the desired ranges for the current analysis column. Could make this stronger by enforcing the min_value and max_value parameters in each widget to correspond to the other so that the chosen max is never less than the chosen min, but initial attempts at this shows strange behavior and isn't worth debugging for now
st.slider(label='Minimum value:', min_value=column_range[0], max_value=column_range[1], key='mg__selected_min_val')
st.slider(label='Maximum value:', min_value=column_range[0], max_value=column_range[1], key='mg__selected_max_val')

# Create a view of the full dataframe that is the selected subset
df_to_plot_selected = df_to_plot_full[(df_to_plot_full['Value'] >= st.session_state['mg__selected_min_val']) & (df_to_plot_full['Value'] <= st.session_state['mg__selected_max_val'])]

# Plot the Plotly figure in Streamlit
fig = go.Figure()
fig.add_trace(go.Scatter(x=df_to_plot_full['Value'], y=df_to_plot_full['Density'], fill='tozeroy', mode='none', fillcolor='yellow', name='Full dataset', hovertemplate=' '))
fig.add_trace(go.Scatter(x=df_to_plot_selected['Value'], y=df_to_plot_selected['Density'], fill='tozeroy', mode='none', fillcolor='red', name='Selection', hoverinfo='skip'))
fig.update_layout(hovermode='x unified', title_text='Density vs. column value')
st.plotly_chart(fig)

# Update the run-before key
if st.session_state['mg__run_before'] == False:
    st.session_state['mg__run_before'] = True
