# export PATH=/home/andrew/apps/micromamba/bin:$PATH
# cd deploy/snowflake
# micromamba create -n basic -f environment.yml
# micromamba run -n basic streamlit run snowflake_compute_pool_instance_types.py

import streamlit as st
import pandas as pd
import plotly.express as px

df = pd.read_excel("snowflake_compute_pool_instance_types.xlsx")

# Sort by memory so the connecting line progresses logically left->right.
df = df.sort_values("memory (gib)")

# Ensure 'fold increase' displays with a single decimal place.
if "fold increase" in df.columns:
    df["fold increase"] = df["fold increase"].astype(float).round(1)

# Create a plotly chart of "credits per hour" vs. "memory (gib)" with hover info showing "name".
fig = px.scatter(
    df,
    x="memory (gib)",
    y="credits per hour",
    hover_data = ["name", "fold increase", "memory (gib)", "credits per hour", "vcpu"],
    title="Snowflake Compute Pool Instance Types: Credits per Hour vs. Memory (GiB)",
    labels={
        "memory (gib)": "Memory (GiB)",
        "credits per hour": "Credits per Hour",
    },
)

# Draw a line connecting the points while retaining markers.
fig.update_traces(mode="lines+markers")
st.plotly_chart(fig)
st.dataframe(df)
st.write("Data source: Snowflake documentation on Compute Pool Instance Types from https://www.snowflake.com/legal-files/CreditConsumptionTable.pdf.")

st.write("Note the two HIGHMEM variants are actually great value for the memory provided.")

st.write("Note also the fold increase refers to cost relative to the XS instance type.")
