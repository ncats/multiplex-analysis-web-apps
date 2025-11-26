import streamlit as st

ST_KEY_PREFIX = "sample_analysis_1.py__"


def run_analysis(param1, param2, results_topdir=None):
    x = param1 * 2
    y = param2.upper()
    return {"x": x, "y": y}


def main():
    
    param1 = st.number_input('Parameter 1', min_value=0, max_value=100, value=50)
    param2 = st.text_input('Parameter 2', value='default text')

    key = ST_KEY_PREFIX + "sample_analysis_results"

    if st.button('Run sample analysis'):
        results = run_analysis(param1, param2)
        st.session_state[key] = results





    if key not in st.session_state:
        st.warning("Sample analysis results are not yet available.")
        return
    
    x = st.session_state[key]["x"]
    y = st.session_state[key]["y"]

    st.write(f'Results:\n x: {x}\n y: {y}')


if __name__ == "__main__":
    main()
