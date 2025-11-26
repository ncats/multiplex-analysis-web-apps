import streamlit as st
import sample_analysis_module

ST_KEY_PREFIX = "sample_analysis_1.py__"


def main():
    
    st.session_state.setdefault(ST_KEY_PREFIX + "parameter_1", 50)
    param1 = st.number_input('Parameter 1', min_value=0, max_value=100, key=ST_KEY_PREFIX + "parameter_1")
    st.session_state.setdefault(ST_KEY_PREFIX + "parameter_2", "default text")
    param2 = st.text_input('Parameter 2', key=ST_KEY_PREFIX + "parameter_2")

    key = ST_KEY_PREFIX + "sample_analysis_results"

    if st.button('Run sample analysis'):
        results = sample_analysis_module.run_analysis(param1, param2)
        st.session_state[key] = results





    if key not in st.session_state:
        st.warning("Sample analysis results are not yet available.")
        return
    
    x = st.session_state[key]["x"]
    y = st.session_state[key]["y"]

    st.write(f'Results:\n x: {x}\n y: {y}')


if __name__ == "__main__":
    main()
