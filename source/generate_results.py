# Import necessary libraries.
import streamlit as st
import os
import framework.utils as framework_utils
import framework.analysis_framework as analysis_framework

# Global variable.
ST_KEY_PREFIX = 'generate_results.py__'


# Main function.
def main():

    # What does this page demonstrate?
    st.write("This page demonstrates generating results potentially asynchronously.")

    # This is a dummy widget useful for testing archive saving/loading; it does nothing.
    key = ST_KEY_PREFIX + 'number_of_results'
    if key not in st.session_state:
        st.session_state[key] = 10
    number_of_results = st.number_input("How many results would you like to generate?", min_value=1, max_value=100, key=key)

    # Output the widget value.
    st.write(f"You have requested to generate {number_of_results} results.")

    # Allow user to set input for a potentially long-running analysis (prime number generation up to this limit).
    key = ST_KEY_PREFIX + 'primes_upper_limit'
    if key not in st.session_state:
        st.session_state[key] = 4000000
    primes_upper_limit = st.number_input("Generate primes up to (inclusive):", min_value=2, key=key)

    # Run primes generation. This demonstrates the general, overall workflow for job submission.
    # Replace something like:
    #   if st.button("Run primes generation"):
    #     find_primes_up_to(limit=primes_upper_limit, results_subdir=os.path.join("results", "primes"))
    # with the following. Note you'll need to add the parameter results_topdir to the function definition and ensure it returns a dictionary of return values. (For this, it is a good idea to write a wrapper around the execution function.) If it writes files to disk, it should do so within the top directory identified by results_topdir. See analysis_functions.find_primes_up_to_limit() and the following block for a full example. Note it's best to wrap pure Python code, like the find_primes_up_to() function--not anything containing Streamlit calls.
    analysis_framework.job_submission(
        job_name="find_primes_up_to",
        inputs={"limit": primes_upper_limit, "results_subdir": os.path.join("results", "primes")},
        analysis_purpose="primes generation",
        st_key_prefix=ST_KEY_PREFIX,
    )
    # If the job has completed (session state key was set to the job results [outputs]), resume below this block.
    key = ST_KEY_PREFIX + 'primes_generation_results'
    if key not in st.session_state:
        st.warning("Primes generation results are not yet available.")
        return
    # Set shortcuts to the job results.
    primes = st.session_state[key]["primes"]
    duration = st.session_state[key]["duration"]

    # Output a note.
    st.write(f"When analysis job completed, it generated {len(primes)} primes in {duration:.2f} seconds.")

    # Obtain the path to the generated results file, stored efficiently in memory.
    key = ST_KEY_PREFIX + 'primes_results_file'
    if (key not in st.session_state) or (not os.path.exists(st.session_state[key])):  # The second condition is needed since the session-stored directory may be old and no longer current. I know this is inefficient (may as well simply have "primes_results_file = os.path.join(utils.session_dir(), "results", "primes", "primes.txt")") but it's a good conceptual example.
        st.session_state[key] = os.path.join(framework_utils.session_dir(), "results", "primes", "primes.txt")
    primes_results_file = st.session_state[key]

    # If this file, primes.txt, exists, write the contents to the screen.
    if os.path.exists(primes_results_file):
        with open(primes_results_file, "r") as f:
            st.text_area("Primes analysis results:", value=f.read(), height=300)


# Call the main function.
if __name__ == "__main__":
    main()
