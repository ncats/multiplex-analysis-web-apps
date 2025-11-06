# Explanation of the framework

## How things work

Whether a job is run sync or async, any job using the framework creates a job-ID-specific subdir of `framework/utils.jobs_dir()`, e.g., `/tmp/multiplex_analysis_web_apps/job_data/1234`.

An `outputs` subdirectory is tacked on and assigned to `results_topdir` in `analysis_functions.py()`, e.g., `/tmp/multiplex_analysis_web_apps/job_data/1234/outputs`.

That `results_topdir` is the top-level directory that should be assigned in any reasonable analysis function that writes files to disk. The files will then be written underneath `/tmp/multiplex_analysis_web_apps/job_data/1234/outputs`.

E.g., with how `generate_results.find_primes_up_to()` is written, it ultimately creates the file `/tmp/multiplex_analysis_web_apps/job_data/1234/outputs/results/primes/primes.txt`.

`analysis_framework.save_job_output_data()` places the serialized version of the dictionary returned by an analysis function into `/tmp/multiplex_analysis_web_apps/job_data/1234/outputs` (i.e., `results_topdir`). It then zips up that entire directory (containing both generated files and the serialized data) and puts it in the job outputs bucket.

Again whether run locally or async (behavior is always the same!), these pending jobs write a `JOB_PENDING` key to the session state. This is a flag to the app that the results ultimately need to be loaded back in.

This loading happens in `monitor_jobs.py`. When the job is complete, that page downloads that data from the job outputs bucket and puts any generated files into the session directory, e.g., `/tmp/multiplex_analysis_web_apps/app_session_data/6789`. E.g., we'd have `/tmp/multiplex_analysis_web_apps/app_session_data/6789/results/primes/primes.txt`.

There's kind of a brief summary in the comment on line 9 of `analysis_functions.py` specifically about using `results_topdir`.

## Bottom line

In general, any reasonable analysis function that writes files to disk should take a single top-level directory saying where to write the results. In the framework, this top-level dir will map to `results_topdir`.

Since that's reasonable, it would be reasonable to modify the SIT or neighborhood profiles or any analysis function to take such a top-level argument, if it doesn't already. I thought the SIT already essentially did that, but I could be wrong, it's been a while.

We should boil any analysis function down to `analysis_functions.find_primes_up_to()` in order to understand how it all works together.

Any analysis function should essentially be wrapped in a function that looks like `find_primes_up_to()`. E.g., you should always be able to put on the Streamlit page:

```python
if st.button("Run primes generation"):
    find_primes_up_to(limit=primes_upper_limit, results_subdir=os.path.join("results", "primes"), another_var=st.session_state["another_var"])
```

The actual definition of that function in `analysis_functions.find_primes_up_to()` must include the `results_topdir` parameter, whether or not it's actually used.

Then, as long as you follow the example in `generate_results.py` to replace the snippet above, everything should "just work."
