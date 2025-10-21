# Containerization Fixes for MAWA (Multiplex Analysis Web Apps)

## Overview
This document details all the changes made to fix issues with running the MAWA application in a Docker containerized environment. The primary issues were related to file path references, volume mounting, and session state management.

## Summary of Issues Fixed

1. **File path references pointing to wrong directories**
2. **Missing volume mounts for persistent data**
3. **Session state initialization timing conflicts**
4. **Framework pages failing due to session ID unavailability**

---

## 1. Docker Compose Volume Mounts

### File: `docker-compose.yml`

**Added volume mounts** to the `frontend` service to provide persistent storage and access to required directories:

```yaml
volumes:
  - ./input:/app/input
  - ./output:/app/output
  - ./saved_streamlit_session_states:/app/saved_streamlit_session_states
  - ./tmp_mawa_data:/tmp/mawa
```

**Purpose:**
- `./input:/app/input` - Allows adding input files while app is running
- `./output:/app/output` - Persists output files and results
- `./saved_streamlit_session_states:/app/saved_streamlit_session_states` - Persists session states
- `./tmp_mawa_data:/tmp/mawa` - Provides space for temporary app data

**Host directories created:**
```bash
mkdir -p output saved_streamlit_session_states tmp_mawa_data
```

---

## 2. Path Reference Fixes

### 2.1 Input Directory Path Fixes

**Files Modified:**
- `source/pages2/datafile_format_unifier.py`
- `source/pages2/open_file.py`
- `source/utils.py`
- `source/pages2/Tool_parameter_selection.py`

**Changes Made:**
```python
# Before (incorrect for container):
directory = os.path.join('.', 'input')
directory = os.path.join('..', 'input')

# After (correct for container):
directory = '/app/input'
```

**Specific Changes:**

1. **datafile_format_unifier.py** (line ~66):
```python
# Constants
directory = '/app/input'  # Changed from os.path.join('.', 'input')
```

2. **open_file.py** (line ~43):
```python
input_dir = '/app/input'  # Changed from os.path.join('..', 'input')
```

3. **utils.py** (line ~318):
```python
input_directory = '/app/input'  # Changed from os.path.join('.', 'input')
```

4. **Tool_parameter_selection.py** (line ~14):
```python
input_directory = '/app/input'  # Changed from os.path.join('.', 'input')
```

### 2.2 Output Directory Path Fixes

**Files Modified:**
- `source/benchmark_collector.py`
- `source/nidap_dashboard_lib.py`

**Changes Made:**

1. **benchmark_collector.py** (lines ~41-47):
```python
# Before:
localdir = os.path.join(utils.session_dir(), 'output')

# After:
localdir = '/app/output'
```

2. **nidap_dashboard_lib.py** (multiple locations):
```python
# Before:
session_state.OutputCSVPaths_U = os.path.join(utils.session_dir(), 'output')
session_state.checkpoint_dir = os.path.join(utils.session_dir(), 'output', 'checkpoints', 'neighborhood_profiles')
output_folder = os.path.join(utils.session_dir(), 'output')

# After:
session_state.OutputCSVPaths_U = '/app/output'
session_state.checkpoint_dir = '/app/output/checkpoints/neighborhood_profiles'
output_folder = '/app/output'
```

---

## 3. Session State Management Fixes

### 3.1 Startup Initialization Fix

**File:** `source/framework/startup.py`

**Issue:** Circular dependency where `utils.session_dir()` was called before session ID was fully available.

**Fix:** Create session directory directly instead of calling `utils.session_dir()`:

```python
# Before:
app_session_id = utils.get_unique_id()
st.session_state[ST_KEY_PREFIX + "app_session_id"] = app_session_id
session_dir = utils.session_dir()  # This could fail

# After:
app_session_id = utils.get_unique_id()
st.session_state[ST_KEY_PREFIX + "app_session_id"] = app_session_id
session_dir = f"/tmp/{utils._app_title_simple()}/app_session_data/{app_session_id}"
os.makedirs(session_dir, exist_ok=True)
```

### 3.2 Framework Pages Session Validation

**Files Modified:**
- `source/framework/manage_sessions.py`
- `source/generate_results.py`

**Added session validation function:**
```python
def check_session_initialized():
    """Check if session is properly initialized, return True if OK, False otherwise."""
    session_id_key = ST_KEY_PREFIX_STARTUP + "app_session_id"
    if session_id_key not in st.session_state:
        st.error("Session not properly initialized. Please refresh the page or return to the main page.")
        st.stop()
        return False
    return True
```

**Added to main functions:**
```python
def main():
    # Check if session is properly initialized
    if not check_session_initialized():
        return
    # ... rest of function
```

### 3.3 Generate Results Robustness

**File:** `source/generate_results.py`

**Added fallback for session directory access:**
```python
# Before:
st.session_state[key] = os.path.join(utils.session_dir(), "results", "primes", "primes.txt")

# After:
session_dir = utils.session_dir()
if session_dir is not None:
    st.session_state[key] = os.path.join(session_dir, "results", "primes", "primes.txt")
else:
    # Fallback to output directory if session dir not available
    st.session_state[key] = "/app/output/results/primes/primes.txt"
    os.makedirs("/app/output/results/primes", exist_ok=True)
```

---

## 4. Dockerfile Updates

**File:** `streamlit/Dockerfile`

**Added input directory copying** (though volume mount is preferred for development):
```dockerfile
# Copy source files to the container
COPY source/ .

# Copy input directory to the container
COPY input/ ./input/

# Copy .git to get commit hash
COPY .git .git
```

---

## 5. Container Structure

### Before Fixes:
```
Container /app/
├── app.py
├── pages2/
└── ... (other source files)
# No access to input directory
# No persistent output
# Session data lost on restart
```

### After Fixes:
```
Container /app/
├── app.py
├── pages2/
├── input/ (mounted from host)
├── output/ (mounted from host)
├── saved_streamlit_session_states/ (mounted from host)
└── ... (other source files)

/tmp/mawa/ (mounted from host)
├── app_session_data/
└── job_data/
```

---

## 6. Benefits of Changes

### 6.1 Development Experience
- ✅ **Add files while app is running** - Input files can be added to `./input/` directory
- ✅ **Persistent output** - Results saved to `./output/` directory persist between container restarts
- ✅ **Session persistence** - Session states saved in `./saved_streamlit_session_states/`
- ✅ **No rebuild required** - Changes to input data don't require container rebuild

### 6.2 Application Stability
- ✅ **Robust session handling** - Framework pages gracefully handle session initialization issues
- ✅ **Fallback paths** - App continues to work even if session directories aren't available
- ✅ **Clear error messages** - Users get helpful feedback when issues occur
- ✅ **No more path errors** - All file references point to correct container locations

### 6.3 Production Readiness
- ✅ **Consistent behavior** - App works the same way regardless of host environment
- ✅ **Data persistence** - Important data survives container restarts
- ✅ **Scalable architecture** - Volume mounts can be replaced with persistent volumes in production
- ✅ **Clean separation** - Host and container concerns are properly separated

---

## 7. Usage Instructions

### 7.1 Starting the Application
```bash
cd /path/to/multiplex-analysis-web-apps
docker compose up -d
```

### 7.2 Adding Input Files
```bash
# Copy files to input directory while app is running
cp /path/to/your/data.csv ./input/
```

### 7.3 Accessing Results
```bash
# View generated output files
ls -la ./output/
```

### 7.4 Stopping the Application
```bash
docker compose down
```

### 7.5 Rebuilding After Code Changes
```bash
docker compose build frontend
docker compose down && docker compose up -d
```

---

## 8. Testing Verification

The following functionality was tested and verified working:

- ✅ **Input file access** - Datafile Unification page can see and process files in input directory
- ✅ **Output file generation** - Results are saved to output directory
- ✅ **Session management** - Manage sessions page works without errors
- ✅ **Framework pages** - Generate results page functions correctly
- ✅ **File persistence** - Data persists across container restarts
- ✅ **Runtime file addition** - New input files can be added while app is running

---

## 9. Files Changed Summary

### Docker Configuration:
- `docker-compose.yml` - Added volume mounts

### Core Application Files:
- `source/pages2/datafile_format_unifier.py` - Fixed input directory path
- `source/pages2/open_file.py` - Fixed input directory path  
- `source/utils.py` - Fixed input directory references
- `source/pages2/Tool_parameter_selection.py` - Fixed input directory path

### Session Management:
- `source/benchmark_collector.py` - Fixed output directory path
- `source/nidap_dashboard_lib.py` - Fixed output directory paths
- `source/framework/startup.py` - Fixed session directory creation
- `source/framework/manage_sessions.py` - Added session validation
- `source/generate_results.py` - Added session validation and fallback paths

### Build Configuration:
- `streamlit/Dockerfile` - Added input directory copying

---

## 10. Maintenance Notes

### Future Development:
- When adding new pages that access file systems, use `/app/input` and `/app/output` paths
- Always check session state availability in framework pages using `check_session_initialized()`
- Test new features with container restarts to ensure data persistence

### Production Deployment:
- Consider replacing volume mounts with persistent volumes for production
- Review security implications of mounted directories
- Monitor disk usage in mounted directories
- Implement backup strategies for persistent data

---

*Document created: October 21, 2025*
*Application: MAWA (Multiplex Analysis Web Apps)*
*Environment: Docker containerized deployment*
