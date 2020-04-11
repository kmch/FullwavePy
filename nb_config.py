# -----------------------------------------------------------------------------
# Load default notebook's callables (defined in fullwavepy/__init__.py)
# -----------------------------------------------------------------------------
from fullwavepy import *

# -----------------------------------------------------------------------------
# Configure matplotlib
# -----------------------------------------------------------------------------

# Set matplotlib's backend ------

# 1. non-interactive plots, display in a notebook cell
%matplotlib inline

# 2. interactivfullwavepy.plots, display in a notebook cell
 # %matplotlib notebook

# Set matplotlib.pyplot's style ---------
plt.style.reload_library()

# Combine styles (right overwrites left wherever they overlap):
plt.style.use(['default', 'ggplot', 'kmc_test'])
# print(plt.style.available)

# -----------------------------------------------------------------------------
# Configure logging
# -----------------------------------------------------------------------------
from fullwavepy.logging_config import *
log_lvl(INFO) # TRACE / DEBUG / INFO / WARNING / ERROR / CRITICAL

# -----------------------------------------------------------------------------
# Other notebook's configuration
# -----------------------------------------------------------------------------

# autocompleting
%config IPCompleter.greedy=True 

# automatically reload modules before execution
#%load_ext autoreload
#%autoreload 2