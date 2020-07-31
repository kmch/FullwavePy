
# -----------------------------------------------------------------------------
# Load default modules defined in fullwavepy/__init__.py
# -----------------------------------------------------------------------------
from fullwavepy import *


# -----------------------------------------------------------------------------
# Set aliases of frequently used jupyter magic commands
# -----------------------------------------------------------------------------
%alias_magic mi matplotlib -p inline
%alias_magic mn matplotlib -p notebook


# -----------------------------------------------------------------------------
# Configure matplotlib
# -----------------------------------------------------------------------------
# SET DEFAULT backend:
# non-interactive plots displayed in a notebook cell
%matplotlib inline
# OTHER OPTIONS: %matplotlib notebook (interactive version)
# -----------------------------------------------------------------------------
# SET DEFAULT style:
plt.style.reload_library()
# a combined style (right overwrites left wherever they overlap):
plt.style.use(['default', 'ggplot', 'kmc_test'])
# OTHER OPTIONS: print(plt.style.available)


# -----------------------------------------------------------------------------
# Configure logging
# -----------------------------------------------------------------------------
# Set up loggers, handlers etc. 
# and load the log_lvl function
from fullwavepy.config.logging import *
# -----------------------------------------------------------------------------
# SET DEFAULT level of log-messages to display
log_lvl(INFO)
# OTHER OPTIONS: TRACE, DEBUG, INFO, WARNING, ERROR, CRITICAL
# or using integers: 0, 10, 20, 30, 40, 50 respectively
# (in order of increasing importance/dicreasing verbosity)


# -----------------------------------------------------------------------------
# Other notebook's configuration
# -----------------------------------------------------------------------------
# Autocompleting
%config IPCompleter.greedy=True 
# -----------------------------------------------------------------------------
# Automatically reload modules before execution
#%load_ext autoreload
#%autoreload 2