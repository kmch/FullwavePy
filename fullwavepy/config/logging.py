"""
This is a configuration file for logs from FullwavePy modules.
The engine is the 'logging' module which is a part of the Python's 
standard library.

Additionally, the 'autologging' module is used as a very convienient wrapper. 
It allows to define a function's logger (it contains both the module's and 
function's name, as defined in formatter), just by decorating this function 
with @logged. It applies also to class methods ('instance methods' of classes,
to be precise) and it's done even more simply, by decorating the class as a whole.
What's more, it provides @traced decorator (works just like @logged)
that allows tracing the execution of the code without a need for 
writing your own decorator (and logging from a function inside the 
decorator definition is challenging anyway).

Notes
-----
Here we set a default level of a root logger, its two handlers (one for each 
type of the standard streams and formatting of all messages.

List of the levels and their
corresponding numerical values:

50 CRITICAL
40 ERROR
30 WARN(ING)
20 INFO 
10 DEBUG
1  TRACE
0  NOTSET

"""
from sys import stdout, stderr
from logging import StreamHandler, Formatter, getLogger,\
  NOTSET, DEBUG, INFO, WARNING, ERROR, CRITICAL
from autologging import TRACE

# -------------------------------------------------------------------------------
# CONVENIENCE FUNCTIONS
# -------------------------------------------------------------------------------
def log_lvl(lvl):
  logger = getLogger()
  logger.setLevel(lvl)
def lll(lvl): # alias of above
  log_lvl(lvl)
# -------------------------------------------------------------------------------
# FILTERS
# -------------------------------------------------------------------------------
class LevelFilter(object):
  """
  Specify a custom logging filter to filter out 
  records with a level you don't need
  
  Notes 
  -----
  This is the format expected by 
  handler.addFilter()
  
  Filters are used primarily to filter records 
  based on more sophisticated criteria than levels
  
  """
  def __init__(self, level):
    self.level = level

  def filter(self, record):
    return record.levelno != self.level
# -------------------------------------------------------------------------------
# FORMATTING OF MESSAGES
# -------------------------------------------------------------------------------
formatter = Formatter("%(levelname)s:%(name)s.%(funcName)s: %(message)s")
# -------------------------------------------------------------------------------
# LOGGERS (ONLY ROOT LOGGER HERE)
# -------------------------------------------------------------------------------
logger = getLogger()
logger.setLevel(INFO)
# -------------------------------------------------------------------------------
# HANDLERS 
# -------------------------------------------------------------------------------
# REDIRECT TO STDERR FOR LVL >= WARN
h1 = StreamHandler(stream=stderr)
h1.setLevel(WARNING)
h1.setFormatter(formatter)
# REDIRECT TO STDOUT FOR LVL >= TRACE
h2 = StreamHandler(stream=stdout)
h2.setLevel(TRACE)
h2.setFormatter(formatter)
# EXCLUDE LEVELS HANDLED BY h1 TO PREVENT REDUNDANCY (DOUBLE MESSAGES)
h2.addFilter(LevelFilter(WARNING))
h2.addFilter(LevelFilter(ERROR))
# PUT TOGETHER
logger.handlers = [h1, h2]
# -------------------------------------------------------------------------------
# MUTE CHATTY MATPLOTLIB'S DEBUG. MESSAGES 
# -------------------------------------------------------------------------------
mpl_logger = getLogger('matplotlib.pyplot')
mpl_logger.setLevel(WARNING)