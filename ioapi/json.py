"""
(c) 2019 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
import numpy as np
from autologging import logged, traced


# -------------------------------------------------------------------------------


@traced
@logged
def save_json(fname, dictionary, **kwargs):
  """
  Export a dictionary to a JSON file.
  
  """
  import json
  with open(fname, 'w') as f:
    json.dump(dictionary, f)


# -------------------------------------------------------------------------------