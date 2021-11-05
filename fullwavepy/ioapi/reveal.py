"""
(c) 2019-2020 Kajetan Chrapkiewicz.
Copywright: Ask for permission writing to k.chrapkiewicz17@imperial.ac.uk.

"""
from autologging import logged, traced
from fullwavepy.generic.parse import kw


@traced
@logged
def save_table(fname, chnls, times, delay=0):
  """
  Save a two-column table (SRC_ID, TIME) 
  of time-picks into a .tbl ASCII file 
  loadable by Reveal.
  
  Parameters
  ----------
  fname : str 
    File name to save including path
    if needed.
  chnls : list
    List of source IDs, e.g.
    [121, 122, 123]
  times : list 
    List of times IN SECONDS
    , must be of the same length 
    as chnls, e.g.
    [4.22, 3.44, 3.99]. It will be 
    converted to whole miliseconds.
  
  Notes
  -----
  Uses shot IDs.
  
  """
  f = open(fname, 'w')
  
  # WRITE HEADER
  f.write('{\n')
  f.write('"version"        : 2,\n')
  f.write('"interptype"     : "f_x",\n')
  f.write('"xkey"           : "SRC_ID",\n')
  f.write('"values"         : "TIME",\n')
  f.write('"options"        : {\n')
  f.write('        "adjust_applied_stat":false\n')
  f.write('    },\n')
  f.write('"f_of_x_picks"   :\n')
  
  f.write('  [\n')
  i = 0
  for chnl, time in zip(chnls, times):
    time = int(time * 1000) # CONVERT TO WHOLE MILISECONDS
    time += delay # DELAY (E.G. FOR BOTTOM MUTE)
    f.write('   [' + str(chnl) + '       , ' + 
            str(time) + '        ]')
    if i < len(chnls) - 1:
      f.write(',\n')
    else:
      f.write('\n')
    i += 1
  
  f.write('  ]\n')
  
  # DON'T FORGET ABOUT THE TRAILING CURLY BRACKET
  f.write('}\n')
  f.close()


# -------------------------------------------------------------------------------

