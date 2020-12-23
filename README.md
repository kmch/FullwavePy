# FullwavePy
Python interface for full-waveform modelling and inversion

## External dependencies
Absolutely necessary:
numpy, matplotlib, logging, autologging

although it's possible (though hard) to get rid of (auto)logging dependency 
by commenting out the logged and traced decorators and replacing
_log methods with simple prints.

Needed in a few places, otherwise can be commented out 
pandas, ipywidgets, ipyvolume, ...  


## Basic usage
```Python
from fullwavepy import *
from fullwavepy.config.logging import * # this needs to be executed in jupyter cell
log_lvl(INFO) # SET DEFAULT level of log-messages to display
p = ProjSyn('test', dt=0.001, ns=1000, dx=50, dims=(101,1,41))
dir(p) # or help(p)
dir(p.inp) # ...
dir(p.out)
```

## Getting help
Open in a browser the following html files.
### Slides
nb_slides.html
### Docs
docs/built/html/index.html
