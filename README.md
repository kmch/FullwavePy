# FullwavePy
Python interface for full-waveform modelling and inversion

## External dependencies
Absolutely necessary:
numpy, matplotlib, autologging
although it's possible to get rid of autologging dependency 
by commenting the logged and traced decorators and replacing
_log methods with simple prints.

Needed in a few places:
pandas, ipywidgets, ipyvolume, ...  


## Basic usage
p = ProjSyn('test', dt=0.001, ns=1000, dx=50, dims=(101,1,41))
dir(p) # or help(p)
dir(p.inp) # ...
dir(p.out)

## Getting help
Open in a browser:
### Slides
nb_slides.html
### Docs
fullwavepy/docs/built/html/index.html
