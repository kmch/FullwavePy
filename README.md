# :white_square_button: FullwavePy
Python interface for full-waveform modelling and inversion. Currently only for the proprietary `FULLWAVE3D` code which can be accessed through academic collaboration with Imperial College London or a membership of the FULLWAVE consortium (Warner et al., 2013). Serving other codes can be easily accomplished through suitable APIs.

## Getting started
### Prerequisite Python packages
#### Required
`numpy`, `matplotlib`, `logging`, `autologging`

It's possible (though hard) to get rid of (auto)`logging` dependency 
by:
- commenting out the *logged* and *traced* decorators and 
- replacing *_log* methods with simple *print*.

#### Optional
`pandas`, `ipywidgets`, ...

These are needed to access advanced features and can be commented out for simple applications.

### Prerequisite binaries
#### Required
`vtr2sgy`, `sgy2vtr`, `SegyPrep`

First two files are a temporary method to handle I/O of SEGY files. The third is the pre-processor of `FULLWAVE3D`.  They have to be included in the path (defined in .bashrc or .profile for Linux), i.e. be callable in the Linux/Mac console
without specifying the path to it. All these are Fortran codes shipped with `FULLWAVE3D` and are not provided here.

#### Optional
`SeismicUnix`

Popular C package for seismic data processing. May be skipped unless data filtering and muting is intended to be applied to the data.

## Basic usage
```Python
from fullwavepy import *
from fullwavepy.config.logging import * # this needs to be executed in jupyter cell
log_lvl(INFO) # SET DEFAULT level of log-messages to display
p = ProjSyn('test', dt=0.001, ns=1000, dx=50, dims=(101,1,41))
help(p)
dir(p)  
dir(p.inp)
dir(p.out)
```

## Getting help

To get more information, open in a browser the following html files:
- *presentations/nb_slides.html* - this may be useful to get the idea of how to use the package, even though
some code blocks are in fact pseudo code.

- *docs/built/html/index.html* - home page of documentation (WIP)

## References
Warner, M. et al. (2013) *Anisotropic 3D full-waveform inversion*, Geophysics, vol. 78, p. R59-R80, <a href="https://library.seg.org/doi/10.1190/geo2012-0338.1">doi:10.1190/geo2012-0338.1</a>

