.. fullwavepy documentation master file, created by
   sphinx-quickstart on Sun May 24 16:50:54 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

FullwavePy: FWI workflow in Jupyter
===================================

Key features
------------
Plotting i/o

Key concepts
------------
FWI is conceptually simple but usually requires a lot of data manipulation
and parameter tuning to make it work successfully. Hence a need for 
convenient, robust and extensible framework that automates the process and 
allows to focus on research questions rather than technicalities.

FullwavePy's project is something defined uniquely by FWI **input** files. 


E.g. a project has exactly one *Runfile* but can have multiple auxiliary, job-specific files (e.g. *Out.log*, *Run.pbs*, *JobInfo.txt*, etc.) 
that have job-file identifiers embedded in their names (actual names *Out0.log*, etc.). Obviously, **output files** are job-specific too, but currently they are **overwritten after each job-run** instead of being endowed with such ids (to save disk space). To preserve them, one needs to create a new project.


The basic building blocks of any project are files. 
They are abstract objects i.e.
they are not bound to any specific I/O implementation.

The binding needed is done under the hood

and the tricky part is to get 
all the necessary data regardless of the io implementation.
This can be challenging e.g. when no headers are associated
with seismic data. An alternative way of getting this data
has to be implemented. Or the user should be informed of the 
limited functionality.

Plotting of the data should be independent of the io.


Requirements
------------
Remember to add parent-path/fullwavepy to your PYTHONPATH.
(parent-path is not enough since we have:
parent-path/fullwavepy/fullwavepy/__init__.py).

If the autogen fails to find one of the modules unlike all 
the other modules, most likely there is a syntax error inside 
this module. Try to import it from python interpreter to to track
down the bug.

.. autosummary::
   :toctree: modules
   
   fullwavepy.dsp.generic
   fullwavepy.dsp.phase
   fullwavepy.dsp.su
   fullwavepy.dsp.wavelet
   fullwavepy.fd.checks
   fullwavepy.fd.stencils
   fullwavepy.generic.decor
   fullwavepy.generic.parse
   fullwavepy.generic.system
   fullwavepy.ioapi.clusters.archer
   fullwavepy.ioapi.clusters.cx1
   fullwavepy.ioapi.clusters.thomas
   fullwavepy.ioapi.fw3d
   fullwavepy.ioapi.generic
   fullwavepy.ioapi.memmap
   fullwavepy.ioapi.proteus
   fullwavepy.ioapi.reveal
   fullwavepy.ioapi.segy
   fullwavepy.ioapi.su
   fullwavepy.ndat.arrays
   fullwavepy.ndat.manifs
   fullwavepy.ndat.points
   fullwavepy.numeric.const
   fullwavepy.numeric.fourier
   fullwavepy.numeric.funcs
   fullwavepy.numeric.generic
   fullwavepy.numeric.operators
   fullwavepy.plot.events
   fullwavepy.plot.generic
   fullwavepy.plot.misc
   fullwavepy.plot.plt1d
   fullwavepy.plot.plt2d
   fullwavepy.plot.plt3d
   fullwavepy.project.files.datalike.generic
   fullwavepy.project.files.datalike.sgy
   fullwavepy.project.files.datalike.ttr
   fullwavepy.project.files.generic
   fullwavepy.project.files.gridded.derivs
   fullwavepy.project.files.gridded.generic
   fullwavepy.project.files.gridded.misc
   fullwavepy.project.files.gridded.models
   fullwavepy.project.files.gridded.surfaces
   fullwavepy.project.files.gridded.wavefields
   fullwavepy.project.files.other.ghost
   fullwavepy.project.files.other.index
   fullwavepy.project.files.text.hed
   fullwavepy.project.files.text.logs
   fullwavepy.project.files.text.misc
   fullwavepy.project.files.text.runfiles
   fullwavepy.project.files.text.srcrec
   fullwavepy.project.files.text.submit
   fullwavepy.project.generic.au
   fullwavepy.project.generic.io
   fullwavepy.project.generic.qc
   fullwavepy.project.lists.basic
   fullwavepy.project.lists.deriv
   fullwavepy.project.lists.extra
   fullwavepy.project.types.basic
   fullwavepy.project.types.deriv
   fullwavepy.project.types.extra
   fullwavepy.seismic.data
   fullwavepy.seismic.metadata
   fullwavepy.seismic.models
   fullwavepy.seismic.srcrec
   fullwavepy.seismic.surfaces
   fullwavepy.seismic.wavefields
   fullwavepy.seismic.wavelets



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
