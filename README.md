# VLBIana


This repository combines a selection of analysing tools for VLBI data, written in Python. It is still in testing mode and mostly meant to distribute these files to colleagues.

By Anne-Kathrin Baczko (baczko@mpifr-bonn.mpg.de).

If you like to use these scripts for any publication, please let me know.

It has only been tested on a small amount of datasets. Hence, let me know if somethings does not work and we find a solution.

For the imported modules to be found the repository has to been added to the local PYTHONPATH variable.

## Requirenments
- astropy, numpy, scikit-image
- The script align_imagesEHTim.py requires the eht-imageing module, which can be obtained from https://github.com/achael/eht-imaging.git
- The script plotSet.py for setting the figure dimensions is taken from https://jwalton.info/Embed-Publication-Matplotlib-Latex/
- Please save pubstyle.mplstyle in a location were your system looks for it, e.g. $PATH or better ~/.config/matplotlib/stylelib for some matplotlib settings to produce nice, publicatoin ready plots. This file is loaded in modules/plot_functions.py. Without it the plots might look weird.

### Ridge-line analysis
- ParselTongue: for information about installation please see https://github.com/kernsuite-debian/parseltongue 
- AIPS: in case you do not want to use AIPSLite (which comes with ParselTongue) but a local installation. Please see http://www.aips.nrao.edu/index.shtml for installation instructions
- VLBIcalib: to be obtained at https://github.com/abaczko/VLBIcalib.git

