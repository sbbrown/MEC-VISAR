# MEC-VISAR
This library analyzes VISAR streak camera images to provide the velocity of a moving surface measured via optical interferometry. Users can analyze images via a built-in GUI or call the analysis class LaunchVISAR in her/his own analysis code.

## Motivation
This project is a pared-down, local version of the standard VISAR code implemented at the Matter in Extreme Conditions (MEC) endstation of the Linac Coherent Light Source (LCLS) at SLAC National Accelerator Laboratory.

VISAR  = Velocity Interferometry System for Any Reflector. VISAR systems send laser light reflected off a moving surface through an interferometer, creating interference patterns that are sampled and streaked in time. These interference patterns or “fringe patterns” shift as the velocity of the surface changes, creating “fringe shifts.” Recording fringe patterns during an experiment provides information about the velocity time history of the surface or “free-surface velocity.” 

The most common application of a VISAR system at MEC is recording the free-surface velocity of a target undergoing shock compression. When a shock propagates through a sample, the rear surface of the sample moves in time. This movement is recorded by the VISAR system. The shock velocity relates to the pressure behind the shock front via Rankine-Hugoniot relations.  

For a full explanation of the setup, mathematical formulation, and analysis of VISAR images, please see https://lcls.slac.stanford.edu/instruments/mec/visar-analysis

NOTE: If you are adapting this script to analyze VISAR images from a facility other than SLAC MEC, you will need to edit the bed specifications in visar_image_analysis.py in the classes MEC_VISAR_image_bed1 and MEC_VISAR_image_bed2 to match the camera timings of your specific facility detectors. 

## Code style
Google python style guide https://github.com/google/styleguide/blob/gh-pages/pyguide.md

## Built with
•	Python 2.7 or newer
•	matplotlib
•	numpy
•	tkinter
•	PIL
•	scipy

## Features
GUI to select VISAR streak and background images, specify necessary analysis parameters, and choose desired plots. 
Beyond performing a routine VISAR analysis, this code contains optional analysis options.
* **implementROI**: instead of analyzing the entire VISAR image, user can select only a limited portion to analyze.
* **FSV_calibration**: user can specify the direction of the initial fringe shift and overwrite the defaults.  
* **FringeJumpCorrection**: user can either specify integer values for the number of fringe shifts added to each bed (i.e. m=1, n=2) OR specify a range (m=0:3, n=0:6). If the range option is utilized, the script plots all the FSVs of all specified fringe shifts on one plot for rapid comparison. 

The analysis can also be called from user-generated code through the class LaunchVISAR in launch_visar_analysis.py

## Code example
While the primary usage is through the GUI, a user can also write her/his own code and call the VISAR analysis through the class LaunchVISAR.
```
# Import VISAR analysis library
from launch_visar_analysis import LaunchVISAR

# Specify filenames of VISAR images to be analyzed
f_v1_ref = VISAR_bed_1_filename_reference_shot.tiff
f_v2_ref = VISAR_bed_2_filename_reference_shot.tiff
f_v1     = VISAR_bed_1_filename_driven_shot.tiff 
f_v2     = VISAR_bed_2_filename_driven_shot.tiff

# Specify required analysis parameters
etalon_bed1 = 11.006  # etalon 1 thickness (mm)
etalon_bed2 = 50.002  # etalon 2 thickness (mm)
FieldOfView = 256     # VISAR field-of-view (micron)

# Specify optional analysis parameters
implementROI          = True    # Implement ROI for analysis (boolean)
FSV_calibration       = True    # Specify direction of fringe shift (boolean)
bed1_dir              = 1       # VISAR bed 1 direction (1= left, -1=right)
bed2_dir              = -1      # VISAR bed 2 direction (1= left, -1=right)
FringeJumpCorrection  = True    # Specify number of fringe jumps (boolean)
bed1_m                = 2       # VISAR bed 1 num. of fringe shifts (int)
bed2_n                = 3       # VISAR bed 2 num. of fringe shifts (int)

# Specify plotting options
plotRawFringe = False   # Show raw VISAR images (boolean)
plotFinalFSV  = True    # Show final FSV lineouts of two beds (boolean)
saveFinalFSV  = True    # Save final FSV lineouts of two beds (boolean)

# Create instance of LaunchVISAR class to perform analysis with specified inputs
AnalysisInstance = LaunchVISAR LaunchVISAR(f_v1_ref, f_v2_ref, f_v1, f_v2, 
etalon_bed1, etalon_bed2, FieldOfView, implementROI, FSV_calibration, bed1_dir, bed2_dir, 
FringeJumpCorrection, bed1_m, bed2_n, plotRawFringe, plotFinalFSV, saveFinalFSV)
```

## Installation
This is a library, so download the scripts (6 total). To run the GUI-based analysis, navigate in a terminal (mac) or command line (windows) and run the GUI script.

FOR MACS ``` python visar_gui.py ```
FOR WINDOWS
```
py visar_gui.py
```
## Credits
This code is based upon the standard VISAR code implemented at the Matter in Extreme Conditions (MEC) endstation. The standard code was written by Akel Hashim based on scripts provided by Bob Nagler. Any and all bugs in this offline version are the responsibility of Shaughnessy Brennan Brown : )

## License
MIT © Shaughnessy Brennan Brown
