#!/usr/bin/env python
'''
Name:          launch_visar_analysis.py
License:       MIT
Version:       2.0
Last modified: 30 Jan. 2019 (SBB)
Authors:       Akel Hashim          (ahashim@slac.stanford.edu)
               Bob Nagler           (bnagler@slac.stanford.edu)
               Shaughnessy Brown    (sbbrown@slac.stanford.edu)
Description:   Creates instance of LaunchVISAR class, which runs a single VISAR
               interferometry analysis on local machine and data. 
               
Notes:         (1) Because this script calls the mecvisarimage class, the 
               timing calibration must be manually changed in the 
               visar_image_analysis.py to be used with new VISAR setups 
               (line 154)
'''

'''Imports -----------------------------------------------------------------'''
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.image as mpimg

import numpy as np

try:    # For Python2
    import PIL
except: # For Python3
    try:
        import pillow
    except:
        print('Import Error: PIL/pillow is needed to run the analysis')
        sys.exit(1)

from glob import glob

from visar_plotting_tools import *
from visar_image_analysis import *

'''Global Functions --------------------------------------------------------''' 
def check(option):
    '''Checks the sign of the option.'''
    
    if option == 1:
        return True
    elif option == 0:
        return False
    
'''Class Definitions -------------------------------------------------------''' 

class LaunchVISAR:
    ''' LaunchVISAR class takes input reference and data images, runs 
        a local VISAR analysis based on the GUI-selected analysis parameters, 
        and plots according to plotting specifications passed from the GUI. 
    '''
    
    # Assign input images, analysis params, and plot specifications to instance
    def __init__(self, f_v1_ref, f_v2_ref, f_v1, f_v2, etalon_bed1, 
                 etalon_bed2, FieldOfView, implementROI, FSV_calibration, 
                 bed1_dir, bed2_dir, FringeJumpCorrection, bed1_m, bed2_n, 
                 plotRawFringe, plotFinalFSV, saveFinalFSV):
        
        # Assign passed reference and data images to self variables
        self.filename_bed1_ref = f_v1_ref
        self.filename_bed2_ref = f_v2_ref
        self.filename_bed1 = f_v1
        self.filename_bed2 = f_v2
        
        # Check that required etalon thicknesses are specified
        if etalon_bed1 and etalon_bed2:
            self.h_v1 = float(etalon_bed1)
            self.h_v2 = float(etalon_bed2)
        else:
            print('\nERROR: No etalon thickness specified for one or both of '\
                  'the beds. VISAR analysis cannot be completed.\n')
            raise ValueError('No etalon thickness specified for one or both of'
                             ' the beds. VISAR analysis cannot be completed.')
        
        # Check if FOV is selected; set to default 256 micron if not specified            
        if FieldOfView:
            self.FOV = float(FieldOfView)
        else:
            print('No Field of View (FOV) specified. Setting default value of'\
                  ' 265 microns.\n')
            self.FOV = 265.0
        
        # Check if ROI selection is selected; use full FOV for analysis if not 
        # checked in GUI
        self.ROI = check(implementROI)
        
        # Check if FSV directions are specified; use defaults (-1 for bed1 and 
        # +1 for bed2) if not specified in GUI
        self.FSV_calib = check(FSV_calibration)
        if self.FSV_calib:
            if bed1_dir:
                self.v1_dir = int(bed1_dir) 
            else:
                print('No calibration direction specified for bed1. Setting ' \
                      'to default value of -1.\n')
                self.v1_dir = -1
            if bed2_dir:
                self.v2_dir = int(bed2_dir) # -1 or 1
            else:
                print('No calibration direction specified for bed2. Setting ' \
                      'to default value of 1.\n')
                self.v2_dir = 1

        # Check if fringe jump correction specified; set to default m=0, n=0 
        # if not specified in GUI for later iteration
        self.fjc = check(FringeJumpCorrection)
        if self.fjc:
            if bed1_m:
                s = bed1_m.split(':')
                if len(s) == 1:
                    self.m = int(bed1_m)
                    self.mh = self.m
                else:
                    self.m = int(s[0])
                    self.mh = int(s[1])
            else:
                print('No fringe jump correction m specified for bed 1. '     \
                      'Setting a default value of m = 0.\n')
                self.m = 0
                self.mh = 0
            if bed2_n:
                s = bed2_n.split(':')
                if len(s) == 1:
                    self.n = int(bed2_n)
                    self.nh = self.n
                else:
                    self.n = int(s[0])
                    self.nh = int(s[1])
            else:
                print('No fringe jump correction n specified for bed 2. '     \
                      'Setting a default value of n = 0.\n')  
                self.n = 0
                self.nh = 0
        
        # Assign plotting specifications passed from GUI       
        self.plotRF  = check(plotRawFringe)     # Plot raw fringe files
        self.plotFSV = check(plotFinalFSV)      # Plot final FSV line-outs
        self.saveFSV = check(saveFinalFSV)      # Save final FSV line-outs .csv

    def __call__(self):
        '''Calls main routine'''
        
        self.main()
    
    def main(self):
        '''Main routine'''
                
        # Initiate a VISAR image instance for both beds -----------------------
        im1 = MEC_VISAR_image_bed1(self.filename_bed1)
        im2 = MEC_VISAR_image_bed2(self.filename_bed2)
        im1_ref = MEC_VISAR_image_bed1(self.filename_bed1_ref)
        im2_ref = MEC_VISAR_image_bed2(self.filename_bed2_ref)
        
            
        # Set meters per pixel for both beds ----------------------------------
        im1._FOV = self.FOV                             # Field of View [um]
        im2._FOV = self.FOV                             # Field of View [um]
        im1._upp = im1._FOV / im1._proc_im.shape[1]     # microns per pixel
        im2._upp = im2._FOV / im2._proc_im.shape[1]     # microns per pixel
        im1._mpp = im1._upp * 1e-6                      # meters per pixel
        im2._mpp = im2._upp * 1e-6                      # meters per pixel
        im1_ref._FOV, im2_ref._FOV = im1._FOV, im2._FOV
        im1_ref._mpp, im2_ref._mpp = im1._mpp, im2._mpp
        
        # Set t0 point for each bed -------------------------------------------
        im1._t0 = 0
        im2._t0 = 0

        # Calculate velocity-per-fringe for both beds -------------------------
        im1._vpf = cal_vpf(self.h_v1)                   # [m/s]
        im2._vpf = cal_vpf(self.h_v2)                   # [m/s]
        print('VPF1 =', im1._vpf, 'm/s')
        print('VPF2 =', im2._vpf, 'm/s\n')
        
        # Initiate raw fringe files
        img1 = np.asarray(PIL.Image.open(self.filename_bed1))
        img2 = np.asarray(PIL.Image.open(self.filename_bed2))
        
        # Set common time base and field of view 
        img_time_base = im1._time_base
        # currently null
        #print(img_time_base)
        #print(type(img_time_base))
        img_FOV = im1._FOV
        
        # Reset images --------------------------------------------------------
        print('Resetting the images...\n')
        im1.reset_im()
        im2.reset_im()
        im1_ref.reset_im()
        im2_ref.reset_im()
        
        
        # Select ROI if GUI checkbox selected ---------------------------------
        if self.ROI:
            print('Selecting the ROI...\n')

            ROI_pixels = select_ROI_plots(img1,img2)
            
            im1._proc_im = im1._proc_im[ROI_pixels[2]:ROI_pixels[3],
                                        ROI_pixels[0]:ROI_pixels[1]]
            im2._proc_im = im2._proc_im[ROI_pixels[6]:ROI_pixels[7],
                                        ROI_pixels[4]:ROI_pixels[5]]
            im1_ref._proc_im = im1_ref._proc_im[ROI_pixels[2]:ROI_pixels[3],
                                                ROI_pixels[0]:ROI_pixels[1]]
            im2_ref._proc_im = im2_ref._proc_im[ROI_pixels[6]:ROI_pixels[7],
                                                ROI_pixels[4]:ROI_pixels[5]]
            
        else:
            ROI_pixels = None
        
        # Plot raw images if GUI checkbox selected ----------------------------
        if self.plotRF:
            plot_raw(img1,img2,img_FOV,img_time_base)
        
        # Resize raw fringe files if ROI; 
        # Must be performed after the raw fringe plotting routine is called
        if ROI_pixels:
            img1 = img1[ROI_pixels[2]:ROI_pixels[3],
                        ROI_pixels[0]:ROI_pixels[1]]
            img2 = img2[ROI_pixels[6]:ROI_pixels[7],
                        ROI_pixels[4]:ROI_pixels[5]]
        
        # Time calibration and FFT --------------------------------------------
        print('Getting the time calibration...\n')
        im1.get_time_cal()
        im2.get_time_cal()
        im1_ref.get_time_cal()
        im2_ref.get_time_cal()
       
        print('Performing the 1D FFT...\n')
        im1.fft_1d()
        im2.fft_1d()
        im1_ref.fft_1d()
        im2_ref.fft_1d()

        ## FFT shift ----------------------------------------------------------
        print('Performing the 1D FFT shift...\n')
        im1.fftshift_1d()
        im2.fftshift_1d()
        im1_ref.fftshift_1d()
        im2_ref.fftshift_1d()
   
        # Find frequency and apply Hanning window -----------------------------
        print('Finding the fringe frequency...\n')
        ff1     = im1.find_fringe_frequency()
        ff2     = im2.find_fringe_frequency()
        ff1_ref = im1_ref.find_fringe_frequency()
        ff2_ref = im2_ref.find_fringe_frequency()
        
        im1.window1d(ff1)
        im2.window1d(ff2)
        im1_ref.window1d(ff1_ref)
        im2_ref.window1d(ff2_ref)

        # Center frequency selection ------------------------------------------
        print('Centering the frequency selection...\n')
        im1.shift_to_center(ff1)
        im2.shift_to_center(ff2)
        im1_ref.shift_to_center(ff1_ref)
        im2_ref.shift_to_center(ff2_ref)

        # Inverse FFT shift ---------------------------------------------------
        print('Performing the 1D IFFT shift...\n')
        im1.ifftshift_1d()
        im2.ifftshift_1d()
        im1_ref.ifftshift_1d()
        im2_ref.ifftshift_1d()

        # Inverse FFT ---------------------------------------------------------
        print('Performing the 1D IFFT...\n')
        im1.ifft_1d()
        im2.ifft_1d()
        im1_ref.ifft_1d()
        im2_ref.ifft_1d()

        # FSV/Phase extraction ------------------------------------------------
        print('Performing the phase extraction...\n')
        im1.get_fsv_refl()
        im2.get_fsv_refl()
        im1_ref.get_fsv_refl()
        im2_ref.get_fsv_refl()

        # Phase unwrapping ----------------------------------------------------
        print('Unwrapping the phase...\n')
        im1.unwrap_phase2()
        im2.unwrap_phase2()
        im1_ref.unwrap_phase2()
        im2_ref.unwrap_phase2()

        # Subtract reference shot ---------------------------------------------
        print('Subtracting the reference shots...\n')
        im1._proc_refl -= im1_ref._proc_refl
        im1._proc_fsv  -= im1_ref._proc_fsv
        im2._proc_refl -= im2_ref._proc_refl
        im2._proc_fsv  -= im2_ref._proc_fsv

        # Calibrate FSV -------------------------------------------------------
        print('Calibrating the FSV...\n')
        if self.FSV_calib:
            im1.cal_fsv(direction = self.v1_dir)
            im2.cal_fsv(direction = self.v2_dir)
        else:
            im1.cal_fsv(direction = -1)
            im2.cal_fsv(direction = 1)

        # Shift to Zero -------------------------------------------------------
        print('Shifting the FSV to zero...\n')
        im1.shift_to_zero()
        im2.shift_to_zero()
        ROI_BT_pixels = None

        # Find breakout time --------------------------------------------------
        print('Finding the breakout time...\n')
        im1.find_breakout_time()
        im2.find_breakout_time()
        
        # Fringe jump correction  ---------------------------------------------
        print('Performing the fringe jump correction...\n')
        
        # If fringe jump specified, apply to data directly
        if self.fjc:
            if self.m == self.mh and self.n == self.nh:
                afjc = fringe_jump_correction(im1,im2,m=self.m,n=self.n)
            else:
                bed1 = []
                for m in range(self.m, self.mh+1):
                    bed1.append(extern_fringe_jump_correction(im1, m))
                bed2 = []
                for n in range(self.n, self.nh+1):
                    bed2.append(extern_fringe_jump_correction(im2, n))
                plot_fsv_lineout(im1,im2,plot=True,save=False,
                                 ROI_info=ROI_pixels,bed1=bed1,bed2=bed2)
                return
        
        # If no fringe jump specified, find using chisquare reduction 
        else:
            afjc = fringe_jump_correction(im1,im2)
            
        im1._proc_fsv, im2._proc_fsv = afjc[0], afjc[1]
        
        # Show and/or save the FSV line-out -----------------------------------
        if self.plotFSV: 
            if self.saveFSV:
                plot_fsv_lineout(im1,im2,plot=True,save=True,
                                     ROI_info=ROI_pixels)
            else:
                plot_fsv_lineout(im1,im2,plot=True,save=False,
                                     ROI_info=ROI_pixels)
        else:
            if self.saveFSV:
                plot_fsv_lineout(im1,im2,plot=False,save=True,
                                 ROI_info=ROI_pixels)
