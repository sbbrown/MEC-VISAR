#!/usr/bin/env python
'''
Name:          launch_visar_gui.py
License:       MIT
Version:       1.0
Last modified: 9 Aug. 2018 (SBB)
Authors:       Bob Nagler           (bnagler@slac.stanford.edu)
               Akel Hashim          (ahashim@slac.stanford.edu)  
               Shaughnessy Brown    (sbbrown@slac.stanford.edu)     
Description:   Runs VISAR analysis on local machine and data
'''

'''Imports -----------------------------------------------------------------'''
from Tkinter import *
from tkFileDialog import askopenfilename

from matplotlib import pyplot
from matplotlib import image as image_p

from numpy import *
from scipy.signal import argrelmax
from scipy.ndimage.filters import gaussian_filter, gaussian_filter1d

from PIL import Image as image

from bisect import bisect_left
import itertools

from visar_bed import *

'''Global Functions --------------------------------------------------------''' 

def cal_vpf(h):
    '''Calculates velocity per fringe for both beds based on respective etalon 
    thicknesses'''
   
    h *= 1e-3                       # Etalon thickness for bed1 in [m]

    # Initate VISAR bed instance
    vpf_i = VISAR_Bed(h)
    
    # Calculate the vpf
    vpf = vpf_i.vpf0()
    
    return vpf
    
def find_breakout_index(fsv_lineout):
    '''Finds the breakout index for a given FSV line-out'''
    
    if fsv_lineout.size == 1024:    # No ROI specified
        sigma = 5
        num_points = 200
        radius = 25
        n = 4
    elif fsv_lineout.size < 1024:   # ROI is specified
        sigma = 5
        num_points = 75
        radius = 25
        n = 2.5
        
    # Apply a Gaussian filter to the line-out
    filtered_data = gaussian_filter1d(fsv_lineout,sigma) 
    
    # Take the derivative of the filtered data
    deriv  = gradient(filtered_data) 
    aderiv = abs(deriv)
    
    # Take the standard deviation of the first 200 points
    stdev  = std(deriv[:num_points]) 
    
    # Locate the indices of the relative maxima spanning 25 data points
    peak_indices = argrelmax(aderiv,order=radius) 
    
    # Peak values for the peak indices
    peaks  = deriv[peak_indices] 

    # Find the first value that is > than n x the standard deviation
    for index,value in enumerate(peaks):
        if abs(value) > n*stdev: 
            breakout, breakout_value = index, value
            break
    
    try:
        breakout_index = peak_indices[0][breakout]
    except:
        breakout_index = 1 
        
    return breakout_index

def fringe_jump_correction(im1,im2,m=None,n=None):
    '''Fringe jump correction routine. Determines the fringe jump 
       correction needed by iterating through different combinations of m,n 
       values and minimizing the chisquare value of resulting FSV line-outs.
    '''
    
    # Create array to iterate over of the form [0,...,10]
    fringe_jumps = linspace(0,10,11)
    
    # Iterate through linspace of possible combinations
    if m == None or n == None:
        chisq_values = {}
        permutations = list(itertools.product(fringe_jumps,repeat=2))
        
        for i in permutations:
            m, n  = i[0], i[1]
            chisq = apply_fringe_jump_correction(im1,im2,m,n)
            chisq_values[i] = chisq

        # Find m,n values for minimum chisq coefficient
        m_n_values = min(chisq_values, key=chisq_values.get)
        m = m_n_values[0]
        n = m_n_values[1]

    # Apply fringe jump corrections using m,n values corresponding to maximum 
    # Pearson coefficient
    for i in arange(im1._proc_fsv.shape[1]):
        for j in arange(int(im1._breakout_indices[i])+1,
                            im1._proc_fsv.shape[0]):
            # Add fringe jump correction to im1
            im1._proc_fsv[j,i] += (m * im1._vpf) 
        
    for i in arange(im2._proc_fsv.shape[1]):
        for j in arange(int(im2._breakout_indices[i])+1,
                            im2._proc_fsv.shape[0]):
            # Add fringe jump correction to im2
            im2._proc_fsv[j,i] += (n * im2._vpf) 
    
    print 'm =', m
    print 'n =', n, '\n'
    return (im1._proc_fsv, im2._proc_fsv)
     
def apply_fringe_jump_correction(im1,im2,m,n):
    '''Iterates through the number of fringe jumps in each bed such that the 
       chisquare value of the FSV of both beds is minimized.'''
    
    xdim1 = im1._proc_fsv.shape[1]
    xdim2 = im2._proc_fsv.shape[1]
    
    # Take the fsv line-out in middle of image
    fsv1_i = im1.fsv_lp(show_plot=False,return_array=True)
    fsv2_i = im2.fsv_lp(show_plot=False,return_array=True)
    
    vpf1 = im1._vpf                     # Velocity per fringe for bed1
    vpf2 = im2._vpf                     # Velocity per fringe for bed2
    
    # Breakout indices of each bed
    breakout_index1 = int(im1._breakout_indices[int(xdim1/2)]) 
    breakout_index2 = int(im2._breakout_indices[int(xdim2/2)]) 
    
    # Separate the fsv of each array at the breakout point
    fsv1_0 = fsv1_i[0:breakout_index1]
    fsv1_1 = fsv1_i[breakout_index1:]
    fsv2_0 = fsv2_i[0:breakout_index2]
    fsv2_1 = fsv2_i[breakout_index2:]
    
    # Apply the fringe jump correction according to the m,n values passed    
    fsv1_1 = fsv1_1 + (m * vpf1)
    fsv2_1 = fsv2_1 + (n * vpf2)

    fsv1 = concatenate((fsv1_0, fsv1_1))
    fsv2 = concatenate((fsv2_0, fsv2_1))
    
    return chisquare(fsv1,fsv2,breakout_index1,breakout_index2)

def chisquare(fsv1,fsv2,bi1,bi2):
    '''Computes the chi square value for the two FSV line-outs'''
    
    chisq = 0
    if fsv1.shape[0] == fsv2.shape[0]:
        for i in arange(fsv1.shape[0]):
            chisq += (fsv1[i] - fsv2[i])**2
    else:
        # Shift the arrays so that both breakout indices match up
        if bi1 > bi2:
            fsv1 = fsv1[bi1-bi2:]
        elif bi2 > bi1:
            fsv2 = fsv2[bi2-bi1:]
        # Truncate the arrays so they are the same length       
        if fsv1.shape[0] > fsv2.shape[0]:
            diff = fsv1.shape[0] - fsv2.shape[0]
            fsv1 = fsv1[:-diff]
        elif fsv2.shape[0] > fsv1.shape[0]:
            diff = fsv2.shape[0] - fsv1.shape[0]
            fsv2 = fsv2[:-diff]
        for i in arange(fsv1.shape[0]):
            chisq += (fsv1[i] - fsv2[i])**2
                       
    return chisq

def extern_fringe_jump_correction(im1, m):
    ''' Creates a plot if range of m and n fringe jump corrections specified
    '''
    
    proc_fsv = copy(im1._proc_fsv)
    for i in arange(im1._proc_fsv.shape[1]):
        for j in arange(int(im1._breakout_indices[i])+1,im1._proc_fsv.shape[0]):
            # Add fringe jump correction to im1
            proc_fsv[j,i] += (m * im1._vpf) 
    return proc_fsv

  
'''Class Definitions -------------------------------------------------------'''
    
class VISAR_image():
    '''Class that analyzes raw visar image

       Calculates free surface velocity based on fringe pattern of raw streak 
       image using the fourier transform method.
       All units are SI unless specifically noted

       class objects/functions:
       
       self._fringe     : array of original fringe pattern
       self._proc_fsv   : array of phase or vpf after reconstruction
       self._proc_refl  : array of target reflectivity after reconstruction
       self._vpf        : velocity per fringe of the visar bed, in m/s
       self._time_coeff : dict of the coeff. of timing differential polynomial
       self._time_base  : sweep time of current image streak; read from tiff

       self.read_fringe_file(filename)  : reads a tiff fringe file from disk
       self.shape()                     : returns shape of _proc_im array
       
       self.fft1d()                     : performs 1D fft of _proc_im over 
                                          x-axis
       self.fftshift1d()                : shifts _proc_im s.t. zero frequency 
                                          is middle of image
       self.find_fringe_frequency()     : returns frequency of the fringes
       
       self.window1d(fringe_freq)       : hanning window of _proc_im around 
                                          fringe_freq
       self.shift_to_center(fringe_freq): shifts _proc_im s.t. fringe_freq in 
                                          center window
       self.ifftshift()                 : performs inverse fftshift; zero 
                                          frequency is now at pixel 0
       self.ifft1d()                    : performs inverse  fft1d of _proc_im
       
       self.get_fsv_refl()              : converts phase (intensity) of 
                                          _proc_im into _proc_vpf (_proc_refl)
       self.fsv_lp(xmin,xmax)           : provides lineout of fsv between xmin 
                                          and xmax; middle 1% is taken if None
    '''

    def __init__(self,filename=None):
        '''Initializes the object.'''

        self._lambda=532e-9 # laser wavelength

        # Read in filename and initialize variables
        filename = self.read_fringe_file(filename)
        if filename == None:
            self._fringe    = None     # Array (real) of original fringe image
            self._proc_im   = None     # Array (complex) of image being proces.
            self._time_base = None     # Time base of Hamamatsu streak camera
        self._proc_fsv  = None     # Array (real) of the fsv being processed
        self._proc_refl = None     # Array (real) of target reflectivity
        self._vpf       = None     # velocity per fringe (real)
        self._breakout  = None     # Array (int) with breakout position. 
                                   #      all element 0 except 1 where breakout
        self._mpp = None           # Micron per pixel in the x axis
        self._t0  = 0              # Pixel of t0; usually rising edge of laser
        self._time_coeff = {}      # Empty dictionary for time coefficients
        self._time_cal   = None    # Time calibration array

    def read_fringe_file(self,filename=None):
        ''' Reads .tiff file from the Hamamtsu streak camera program

            Parses tags of the tiff for the time range
        '''
        
        # Determine filename if None was passed
        if filename == None:
            Tk().withdraw()
            filename_new = askopenfilename(title='Select data .tiff', 
                                           filetypes=[('tif files','*.tif'),
                                                      ('All files', '*')])
        else:
            filename_new = filename
        
        # Open image to be processed
        if filename_new=='' or filename_new==():
            return None
        else:
            im  = image.open(filename_new)
            iml = reshape(list(im.getdata()),(im.size[1],im.size[0]))
            self._fringe  = where(iml<0,iml+65536,iml)
            self._proc_im = self._fringe
            
            # Parse tag data for timing
            info = im.tag.tagdata[270]
            r_i  = info.find('Time Range')
            r_e  = info.find('ns',r_i)
            
            # Set time base
            self._time_base = info[r_i+12:r_e+2]
            if self._time_base not in ('0.5 ns','1 ns','2 ns','5 ns','10 ns',
            '20 ns','50 ns','100 ns','200 ns','500 ns'):
                self._time_base = "None"
        
        return filename_new

    def reset_im(self):
        ''' Reset proc_im to the original fringe image'''
        
        self._proc_im = self._fringe      
        
    def get_time_cal(self):
        '''Calculates time of each pixel and places in array self._time_cal
            1) builds difference polynomial, if the time range in dict
            2) builds array of the time for each pixel
            3) reference to t0
            4) time cal array saved in self._time_cal
        '''
        # Check if time base assigned from tags
        if self._time_base == None:
            print 'No time base defined. aborting'
            return
        if self._time_base not in self._time_coeff:
            print 'No coefficients found for time base of '+ self._time_base
            return
        
        # Assign shape of image to local variable
        sizex = self.shape()[1]
        sizey = self.shape()[0]
        
        # Create polynomial coefficients
        f = poly1d(self._time_coeff[self._time_base])
        self._time_cal = zeros(sizey)
        
        for index in range(1,sizey):
            self._time_cal[index] = self._time_cal[index-1]+f(index)
        
        self._time_cal = self._time_cal 
        
    def shape(self):
        '''Returns the shape of the _proc_im'''
        
        return shape(self._proc_im)
   
    def fft_1d(self,show_plot=True):
        '''Performs one-dimensional fft on processed image over the x-axis'''
        
        self._proc_im = fft.fft(self._proc_im, axis=1)

    def fftshift_1d(self):
        '''Shifts zero freq to middle of array, 1d shift'''
        
        self._proc_im = fft.fftshift(self._proc_im,axes=1)

    def find_fringe_frequency(self):
        '''Finds pixel value corresponding to fringe frequency

           Should be run on fft of the fringes, which is shifted to get zero 
           frequency at the center of the image. This shifted fft is projected 
           on the x axis. This gives a maximum in the center (corresponding to 
           the constant background) and two side maxima, corresponding to the 
           fringe frequency. The pixel of the sidemaxima is found by finding 
           the index of the maximum starting from a number of pixels to the 
           right of the center. This number of pixels is the offset. 
        '''
        
        xpr = sum(self._proc_im,axis=0) # Projection on x axis
        center_pix = size(xpr)/2.0
        offset = 5
        return center_pix+offset+argmax(abs(xpr[center_pix+offset:]))

    def window1d(self,fringe_freq,width=None):
        '''Windows shifted 1dfft around the passed fringe_freq.

           Window is simple hanning window, centered at the fringe_freq.
           The width is taken as distance between the fringe freq and the
           zero frequency. The window is flat in the y direction.
        '''
        
        if width == None:
            width = int(fringe_freq-shape(self._proc_im)[1]/2)
        
        width  = 2*(width/2)+1      # add 1 to ensure the width is odd
        wind_h = hanning(width)
        
        wind1d = zeros(self.shape()[1])
        wind1d[fringe_freq-width/2:fringe_freq+width/2+1]=wind_h
        wind2d = zeros(self.shape())
        wind2d[:,:]   = wind1d
        self._proc_im = self._proc_im*wind2d

    def shift_to_center(self,fringe_freq):
        '''Shifts fringe_freq to center of image, only in x direction.

           The function shifts the fringe_freq to the center of the image.
           It assumes proc_im is an fft that has been windowed and that the 
           shift will always be to lower indices. The shift drops the values
           that would be shifted to negative indices (left side of image) and 
           pads the right side of the image with zeros to get the same size. 
           This should be the same as the cyclic shift; outside of a small 
           window there should only be zeros.
        '''
        
        # Amount the image will be shifted
        try:
            d_shift = fringe_freq-self.shape()[1]/2 
            self._proc_im[:,
                          0:self.shape()[1]-d_shift]=self._proc_im[:,d_shift:]
        except:
            d_shift = int(fringe_freq-self.shape()[1]/2) 
            self._proc_im[:,
                          0:self.shape()[1]-d_shift]=self._proc_im[:,d_shift:]
            
    def ifftshift_1d(self):
        
        '''Performs 1dfft shift; the zero frequency is again at position 0'''
        self._proc_im = fft.ifftshift(self._proc_im,axes=1)

    def ifft_1d(self):
        
        '''Performs 1Difft of  _proc_im, over x-axis'''
        self._proc_im = fft.ifft(self._proc_im,axis=1)

    def get_fsv_refl(self):
        '''Takes _proc_im and puts phase->_proc_vpf and intens.->_proc_refl'''
        
        self._proc_fsv  = angle(self._proc_im)
        self._proc_refl = abs(self._proc_im)
        
    def unwrap_phase2(self):
        '''Unwrap phase routine using numpy.unwrap
           Unwrap phase of each _proc_fsv array in space.'''
        
        # Unwrap the phase in space
        for i in arange(self._proc_fsv.shape[0]):            
            self._proc_fsv[i,:] = unwrap(self._proc_fsv[i,:])

        # Unwrap the phase in time
        for i in arange(self._proc_fsv.shape[1]):            
            self._proc_fsv[:,i] = unwrap(self._proc_fsv[:,i])        

    def fsv_lp(self,xmin=None,xmax=None,in_pix=False,show_plot=True,
               return_array=False):
        '''Takes the avg line profile of the fsv over y

           Intgrates fsv over x (space axis) and plots fsv vs y (time axis)
           If no values are passed, the middle 1% of the image is taken. 
           If show_plot=False, the plot isn't shown (for use with subplot).
           If in_pix=True, time_calibration isn't used and the x-axis of the
           plot is in pixels.
        '''
        # Set boundaries if none were passed (or only one was passed)
        if (xmin == None or xmax == None):
            xmin = int(self.shape()[1]*0.495)
            xmax = int(xmin/0.495*0.505)

        # Create 1d array, y(the fsv) and the x (time)
        lp = mean(self._proc_fsv[:,xmin:xmax],1)

        if show_plot :
            pyplot.show()

        if return_array:
            return lp

    def cal_fsv(self,direction=-1):
        '''Calibrates fsv with self._vpf after reconstruction and phase 
           unwrapping

           Program should be run after all the phase reconstruction and 
           unwrapping completed. Function is simple multiplication with 
           self._vpf (velocity per fringe).
           If the direction=1, increasing phase means increasing velocity
           If the direction=-1, increasing phase means decreasing velocity
           
           The standard direction is -1, indicating fringes shift to right on
           breakout; this is the standard at MEC.
        '''
        
        if self._vpf == None:
            print 'No velocity per fringe defined; cannot calibrate phase.'
        else:
            self._proc_fsv = self._proc_fsv*direction*self._vpf/(2*pi)

    def find_breakout_time(self,sigma=None):
        '''Determines where the breakout is in reconstructed ._proc_fsv. 
           There is potentially a 2pi fringe shift uncertainty.'''
        
        xdim = self._proc_fsv.shape[1]
        ydim = self._proc_fsv.shape[0]

        # Check is ROI was used
        if ydim == 1024:
            ROI = False
        elif ydim < 1024:
            ROI = True

        if ROI == False:
            up = 20
            down = 60
            n = 60
        elif ROI == True:
            up = 10
            down = 30
            n = 20

        self._breakout_indices = list()

        center_breakout_index = find_breakout_index(self._proc_fsv
                                                    [:,int(xdim/2)])
        self._breakout_indices.append(center_breakout_index)

        # Iterating from the center to the left
        for i in arange(int(xdim/2)-1,-1,-1): 
            
            # Find the breakout index
            breakout_index = find_breakout_index(self._proc_fsv[:,i]) 

            # Test if the breakout index is much smaller than the first entry
            if (self._breakout_indices[0] - breakout_index) > up: 
                breakout_index2 = find_breakout_index(self._proc_fsv
                                                      [breakout_index+1:,i])
                breakout_index += breakout_index2
                
                # Test if it is within n of the first entry
                if abs(breakout_index - self._breakout_indices[0]) > n: 
                    
                    # Repeat the first entry
                    self._breakout_indices.insert(0,self._breakout_indices[0]) 
                else:
                    self._breakout_indices.insert(0,breakout_index)

            # Test if the breakout index is much larger than the first entry
            elif (breakout_index - self._breakout_indices[0]) > down: 
                breakout_index2 = find_breakout_index(self._proc_fsv
                                                      [:breakout_index,i])
                breakout_index = breakout_index2
                
                # Test if it is within n of the first entry
                if abs(breakout_index - self._breakout_indices[0]) > n: 
                    
                    # Repeat the first entry
                    self._breakout_indices.insert(0,self._breakout_indices[0]) 
                else:
                    self._breakout_indices.insert(0,breakout_index)

            else:
                self._breakout_indices.insert(0,breakout_index)

        # Find the breakout index
        for i in arange(int(xdim/2)+1,self._proc_fsv.shape[1]):
            breakout_index = find_breakout_index(self._proc_fsv[:,i]) 

            # Test if the breakout index is much smaller than the last entry
            if (self._breakout_indices[-1] - breakout_index) > up: 
                breakout_index2 = find_breakout_index(self._proc_fsv
                                                      [breakout_index+1:,i])
                breakout_index += breakout_index2
                
                # Test if it is within n of the last entry
                if abs(breakout_index - self._breakout_indices[-1]) > n: 
                    
                    # Repeat the last entry
                    self._breakout_indices.append(self._breakout_indices[-1]) 
                else:
                    self._breakout_indices.append(breakout_index)

            # Test if the breakout index is much larger than the last entry
            elif (breakout_index - self._breakout_indices[-1]) > down: 
                breakout_index2 = find_breakout_index(self._proc_fsv
                                                      [:breakout_index,i])
                breakout_index = breakout_index2
                
                # Test if it is within n of the last entry
                if abs(breakout_index - self._breakout_indices[-1]) > n: 
                    
                    # Repeat the last entry
                    self._breakout_indices.append(self._breakout_indices[-1]) 
                else:
                    self._breakout_indices.append(breakout_index)

            else:
                self._breakout_indices.append(breakout_index)

        self._breakout_indices = asarray(self._breakout_indices)
                
    def shift_to_zero(self):
        '''Performs background subraction of average of the first 50 points s.t
           FSV starts at ~0.
        '''
        
        for i in arange(self._proc_fsv.shape[1]):
            fsv_background = 0
            
            for j in arange(50):
                fsv_background += self._proc_fsv[j,i]
                
            # Find average of first 50 points
            fsv_background /= 50
            
            # Perform subtraction
            self._proc_fsv[:,i] -= fsv_background
       
            
class MEC_VISAR_image_bed1(VISAR_image):
    '''Class that implements a visar image of bed 1 of the SLAC MEC VISAR.

       Uses time calibration of bed1 of the MEC visar. Otherwise identical to 
       parent class VISAR_image
    '''

    def __init__(self,filename=None):
        VISAR_image.__init__(self,filename=filename)
        ### definition of the time calibration coefficient

        self._time_coeff['0.5 ns']=[5.44528e-04,-1.85189e-07,8.31059e-11]
        self._time_coeff['0.5 ns'].reverse()

        self._time_coeff['1 ns']=[1.17096e-03,-4.73566e-07,2.96710e-10]
        self._time_coeff['1 ns'].reverse()

        self._time_coeff['2 ns']=[1.74166e-03,-4.02184e-07,2.09947e-10]
        self._time_coeff['2 ns'].reverse()

        self._time_coeff['5 ns']=[5.74485e-03,-1.97301e-06,1.39066e-09]
        self._time_coeff['5 ns'].reverse()

        self._time_coeff['10 ns']=[9.84870e-03,-3.53869e-06,2.66292e-09]
        self._time_coeff['10 ns'].reverse()

        self._time_coeff['20 ns']=[1.91357e-02,-4.75469e-06,6.47073e-09]
        self._time_coeff['20 ns'].reverse()

        self._time_coeff['50 ns']=[4.97666e-02,-1.37149e-05,7.42748e-09]
        self._time_coeff['50 ns'].reverse()

        self._time_coeff['100 ns']=[1.01233e-01,-2.69505e-05,1.01755e-08]
        self._time_coeff['100 ns'].reverse()

        self._time_coeff['200 ns']=[1.87951e-01,-3.00835e-05,1.51207e-08]
        self._time_coeff['200 ns'].reverse()

        self._time_coeff['500 ns']=[4.57900e-01,-2.16852e-05,2.43540e-08]
        self._time_coeff['500 ns'].reverse()
                
                
class MEC_VISAR_image_bed2(VISAR_image):
    '''Class that implements a visar image of bed 2 of the SLAC MEC VISAR.

       Uses time calibration of bed1 of the MEC visar. Otherwise identical to 
       parent class VISAR_image
    '''

    def __init__(self,filename=None):
        VISAR_image.__init__(self,filename=filename)

        self._time_coeff['0.5 ns']=[5.02968e-04,-1.26507e-07,7.05022e-11]
        self._time_coeff['0.5 ns'].reverse()

        self._time_coeff['1 ns']=[1.22612e-03,-5.72856e-07,3.33411e-10]
        self._time_coeff['1 ns'].reverse()

        self._time_coeff['2 ns']=[1.91378e-03,9.13826e-10,-3.29732e-10]
        self._time_coeff['2 ns'].reverse()

        self._time_coeff['5 ns']=[5.60041e-03,1.30464e-06,-2.36493e-09]
        self._time_coeff['5 ns'].reverse()

        self._time_coeff['10 ns']=[9.04307e-03,4.10175e-06,-3.82061e-09]
        self._time_coeff['10 ns'].reverse()

        self._time_coeff['20 ns']=[1.52475e-02,9.14541e-06,-5.96076e-09]
        self._time_coeff['20 ns'].reverse()

        self._time_coeff['50 ns']=[5.36197e-02,-1.65190e-05,9.97896e-09]
        self._time_coeff['50 ns'].reverse()

        self._time_coeff['100 ns']=[9.31495e-02,4.38192e-06,3.60153e-10]
        self._time_coeff['100 ns'].reverse()

        self._time_coeff['200 ns']=[1.92754e-01,2.08386e-05,8.00533e-10]
        self._time_coeff['200 ns'].reverse()

        self._time_coeff['500 ns']=[4.59835e-01,9.89071e-05,-4.39776e-08]
        self._time_coeff['500 ns'].reverse()