#!/usr/bin/env python
'''
Name:          visar_ROI_selection.py
License:       MIT
Version:       1.0
Last modified: 8 Aug. 2018 (SBB)
Authors:       Akel Hashim          (ahashim@slac.stanford.edu)
               Shaughnessy Brown    (sbbrown@slac.stanford.edu)
               Bob Nagler           (bnagler@slac.stanford.edu)
Description:   Creates ROI object that specifies rectagular selection on image 
               from user mouse actions
'''

'''Imports -----------------------------------------------------------------'''
import matplotlib.pyplot as plt

'''Class Definitions -------------------------------------------------------'''
class ROICoordinates:
    def __init__(self):
        self.x0 = None
        self.y0 = None
        self.xf = None
        self.yf = None
     
    def connect(self):
        '''Connects events.'''
        
        plt.connect('button_press_event', self.on_click)
        plt.connect('button_release_event', self.on_release)        
        
    def on_click(self,event):
        '''Records x0 and y0 coords when clicked'''
        
        self.x0, self.y0 = event.xdata, event.ydata
        if event.button == 1:
            if event.inaxes is not None:
                print('\nData coordinates: x0 = %f, y0 = %f'                  \
                      % (event.xdata, event.ydata))
                                          
    def on_release(self,event):
        '''Gets xf and yf coordinates when released'''
        
        self.xf, self.yf = event.xdata, event.ydata
        if event.button == 1:
            if event.inaxes is not None:
                print('\nData coordinates: xf = %f, yf = %f \n'               \
                      % (event.xdata, event.ydata))
