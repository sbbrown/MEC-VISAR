#!/usr/bin/env python
'''
Name:          visar_plotting_tools.py
License:       MIT
Version:       2.0
Last modified: 30 Jan. 2019 (SBB)
Authors:       Akel Hashim          (ahashim@slac.stanford.edu)
               Bob Nagler           (bnagler@slac.stanford.edu)
               Shaughnessy Brown    (sbbrown@slac.stanford.edu)
Description:   Plotting fuctions that display graphs user requested in GUI
'''

'''Imports -----------------------------------------------------------------'''
import time

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl

import numpy as np

try:    # For Python2
    import PIL
except: # For Python3
    try:
        import pillow
    except:
        print('Import Error: PIL/pillow is needed to run the analysis')
        sys.exit(1)

import threading
from textwrap import wrap

from visar_image_analysis import *
#from mecvisarimage import *
from visar_ROI_selection import ROICoordinates


'''Global Functions --------------------------------------------------------'''

def determine_xticks_labels(FOV,xdim):
    '''Determines xticks and xtick labels for plotting'''
    
    xticks = list(np.linspace(0,xdim,6))
    xtick_labels = [str(int(FOV * x/xdim)) for x in xticks]
    
    return (xticks,xtick_labels)

def determine_yticks_labels(time_base,ydim):
    '''Determines yticks and ytick labels for plotting'''
    
    for i,item in enumerate(time_base):
        if time_base[i].isdigit():
            continue
        else:
            stop = i
            break
    
    #print(time_base)
    #print(type(time_base))
    sweep_window = float(time_base[:stop])
    
    yticks = list(np.linspace(0,ydim-1,6))
    ytick_labels = [str(int(sweep_window * y/(ydim-1))) for y in yticks]
    
    return (yticks,ytick_labels,sweep_window)
        
def ROI_ticks_labels(time_base,FOV,ROI_info,xdim,ydim):
    '''Determines ticks and labels for plotting ROI'''
    
    for i,item in enumerate(time_base):
        if time_base[i].isdigit():
            continue
        else:
            stop = i
            break
    
    #print(time_base)
    #print(type(time_base))
    sweep_window = float(time_base[:stop])
    
    # Assign locations selected by ROI script
    x0 = ROI_info[0]
    xf = ROI_info[1]
    y0 = ROI_info[2]
    yf = ROI_info[3]

    ROI_xticks_pixels = list(np.linspace(0,xdim-1,6))
    ROI_xticks_labels = [str(int(i)) for i in np.linspace(FOV*x0/1344,
                         FOV*xf/1344,6)]

    ROI_yticks_pixels = list(np.linspace(0,ydim-1,4))
    ROI_yticks_labels = [str(round(i,1)) for i in np.linspace(
                         sweep_window*y0/1024,sweep_window*yf/1024,4)]
    
    return (ROI_xticks_pixels, ROI_xticks_labels, 
            ROI_yticks_pixels, ROI_yticks_labels)
    
ROI_pixels=tuple(8*[None])
ROI_BT_pixels=tuple(8*[None])

def select_ROI_plots(img1,img2):
    '''Plots the raw fringe images for the ROI selection and returns 
       coordinates selected by click and release.'''
    global ROI_pixels
    
    print('\nWARNING: If Interactive Plotting selected prior, ROI selection ' \
          'plots may not function properly. If this occurs, restart mecana.\n')
          
    title = 'Select the ROI using Zoom tool, dragging from top left to       '\
            'bottom right of desired window. Press HOME to view the entire   '\
            'region. Close both windows when finished.'

    # Display data image 1
    fig1, ax1 = plt.subplots(figsize=(9,7))
    fig1.suptitle('\n'.join(wrap('Bed 1: ' + title,70)), fontsize = 15)
    imgplot = ax1.imshow(img1, cmap='Greys_r')
    
    # Select ROI of data image 2
    ROI_V1 = ROICoordinates()
    ROI_V1.connect()
    ROI1_done = threading.Event()
    fig1.canvas.mpl_connect('close_event', lambda evt: ROI1_done.set())
    fig1.show()
    if ROI_pixels != None:
        fig1.canvas.toolbar.push_current()
        ax1.set_xlim(ROI_pixels[0], ROI_pixels[1])
        ax1.set_ylim(ROI_pixels[2], ROI_pixels[3])
    
    # Display data image 1
    fig2, ax2 = plt.subplots(figsize=(9,7))
    fig2.suptitle('\n'.join(wrap('Bed 2: ' + title,70)), fontsize = 15)
    imgplot = ax2.imshow(img2, cmap='Greys_r')
    
    # Select ROI of data image 2
    ROI_V2 = ROICoordinates()
    ROI_V2.connect()
    ROI2_done = threading.Event()
    fig2.canvas.mpl_connect('close_event', lambda evt: ROI2_done.set())
    fig2.show()
    if ROI_pixels != None:
        fig2.canvas.toolbar.push_current()
        ax2.set_xlim(ROI_pixels[4], ROI_pixels[5])
        ax2.set_ylim(ROI_pixels[6], ROI_pixels[7])

    while (not ROI1_done.isSet() or not ROI2_done.isSet()):
        try:
            plt.pause(0.5)
        except:
            pass

    ROI_pixels = (int(ROI_V1.x0) if ROI_V1.x0 != None else ROI_pixels[0],
                  int(ROI_V1.xf) if ROI_V1.xf != None else ROI_pixels[1],
                  int(ROI_V1.y0) if ROI_V1.y0 != None else ROI_pixels[2],
                  int(ROI_V1.yf) if ROI_V1.yf != None else ROI_pixels[3],
                  int(ROI_V2.x0) if ROI_V2.x0 != None else ROI_pixels[4],
                  int(ROI_V2.xf) if ROI_V2.xf != None else ROI_pixels[5],
                  int(ROI_V2.y0) if ROI_V2.y0 != None else ROI_pixels[6],
                  int(ROI_V2.yf) if ROI_V2.yf != None else ROI_pixels[7])
    
    return ROI_pixels
    
def plot_raw(img1,img2,img_FOV,img_time_base):
    '''Plots the raw fringe images'''
    
    xdim = img1.shape[1]
    ydim = img1.shape[0]
    
    xtick_info = determine_xticks_labels(img_FOV,xdim)
    ytick_info = determine_yticks_labels(img_time_base,ydim)
    
    intensity_lineout = ytick_info[1][1]
    intensity_label = 'S(f,t)|t = ' + intensity_lineout + ' ns'
    
    mpl.rc('font',family='Times New Roman')
    
    fig = plt.figure(figsize = (15,10))
    fig.canvas.set_window_title('Raw Fringe Files & Intensity Line-Out')
    
    gs = gridspec.GridSpec(1, 2, height_ratios=[1]) 
    # gs = gridspec.GridSpec(1, 2, height_ratios=[2, 1]) 
    
    a1 = fig.add_subplot(gs[0])
    a1.set_title('VISAR Bed 1', fontsize = 20)
    plt.xlabel('x ($\mu$m)', fontsize = 20)
    plt.xticks(xtick_info[0],xtick_info[1], fontsize = 20)
    plt.ylabel('Time (ns)', fontsize = 20)
    plt.yticks(ytick_info[0],ytick_info[1], fontsize = 20)
    imgshow1 = plt.imshow(log10(img1), cmap = 'Greys', interpolation='none')


    a2 = fig.add_subplot(gs[1])
    a2.set_title('VISAR Bed 2', fontsize = 20)
    plt.xlabel('x ($\mu$m)', fontsize = 20)
    plt.xticks(xtick_info[0],xtick_info[1], fontsize = 20)
    plt.ylabel('Time (ns)', fontsize = 20)
    plt.yticks(ytick_info[0],ytick_info[1], fontsize = 20)
    imgshow2 = plt.imshow(log10(img2), cmap = 'Greys', interpolation='none')

    plt.tight_layout()

    plt.show(block=False)
    
def plot_fsv_lineout(im1,im2, plot=True, save=False, ROI_info=None, 
                     bed1=None, bed2=None):
    '''Plots two FSVs such that their breakout times overlap. Saves the FSV 
       lineouts to csv if selected in GUI'''

    xdim1 = im1._proc_im.shape[1]
    ydim1 = im1._proc_im.shape[0]
    xdim2 = im2._proc_im.shape[1]
    ydim2 = im2._proc_im.shape[0]
    
    if ydim1 < 1024:
        ROI = True        
    else:
        ROI = False

    t0_v1 = im1._t0
    t0_v2 = im2._t0

    # Build the two 1D FSV arrays with proper labels (ROI vs. no ROI)
    if bed1:
        y1 = []
        for v in bed1:
            im1._proc_fsv = v;
            y1.append(im1.fsv_lp(show_plot=False,return_array=True))
        pass
    else:
        y1 = [im1.fsv_lp(show_plot=False,return_array=True)]
    if bed2:
        y2 = []
        for v in bed2:
            im2._proc_fsv = v;
            y2.append(im2.fsv_lp(show_plot=False,return_array=True))
        pass
    else:
        y2 = [im2.fsv_lp(show_plot=False,return_array=True)]

    if ROI == False:  
        xtick_info = determine_xticks_labels(im1._FOV,xdim1)
        ytick_info = determine_yticks_labels(im1._time_base,ydim1)
        time_window = ytick_info[2]
        
        x1 = linspace(0,time_window,ydim1) - t0_v1
        x2 = linspace(0,time_window,ydim2) - t0_v2
        
    elif ROI == True:        
        ROI_tick_info_V1 = ROI_ticks_labels(im1._time_base,im1._FOV,
                                            ROI_info[:4],xdim1,ydim1)
        ROI_tick_info_V2 = ROI_ticks_labels(im2._time_base,im2._FOV,
                                            ROI_info[4:],xdim2,ydim2)
        
        x1 = linspace(float(ROI_tick_info_V1[3][0]),
                            float(ROI_tick_info_V1[3][-1]),ydim1) - t0_v1
        x2 = linspace(float(ROI_tick_info_V2[3][0]),
                            float(ROI_tick_info_V2[3][-1]),ydim2) - t0_v2
     
    # Specify plotting options
    mpl.rc('font',family='Times New Roman')
    
    fig = plt.figure(figsize = (15,8))
    fig.canvas.set_window_title('Free Surface Velocity')
    
    a = fig.add_subplot(111)
    a.set_title('Free Surface Velocity Line-out', fontsize = 20)
    plt.xlabel('Time (ns)', fontsize = 20)
    plt.xticks(fontsize = 20)
    plt.ylabel('Free Surface Velocity (m/s)', fontsize = 20)
    plt.yticks(fontsize = 20)
    
    # Add FSVs to plot
    for y in y1:
        bed_1, = plt.plot(x1, y, 'r-', label = 'Bed 1')
    for y in y2:
        bed_2, = plt.plot(x2, y, 'b-', label = 'Bed 2')
    
        plt.legend(['bed 1', 'bed 2'], prop={'size':20}, loc='upper left')
    ax = plt.gca()
    leg = ax.get_legend()
    leg.legendHandles[0].set_color('red')
    leg.legendHandles[1].set_color('blue')
    
    
    
    # Save to .csv if checked in GUI
    if save == True:
        z1 = np.c_[x1,y1[0]]
        for i in range(1, len(y1)):
            z1 = np.c_[z1,y1[i]]
        z2 = np.c_[x2,y2[0]]
        for i in range(1, len(y2)):
            z2 = np.c_[z2,y2[i]]
        
        np.savetxt('FSV1.csv', z1, delimiter=",") # Save in the local directory
        np.savetxt('FSV2.csv', z2, delimiter=",") # Save in the local directory
        
        print('FSV line-outs saved to FSV1.csv and FSV2.csv in local dir. \n')
        
    if plot == True:
        plt.tight_layout()
     
        plt.show()
