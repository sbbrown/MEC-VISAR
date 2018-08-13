#!/usr/bin/env python
'''
Name:          visar_gui.py
License:       MIT
Version:       1.0
Last modified: 7 Aug. 2018 (SBB)
Authors:       Akel Hashim          (ahashim@slac.stanford.edu)
               Bob Nagler           (bnagler@slac.stanford.edu)
               Shaughnessy Brown    (sbbrown@slac.stanford.edu)
Description:   Creates GUI for user to specify VISAR files to be analyzed, 
               key interferometry metrics (i.e. etalons, field-of-view), and 
               desired output plotting. Script loads images into arrays, passes
               arrays and specified interferometry metrics to 
               launch_visar_analysis
'''

'''Imports -----------------------------------------------------------------'''
import sys
import os

try:    # For Python2
    import Tkinter as tk
    import ttk
except: # For Python3
    try:
        import tkinter as tk
         
    except:
        print "Import Error: Tkinter is needed to run the GUI"
        sys.exit(1)
     
from Tkinter import Tk
from tkFileDialog import askopenfilename
        
from launch_visar_analysis import LaunchVISAR

'''Class Definitions -------------------------------------------------------''' 

class VisarGUI:
    '''Interactive GUI for VISAR Analysis. GUI is organized in 6 frames.
        1. Canvas, body frame, and scrollbar
        2. Selecting files to analyze
            Reference images
            Data images
        3. Querying required analysis parameters
            Etalon thicknesses
            Field of view
        4. Querying optional analysis parameters
            Select ROIs
            Direction of FSV calibration
            FSV direction, Bed 1
            FSV direction, Bed 2
            Fringe jump corrections
            Fringe jumps, Bed 1
            Fringe jumps, Bed 2
        5. Specifying plotting options
            Y/N Enable interactive
            # Y/N View raw images
            # Y/N View final FSV lineout
            # Y/N Save final FSV lineout
        6. Analyze & Quit buttons
    '''
   
    def __init__(self,master):
        self.root = master
        self.home = os.path.expanduser('~/Desktop')
        master.title('VISAR Analysis')
        master.resizable(True, True)
        master.wm_geometry("800x850")
        
        # First frame for canvas, body frame, and scrollbar ----------------- 1
        self.canvas = tk.Canvas(master)
        self.canvas.pack(side='left', expand=1, fill='both')
        
        self.frame_body = tk.Frame(self.canvas)
        self.frame_body.grid(row=0,column=0)
        
        self.sb = tk.Scrollbar(master, orient='vertical', 
                               command=self.canvas.yview)
        self.sb.pack(side='right', fill='y')
        
        self.canvas.configure(yscrollcommand=self.sb.set)
        self.canvas.create_window((0,0),window=self.frame_body,anchor='nw')
        self.frame_body.bind('<Configure>', self.ScrollFunction)
        
        # Second frame for selecting files to analyze ----------------------- 2
        self.frame_content = tk.Frame(self.frame_body)
        self.frame_content.grid(row=0, column=0, sticky='news') 
        ttk.Label(self.frame_content, text='Input Files:', 
                  font=('Arial', 18, 'bold', 'underline')).                   \
                  grid(row=0, column=0, padx=5, pady=5, sticky='w') 
        
        ### Reference images -------------------------
        ttk.Label(self.frame_content, text='Reference Images:', 
                  font=('Arial', 12, 'bold')).                                \
                  grid(row=1, column=0, padx=5, sticky='w')
        
        # Reference image, Bed 1
        ttk.Label(self.frame_content, text='Bed1:', 
                  font=('Arial', 12)).grid(row=1, column=1, padx=5, sticky='w')
        self.file_bed1_ref = tk.Button(self.frame_content, 
                                       text='select bed1 ref', width=15,
                                       justify='center', fg='blue', 
                                       font=('Arial', 12), 
                                       command = self.select_bed1_ref)          
        
        self.file_bed1_ref.grid(row=1, column=2, padx=5, sticky='w')           
        
        # Reference image, Bed 2
        ttk.Label(self.frame_content, text='Bed2:',
                  font=('Arial', 12)).grid(row=1, column=5, padx=5, sticky='w')
        self.file_bed2_ref = tk.Button(self.frame_content, 
                                       text='select bed2 ref', width=15, 
                                       justify='center', fg='blue', 
                                       font=('Arial', 12), 
                                       command = self.select_bed2_ref)          
        self.file_bed2_ref.grid(row=1, column=6, padx=5, sticky='w')
   
        ### Data images -------------------------
        ttk.Label(self.frame_content, text='Data Images:',
                  font=('Arial', 12, 'bold')).                                \
                  grid(row=3, column=0, padx=5,  sticky='w')

        # Data image, Bed 1
        ttk.Label(self.frame_content, text='Bed1:', 
                  font=('Arial', 12)).grid(row=3, column=1, padx=5, sticky='w')
        self.file_bed1 = tk.Button(self.frame_content, 
                                   text='select bed1', width=15, 
                                   justify='center', fg='blue', 
                                   font=('Arial', 12), 
                                   command = self.select_bed1_data)             
        self.file_bed1.grid(row=3, column=2, padx=5, sticky='w')              
        
        # Data image, Bed 2
        ttk.Label(self.frame_content, text='Bed2:', 
                  font=('Arial', 12)).grid(row=3, column=5, padx=5, sticky='w')
        self.file_bed2 = tk.Button(self.frame_content, 
                                   text='select bed2', width=15, 
                                   justify='center', fg='blue', 
                                   font=('Arial', 12), 
                                   command = self.select_bed2_data)             
        self.file_bed2.grid(row=3, column=6, padx=5, sticky='w')              
        
        # Third frame for querying required analysis parameters ------------- 3
        self.frame_parameters = tk.Frame(self.frame_body)
        self.frame_parameters.grid(row=1, column=0, sticky='news')
        ttk.Label(self.frame_parameters, text='Input Parameters:', 
                  font=('Arial', 18, 'bold', 'underline')).                   \
                  grid(row=0, column=0, padx=5, pady=5, sticky='w') 
        
        ### Etalon thicknesses -------------------------
        ttk.Label(self.frame_parameters, text='Etalon Thicknesses:', 
                  font=('Arial', 12, 'bold')).                                \
                  grid(row=1, column=0, padx=5, sticky='w')
        
        # Etalon thickness, Bed 1
        ttk.Label(self.frame_parameters, text='Bed1:', 
                  font=('Arial', 12)).grid(row=1, column=1, padx=5, sticky='w')
        ttk.Label(self.frame_parameters, text='mm', 
                  font=('Arial', 12)).grid(row=1, column=3, sticky='w')
        self.entry_h_v1 = ttk.Entry(self.frame_parameters, width=12, 
                                    justify='center', font = ('Arial', 12))     
        self.entry_h_v1.grid(row=1, column=2, padx=5, sticky='w')
        
        # Etalon thickness, Bed 2
        ttk.Label(self.frame_parameters, text='Bed2:', 
                  font=('Arial', 12)).grid(row=1, column=5, padx=5, sticky='w')
        ttk.Label(self.frame_parameters, text='mm', 
                  font=('Arial', 12)).grid(row=1, column=7, sticky='w')
        self.entry_h_v2 = ttk.Entry(self.frame_parameters, width=12, 
                                    justify='center', font = ('Arial', 12))
        self.entry_h_v2.grid(row=1, column=6, padx=5, sticky='w')
        
        ### Field of view ------------------------------
        ttk.Label(self.frame_parameters, text='Field of View (FOV):', 
                  font=('Arial', 12, 'bold')).                                \
                  grid(row=2, column=0, padx=5, sticky='w')
        ttk.Label(self.frame_parameters, text='microns', 
                  font=('Arial', 12)).grid(row=2, column=3, sticky='w')
        self.entry_FOV = ttk.Entry(self.frame_parameters, width=12, 
                                   justify='center', font = ('Arial', 12))
        self.entry_FOV.grid(row=2, column=2, padx=5, sticky='w')
        
        # Fourth frame for querying optional analysis parameters ------------ 4
        #### Note: Slave to root window, not frame_options
        self.frame_options = tk.Frame(self.frame_body)
        self.frame_options.grid(row=3, column=0, sticky='news') 
        ttk.Label(self.frame_body, text='Optional Analysis Tools:', 
                  font=('Arial', 18, 'bold', 'underline')).                   \
                  grid(row=2, column=0, padx=5, pady=5, sticky='w')
      
        ### Select ROIs --------------------------------
        self.frame_ROI = tk.Frame(self.frame_options)
        self.frame_ROI.grid(row=1,column=0, sticky='news')
        
        # Y/N Implement ROI
        self.ROI = tk.IntVar()
        self.chk_ROI = tk.Checkbutton(self.frame_ROI, variable=self.ROI, 
                                      text='Implement ROI for entire image', 
                                      font = ('Arial', 12, 'bold'))            
        self.chk_ROI.grid(row=1, column=0, padx=5, sticky='w')
        
        ### Direction of FSV calibration ---------------
        self.frame_FSV_calib = tk.Frame(self.frame_options)
        self.frame_FSV_calib.grid(row=2, column=0, sticky='news')
        
        # Y/N Specify FSV direction
        self.FSV_calib = tk.IntVar()
        self.chk_FSVcalib = tk.Checkbutton(self.frame_FSV_calib, 
                                           variable=self.FSV_calib, 
                                           text='FSV calibration direction:', 
                                           font = ('Arial', 12, 'bold'), 
                                           command=lambda: 
                                               self.check_FSV_calib())
        self.chk_FSVcalib.grid(row=0, column=0, padx=5, sticky='w')

        # FSV direction, Bed 1
        ttk.Label(self.frame_FSV_calib, text='Bed1:', 
                  font=('Arial', 12)).grid(row=0, column=1, padx=8, sticky='w')
        ttk.Label(self.frame_FSV_calib, text='(def. = -1)', 
                  font=('Arial', 12)).grid(row=0, column=3, padx=5, sticky='w')
        self.entry_FSVdir_v1 = ttk.Entry(self.frame_FSV_calib, 
                                         state='disabled', width=6, 
                                         justify='center', 
                                         font = ('Arial', 12))                  
        self.entry_FSVdir_v1.grid(row=0, column=2, sticky='e')
        
        # FSV direction, Bed 2
        ttk.Label(self.frame_FSV_calib, text='Bed2:', 
                  font=('Arial', 12)).grid(row=0, column=4, padx=8, sticky='w')
        ttk.Label(self.frame_FSV_calib, text='(def. = 1)', 
                  font=('Arial', 12)).grid(row=0, column=6, padx=5, sticky='w')
        self.entry_FSVdir_v2 = ttk.Entry(self.frame_FSV_calib, 
                                         state='disabled', width = 6, 
                                         justify='center', 
                                         font = ('Arial', 12))
        self.entry_FSVdir_v2.grid(row=0, column=5, sticky='e')
        
        ### Fringe jump corrections --------------------
        self.FJC = tk.IntVar()
        self.chk_FJC = tk.Checkbutton(self.frame_options, variable=self.FJC, 
                                      text='Manually specify fringe jump '
                                      'corrections for each bed. Inputs are '
                                      'integers or integer ranges with ":"', 
                                      wraplength=800, justify='left', 
                                      font=('Arial', 12, 'bold'), 
                                      command=lambda: self.check_FJC())
        self.chk_FJC.grid(row=4, column=0, padx=5, sticky='w')
        
        self.frame_fjc = tk.Frame(self.frame_options)
        self.frame_fjc.grid(row=5, column=0, sticky='news')
        
        ttk.Label(self.frame_fjc, text='Note: If left unchecked, the analysis '
                  'routine will automatically attempt to find a match for the '
                  'two beds.', wraplength=800, justify='left',
                  font=('Arial',12)).                                         \
                  grid(row=0, column=0, columnspan=12, padx=30, sticky='w')  

        ttk.Label(self.frame_fjc, text='If either input is a range, no data '
                  'will be saved and the final FSV lineout will be shown.',
                  wraplength=800, justify='left', font=('Arial',12)).         \
                  grid(row=1, column=0, columnspan=12, padx=30, sticky='w')  
        
        # Fringe jumps, Bed 1
        ttk.Label(self.frame_fjc, text='Bed1:', font=('Arial', 12)).          \
        grid(row=2, column=0, padx=30, sticky='w')                            
        ttk.Label(self.frame_fjc, text='m = ', font=('Arial', 12)).           \
        grid(row=2, column=0, sticky='e')                                     
        self.entry_m = ttk.Entry(self.frame_fjc, state='disabled', width=6, 
                                 justify='center', font = ('Arial', 12))
        self.entry_m.grid(row=2, column=1, sticky='w')
        
        # Fringe jumps, Bed 2
        ttk.Label(self.frame_fjc, text='Bed2:', font=('Arial', 12)).          \
        grid(row=2, column=4, padx=30, sticky='w')
        ttk.Label(self.frame_fjc, text='n = ', font=('Arial', 12)).           \
        grid(row=2, column=4, sticky='e')
        self.entry_n = ttk.Entry(self.frame_fjc, state='disabled', width=6, 
                                 justify='center', font = ('Arial', 12))        
        self.entry_n.grid(row=2, column=5, sticky='w')
        
        # Fifth frame for specifying plotting options ----------------------- 5
        self.frame_plotting = tk.Frame(self.frame_options)
        self.frame_plotting.grid(row=6, column=0, sticky='news')
        ttk.Label(self.frame_plotting, text='Plotting options:', 
                  font=('Arial', 12, 'bold')).                                \
                  grid(row=0, column=0, padx=5, sticky='w')
       
        # View raw images
        self.plot_RF = tk.IntVar()
        self.chk_RF = tk.Checkbutton(self.frame_plotting, 
                                     variable=self.plot_RF, 
                                     text='View raw images of both beds.', 
                                     wraplength=800, justify='left', 
                                     font=('Arial', 12))
        self.chk_RF.grid(row=2, column=0, padx=30, sticky='w')                
        
        # View final free-surface velocity lineout
        self.plot_FSV = tk.IntVar()
        self.chk_FinalFSV = tk.Checkbutton(self.frame_plotting, 
                                           variable=self.plot_FSV, 
                                           text='View the final FSV line-out.', 
                                           wraplength=800, justify='left', 
                                           font=('Arial', 12))
        self.chk_FinalFSV.grid(row=8, column=0, padx=30, sticky='w')
        
        # Save final free-surface velocity lineout
        self.save_FSV = tk.IntVar()
        self.chk_saveFSV = tk.Checkbutton(self.frame_plotting, 
                                          variable=self.save_FSV, 
                                          text='Save the final FSV line-out', 
                                          wraplength=800, justify='left', 
                                          font=('Arial', 12))                   
        self.chk_saveFSV.grid(row=11, column=0, padx=30, sticky='w')
        
        # Analyze & Quit buttons--------------------------------------------- 6
        self.frame_buttons = tk.Frame(self.frame_plotting)
        self.frame_buttons.grid(row=0, column=1, rowspan=12, sticky='news')
        self.frame_buttons.grid_columnconfigure(0, weight=1)
        self.frame_buttons.grid_rowconfigure(0, weight=1)
        
        self.analyze_button = tk.Button(self.frame_buttons, text="Analyze", 
                                        width=10, justify='center', fg='blue', 
                                        font=('Arial', 24, 'bold'), 
                                        command = self.Analyze)                 
        self.analyze_button.grid(row=0, column=0, padx=100, pady=10,sticky='s')
        
        self.quit_button = tk.Button(self.frame_buttons, text="Quit", 
                                     width=10, justify='center', fg='red', 
                                     font=('Arial', 24, 'bold'), 
                                     command = self.Quit)
        self.quit_button.grid(row=1, column=0, padx=100, pady=20, stick='n')    
          
    def ScrollFunction(self,event):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"), height=800)    
               
    def Analyze(self):
        '''Assigns variables queried from GUI and passes to analysis script 
        via instance of class LaunchVISAR'''
        
        print "\nAnalyzing the VISAR data now...\n"
        
        # Set variables to pass to launch_visar script
        f_v1_ref = self.filename_bed1_ref
        f_v2_ref = self.filename_bed2_ref
        f_v1     = self.filename_bed1
        f_v2     = self.filename_bed2
        
        h_v1 = self.entry_h_v1.get()
        h_v2 = self.entry_h_v2.get()
        
        FOV = self.entry_FOV.get()
        ROI = self.ROI.get()
        
        FSV_calib = self.FSV_calib.get()
        FSVdir_v1 = self.entry_FSVdir_v1.get()
        FSVdir_v2 = self.entry_FSVdir_v2.get()
        
        FJC = self.FJC.get()
        m   = self.entry_m.get()
        n   = self.entry_n.get()
        
#        iplot   = self.iplot.get()
        plotRF  = self.plot_RF.get()
        plotFSV = self.plot_FSV.get()
        saveFSV = self.save_FSV.get()
        
        print '\nBeginning VISAR analysis now...\n'
        Analyze = LaunchVISAR(f_v1_ref, f_v2_ref, f_v1, f_v2, h_v1, h_v2, FOV, 
                              ROI, FSV_calib, FSVdir_v1, FSVdir_v2, FJC, m, n, 
                              plotRF, plotFSV, saveFSV)
        Analyze()
        print 'The VISAR analysis has finished for the selected runs!'
          
    def Quit(self):
        
        print '\nQuitting the application! \n'
        self.root.destroy()
        sys.exit(0)
        
    def select_bed1_ref(self):
        '''Selects the reference VISAR image of bed1'''
        
        Tk().withdraw()
        self.filename_bed1_ref = askopenfilename(title = 'Select reference'
                                                     ' data for bed 1:', 
                                                     initialdir = self.home)
        print 'Selected file ' + self.filename_bed1_ref + ' for bed 1.\n'
        ttk.Label(self.frame_content, text=self.filename_bed1_ref, 
                  font=('Arial', 8)).grid(row=2, column=1, padx=5, sticky='w')
        
    def select_bed2_ref(self):
        '''Selects the reference VISAR image of bed2'''
        
        Tk().withdraw()
        self.filename_bed2_ref = askopenfilename(title = 'Select reference'
                                                     ' data for bed 2:',
                                                     initialdir = self.home)
        print 'Selected file ' + self.filename_bed1_ref + ' for bed 2.\n'
        ttk.Label(self.frame_content, text=self.filename_bed2_ref, 
                  font=('Arial', 8)).grid(row=2, column=5, padx=5, sticky='w')

    def select_bed1_data(self):
        '''Selects the VISAR image of bed1'''
        
        Tk().withdraw()
        self.filename_bed1 = askopenfilename(title = 'Select shot data for'
                                                 ' bed 1:',
                                                 initialdir = self.home)
        print 'Selected file ' + self.filename_bed1 + ' for bed 1.\n'
        ttk.Label(self.frame_content, text=self.filename_bed1, 
                  font=('Arial', 8)).grid(row=4, column=1, padx=5, sticky='w')
        
    def select_bed2_data(self):
        '''Selects the VISAR image of bed2'''
        
        Tk().withdraw()
        self.filename_bed2 = askopenfilename(title = 'Select shot data for'
                                                 ' bed 2:',
                                                 initialdir = self.home)
        print 'Selected file ' + self.filename_bed2 + ' for bed 2.\n'        
        ttk.Label(self.frame_content, text=self.filename_bed2, 
                  font=('Arial', 8)).grid(row=4, column=5, padx=5, sticky='w')
            
    def check_FSV_calib(self):
        '''Checks if free surface velocity direction is selected'''
        
        if self.FSV_calib.get() == 0:
            self.entry_FSVdir_v1.configure(state='disabled')
            self.entry_FSVdir_v2.configure(state='disabled')
        elif self.FSV_calib.get() == 1:
            self.entry_FSVdir_v1.configure(state='normal')
            self.entry_FSVdir_v2.configure(state='normal')
            
    def check_FJC(self):
        '''Checks if fringe jump correction is selected'''
        
        if self.FJC.get() == 0:
            self.entry_m.configure(state='disabled')
            self.entry_n.configure(state='disabled')
        elif self.FJC.get() == 1:
            self.entry_m.configure(state='normal')
            self.entry_n.configure(state='normal')
        
'''Main --------------------------------------------------------------------'''        
def main():            
    root = tk.Tk()
    GUI = VisarGUI(root)
    root.mainloop()
    
if __name__ == "__main__": 
    main()