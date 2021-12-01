#############################
#G Taylor July 2019
#Takes sdt file output from the Beckr&Hickl card running in FiFo Image/Scan Sync In mode
#Cross correlates an IRF with each pixel and produces depth profiles
#
#
#See requirements.txt for dependencies
#
#############################
import numpy as np
import tkinter as tk
from tkinter import Tk, Frame, ttk, messagebox
from tkinter.filedialog import askopenfilename, askdirectory
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from sdtfile import SdtFile
from scipy.signal import correlate, savgol_filter
from scipy.special import erfc
from scipy.optimize import curve_fit
from math import sqrt
import csv
import plotly.graph_objects as go
import time

class LIDARDataPx(Tk):
    def __init__(self, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)
        Tk.wm_title(self, "LIDAR Data Processing")
        #Tk.iconbitmap(self, default="UD_Smiley.ico")
        #Data holders
        self.filename=''
        self.IRF_file = ''
        self.noise_file=''
        self.IRF_var = tk.IntVar()
        self.SBR_var=tk.IntVar()
        self.fit_var=tk.IntVar()
        self.start_gate = tk.StringVar()
        self.start_gate.set('0')
        self.end_gate = tk.StringVar()
        self.end_gate.set('0')
        self.cut_off_max = tk.StringVar()
        self.cut_off_max.set('0')
        self.scene_size = tk.StringVar()
        self.scene_size.set('0')
        self.aspect_togg=tk.IntVar()
        self.remove_empty_pix=tk.IntVar()
        self.bin_toggle=tk.IntVar()
        self.bin_factor=tk.StringVar()
        self.bin_factor.set('0')

        container=ttk.Frame(self)
        container.pack(side='top', fill='both', expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames={}
        self.frames[MainPage]=MainPage(container, self)
        self.frames[MainPage].grid(row=0, column=0, sticky='nsew')

        self.show_frame(MainPage)

    def show_frame(self, cont):
        #Raises the chosen page to the top.
        frame = self.frames[cont]
        frame.tkraise()

class MainPage(ttk.Frame):
    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        load_file_butt=ttk.Button(self, text='Choose sdt/npy file', command=lambda:self.load_file(controller))
        load_file_butt.grid(row=1, column=1)

        ttk.Label(self, text="Gating start/end: (leave end at 0 for no gating)").grid(row=2,column=1)
        gate_start = ttk.Entry(self, textvariable=controller.start_gate)
        gate_start.grid(row=3,column=1)
        gate_end = ttk.Entry(self, textvariable=controller.end_gate)
        gate_end.grid(row=4,column=1)

        ttk.Label(self, text="Largest value cutoff? (um, leave at 0 for no cut off):").grid(row=3, column=2)
        cut_off_max = ttk.Entry(self, textvariable=controller.cut_off_max)
        cut_off_max.grid(row=4, column=2)

        IRF_check=ttk.Checkbutton(self, text="Automatically determine IRF?", variable=controller.IRF_var)
        IRF_check.grid(row=5,column=1)

        SBR_check=ttk.Checkbutton(self, text="Subtract back reflections?", variable=controller.SBR_var)
        SBR_check.grid(row=6, column=1)

        plot_butt=ttk.Button(self, text='Plot it', command=lambda:self.plot_it(controller))
        plot_butt.grid(row=8, column=1, columnspan=2)

        fit_butt=ttk.Checkbutton(self, text='Fit instead of X-corr?', variable=controller.fit_var)
        fit_butt.grid(row=7, column=1)

        aspect_correct_toggle=ttk.Checkbutton(self, text='Correct aspect in graph? Input scene size (mm):', variable=controller.aspect_togg)
        aspect_correct_toggle.grid(row=1, column=2)

        scene_size = ttk.Entry(self, textvariable=controller.scene_size)
        scene_size.grid(row=2, column=2)

        remove_empty_pix_toggle = ttk.Checkbutton(self, text='Remove empty pixels from plots?', variable=controller.remove_empty_pix)
        remove_empty_pix_toggle.grid(row=5, column=2)

        binning_toggle=ttk.Checkbutton(self, text='Bin the data? Input bin factor:', variable=controller.bin_toggle)
        binning_toggle.grid(row=6, column=2)

        bin_factor=ttk.Entry(self, textvariable=controller.bin_factor)
        bin_factor.grid(row=7, column=2)


    def load_file(self, controller):
        controller.filename = askopenfilename(initialdir="C:\\", title="Choose an sdt/npy file")
        
    def plot_it(self, controller):
        if controller.filename=='':
            messagebox.showerror('Error', 'No sdt file loaded!')
        else:
            #Checks if .sdt or .npy
            if controller.filename[-3:] == 'sdt':
                #Takes main file with all pixels#
                sdt_file=SdtFile(controller.filename)
                #Pulls the TAC paramters from the sdt file. For some reason the 'times' array from the sdt file is not the correct times - poss software not updated for new card.
                adc_re=sdt_file.measure_info[0].adc_re
                tac_r=sdt_file.measure_info[0].tac_r
                tac_g=sdt_file.measure_info[0].tac_g
                image_data=sdt_file.data[0]
                image_size_x=image_data.shape[0]
                image_size_y=image_data.shape[1]
            elif controller.filename[-3:] == 'npy':
                temp_data=np.load(controller.filename, allow_pickle=True)
                image_size_x=temp_data.shape[0]
                image_size_y=temp_data.shape[1]
                image_data=np.ndarray((image_size_x, image_size_y), dtype='object')
                for i in range(image_size_y):
                    for j in range(image_size_x):
                        image_data[i][j]=temp_data[i][j][0]
                #Need to update to pull TAC parameters from .set file, or bundle them with the numpy data
                adc_re=4096
                tac_r=2.5016787e-8
                tac_g=15
            else:
                messagebox.showerror('Error', 'Invalid filetype!')
            #Gets image size and cimputes dt/times arrays
            
            dt=tac_r/tac_g/adc_re
            times=np.arange(0,int(adc_re))*dt
            #Binning
            if controller.bin_toggle.get() == 1:
                processed_data = np.ndarray((image_size_x, image_size_y), dtype='object') #new array for data as new shape
                bin_factor=int(controller.bin_factor.get())
                image_data=image_data.reshape(*image_data.shape[:2], -1, bin_factor).sum(axis=-1)
                dt=dt*bin_factor
                times=np.arange(0,int(adc_re/bin_factor))*dt

            #sets end point to end if gating not required
            if controller.end_gate.get()=='0':
                end_point=len(image_data[0][0])
                start_point=0
            else:
                start_point=round(int((int(controller.start_gate.get())*1e-12)/dt)) #converts back from ps to bins
                end_point=round(int((int(controller.end_gate.get())*1e-12)/dt))
                image_data=image_data[:,:,start_point:end_point]
                times=times[start_point:end_point]

            #If removing the BR/Noise this takes a seperate sdt file#
            if controller.SBR_var.get() == 1:
                if controller.end_gate.get() != '0':
                    messagebox.showerror('Error', 'Cannot use noise removal with gating')
                    controller.SBR_var.set('0')
                else:
                    controller.noise_file = askopenfilename(title="Choose a noise file")
                    noise_file_sdt = SdtFile(controller.noise_file)
                    noise_data=noise_file_sdt.data[0][0][start_point:end_point]
                    #add binning here if we use this

            #Processes the max counts - also finds brightest pixel for IRF if desired#
            #Also notes any pixels with no real counts in it and notes the ID's for post px-ing
            max_arr = np.zeros((image_size_y, image_size_x))
            max_count_all=0
            pixel=(0,0)
            for i in range(image_size_y):
                for j in range(image_size_x):
                    if controller.SBR_var.get() == 1:
                        image_data[i][j]=image_data[i][j]-noise_data
                    max_count=np.amax(image_data[i][j])
                    max_arr[i][j]=max_count
                    if max_count > max_count_all:
                        max_count_all=max_count
                        pixel=(i,j)
            empty_pixels=np.transpose(np.nonzero(max_arr < 2))
            #plots the brightest and fits it to see where we're at
            centre, scale = fit_exp_mod_gauss(times, image_data[pixel[0]][pixel[1]], dt, plotting=True)
            #Checks if you're happy with gating?
            MsgBox = tk.messagebox.askquestion ('Proceed?','Are you happy with the gating?',icon = 'warning')
            if MsgBox == 'yes':
                #Takes brightest pixel as IRF or takes external sdt file#
                if controller.IRF_var.get() == 1:
                    IRF_pix = image_data[pixel[0]][pixel[1]]
                elif controller.fit_var.get() == 1:
                    pass
                else:
                    controller.IRF_file = askopenfilename(title="Choose IRF sdt file")
                    IRF_file_sdt = SdtFile(controller.IRF_file)
                    max_irf = np.argmax(IRF_file_sdt.data[0][0])
                    half_gate = ((int(controller.end_gate.get())-int(controller.start_gate.get()))/2)
                    start_irf = int(max_irf-half_gate)
                    stop_irf = int(max_irf+half_gate)
                    if controller.bin_toggle.get() == 1:
                        IRF_pix = IRF_file_sdt.data[0][0]
                        IRF_pix=IRF_pix.reshape(int(adc_re/bin_factor),-1).sum(axis=1)
                        IRF_pix=IRF_pix[start_irf:stop_irf]
                    else:
                        IRF_pix = IRF_file_sdt.data[0][0][start_irf:stop_irf] 
                    #Add binning here too
                    plt.plot(IRF_pix)
                    plt.show()    

                #prepare an empty 2d arr
                img_arr = np.zeros((image_size_y, image_size_x))    

                #either fit or X-corr
                if controller.fit_var.get() == 1:
                    for i in range(image_size_y):
                        for j in range(image_size_x):
                            try:
                                centre, scale = fit_exp_mod_gauss(times, image_data[i][j], dt)
                                img_arr[i,j]=centre
                            except TypeError:
                                img_arr[i,j]=float('nan')
                            except RuntimeError:
                                img_arr[i,j]=float('nan')                
                else:
                    t0=time.time()
                    #cross correlates the data#
                    for i in range(image_size_y):
                        for j in range(image_size_x):
                            corr = cross_correlate(image_data[i][j] , IRF_pix)
                            max_val = np.argmax(corr)
                            img_arr[i,j]=max_val       
                
                #Scale and convert to mm
                img_arr= -img_arr - np.nanmean(-img_arr)
                img_arr=img_arr*dt*1e6*3e8*0.5
                #Bit of data manipulation for any pixels that are huge
                if controller.cut_off_max.get() != '0':
                    super_threshold_indices = img_arr > float(controller.cut_off_max.get())
                    img_arr[super_threshold_indices] = np.mean(img_arr) 

                #If pixels have no counts (or just noise) then leave them out.
                if controller.remove_empty_pix.get() == 1:
                    for i in empty_pixels:
                        img_arr[i[0], i[1]] = np.nan

                ####PLOTS####
                #3d#
               # fig1=plt.figure()
               # f1_ax=fig1.gca(projection='3d')
               # x_data=np.arange(image_size_x)
               # y_data=np.arange(image_size_y)
               # x_data, y_data=np.meshgrid(x_data, y_data)
               # surf=f1_ax.plot_surface(x_data, y_data, img_arr, cmap=cm.jet)
               # fig1.colorbar(surf, shrink=0.5, aspect=5)
               # fig1.suptitle('3D')
                
                #2d#
               # fig2=plt.figure()
               # flat_plot=plt.imshow(img_arr, cmap=cm.jet)
               # fig2.colorbar(flat_plot, shrink=0.5, aspect=5)
               # fig2.suptitle('2D')
                
                #counts#
                tfin=time.time()-t0
                print('Time to process: {} seconds'.format(tfin))
                fig3=plt.figure()            
                cnt_map=plt.imshow(max_arr, cmap=cm.jet, origin='lower')
                fig3.colorbar(cnt_map, shrink=0.5, aspect=5)
                print(np.mean(max_arr))
                #fig3.suptitle('Counts Map')    

                fig5=go.Figure(data=[go.Heatmap(z=img_arr, colorscale='Jet', colorbar=dict(thickness=80,
                               ticklen=3, tickcolor='black',
                               tickfont=dict(size=36, color='black')))])
                fig5.update_layout(width = 1500, height = 1500, font=dict(family="Arial",size=36,color="black"))
                fig5.show()    

                if controller.aspect_togg.get() == 1:
                    if controller.scene_size.get() == '':
                        messagebox.showerror('Error', 'You have not entered a scene size!')
                        aspect_dict=dict(x=1,y=1,z=1)
                    else:
                        scene_size=int(controller.scene_size.get()) #mm
                        aspect_dict=dict(x=1, y=1, z=np.nanmax(img_arr)/(scene_size*1e3))
                else:
                    aspect_dict=dict(x=1,y=1,z=1)
            
                fig6=go.Figure(data=[go.Surface(z=img_arr, colorscale='Jet', colorbar=dict(thickness=40,
                               ticklen=3, tickcolor='black',
                               tickfont=dict(size=9, color='black')))])
                fig6.update_layout(font=dict(family="Arial",size=9,color="black"), scene=dict(aspectmode='manual',
                               aspectratio=aspect_dict))
                fig6.show()
                
                #plots one histogram from a bright pixel to check for gating errors etc and check fit
                #fig4=plt.figure()
                #plt.plot(times[start_point:end_point], sdt_file.data[0][pixel[0]][pixel[1]][start_point:end_point],'o')
                #fig4.suptitle('Brightest histogram')
                
                plt.show()
            else:
                pass


#X-corr function
def cross_correlate(data_1, data_2):
    #corr = np.fft.ifft(np.fft.fft(data_1)*np.conj(np.fft.fft(data_2)))
    corr = correlate(data_1, data_2, mode='same')
    #plt.plot(corr)
    #plt.show()
    return corr

def fit_exp_mod_gauss(times, counts, dt, plotting=False):
    #convert times to ps
    times=times/1e-12
    #use savitzky-golay poly to smooth data for guess
    smoothed_counts=savgol_filter(counts, 5, 4)
    #normalise by maxcounts then shift max counts bin to 0 - needed?
    max_ind = np.argmax(smoothed_counts)
    max_val = smoothed_counts[max_ind]
    time_shift=times[max_ind]
    times=times-time_shift
    smoothed_counts=smoothed_counts/max_val
    if max_val == 0:
        max_val = 1
    counts=counts/max_val

    try:
        dis_width=((np.where(smoothed_counts>=0.5))[0][-1]-(np.where(smoothed_counts>=0.5)[0][0]))
    except IndexError:
        dis_width=50
    #fwhm_guess =dis_width*dt
    #print(fwhm_guess)
    #work out where to start/end fittingbased on width
    start_p=max_ind-(5*dis_width)
    end_p=max_ind+(5*dis_width)
    sliced_times=times[start_p:end_p]
    if len(sliced_times) < 6:
        sliced_times = times
        start_p = 0
        end_p = len(times)
    try:
        popt, pcov = curve_fit(exp_mod_gauss, sliced_times, counts[start_p:end_p], p0=[1,-1, 1, 0])#,  bounds=((-10, 0, -np.inf),(10, 2000, np.inf)))
    except ValueError:
        print(counts[start_p:end_p])
    fitted_func = exp_mod_gauss(sliced_times, popt[0], popt[1], popt[2], popt[3])
    centre=sliced_times[np.argmax(fitted_func)]+time_shift
    scale=np.amax(fitted_func)*max_val

    if plotting==True:
        plt.plot(times+time_shift, counts*max_val,'o', markersize=1)
        plt.plot(sliced_times+time_shift,fitted_func*max_val/max(fitted_func),'-')
        plt.title('Datpoints and EMG fit')
        plt.legend(['Data', 'Fit'])
        plt.show()

    return centre, scale

def exp_mod_gauss(x, b, m, s, l):
    y = b*(0.5*l*np.exp(0.5*l*(2*m+l*s*s-2*x))*erfc((m+l*s*s-x)/(np.sqrt(2)*s)))
    return y
    #l=Lambda, s=Sigma, m=Mu, #b=scaling

def main():
    app = LIDARDataPx()
    app.geometry("650x250")
    app.mainloop()

if __name__ == '__main__':
    main()

