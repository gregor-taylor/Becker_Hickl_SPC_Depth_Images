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

class LIDARDataPx(Tk):
    def __init__(self, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)
        Tk.wm_title(self, "Data Px")
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
        load_file_butt=ttk.Button(self, text='Choose sdt file', command=lambda:self.load_file(controller))
        load_file_butt.grid(row=1, column=1)

        ttk.Label(self, text="Gating start/end: (leave end at 0 for no gating)").grid(row=2,column=1)
        gate_start = ttk.Entry(self, textvariable=controller.start_gate)
        gate_start.grid(row=3,column=1)
        gate_end = ttk.Entry(self, textvariable=controller.end_gate)
        gate_end.grid(row=4,column=1)

        IRF_check=ttk.Checkbutton(self, text="Automatically determine IRF?", variable=controller.IRF_var)
        IRF_check.grid(row=5,column=1)

        SBR_check=ttk.Checkbutton(self, text="Subtract back reflections?", variable=controller.SBR_var)
        SBR_check.grid(row=6, column=1)

        plot_butt=ttk.Button(self, text='Plot it', command=lambda:self.plot_it(controller))
        plot_butt.grid(row=7, column=1)

        fit_butt=ttk.Checkbutton(self, text='Fit instead of X-corr?', variable=controller.fit_var)
        fit_butt.grid(row=8, column=1)

    def load_file(self, controller):
        controller.filename = askopenfilename(initialdir="C:\\", title="Choose an sdt file")
        
    def plot_it(self, controller):
        if controller.filename=='':
            messagebox.showerror('Error', 'No sdt file loaded!')
        else:
            #Takes main file with all pixels#
            sdt_file=SdtFile(controller.filename)
            image_size_x=sdt_file.data[0].shape[0]
            image_size_y=sdt_file.data[0].shape[1]
            #Pulls the TAC paramters from the sdt file. For some reason the 'times' array from the sdt file is not the correct times - poss software not updated for new card.
            adc_re=sdt_file.measure_info[0].adc_re
            tac_r=sdt_file.measure_info[0].tac_r
            tac_g=sdt_file.measure_info[0].tac_g
            dt=tac_r/tac_g/adc_re
            times=range(0,int(sdt_file.measure_info[0].adc_re))*dt

            #sets end point to end if gating not required
            if controller.end_gate.get()=='0':
                end_point=len(sdt_file.data[0][0][0])
                start_point=0
            else:
                start_point=round(int((int(controller.start_gate.get())*1e-12)/dt)) #converts back from ps to bins
                end_point=round(int((int(controller.end_gate.get())*1e-12)/dt))

            #If removing the BR/Noise this takes a seperate sdt file#
            if controller.SBR_var.get() == 1:
                if controller.end_gate.get() != '0':
                    messagebox.showerror('Error', 'Cannot use noise removal with gating')
                    controller.SBR_var.set('0')
                else:
                    controller.noise_file = askopenfilename(title="Choose a noise file")
                    noise_file_sdt = SdtFile(controller.noise_file)
                    noise_data=noise_file_sdt.data[0][0][start_point:end_point]

            #Processes the max counts - also finds brightest pixel for IRF if desired#
            max_arr = np.zeros((image_size_y, image_size_x))
            max_count_all=0
            pixel=(0,0)
            for i in range(image_size_y):
                for j in range(image_size_x):
                    if controller.SBR_var.get() == 1:
                        sdt_file.data[0][i][j]=sdt_file.data[0][i][j]-noise_data
                    max_count=np.amax(sdt_file.data[0][i][j][start_point:end_point])
                    max_arr[i][j]=max_count
                    if max_count > max_count_all:
                        max_count_all=max_count
                        pixel=(i,j)
            #plots the brightest and fits it to see where we're at
            centre, scale = fit_exp_mod_gauss(times[start_point:end_point], sdt_file.data[0][pixel[0]][pixel[1]][start_point:end_point], dt, plotting=True)

            #Takes brightest pixel as IRF or takes external sdt file#
            if controller.IRF_var.get() == 1:
                IRF_pix = sdt_file.data[0][pixel[0]][pixel[1]][start_point:end_point]
            elif controller.fit_var.get() == 1:
                pass
            else:
                controller.IRF_file = askopenfilename(title="Choose IRF sdt file")
                IRF_file_sdt = SdtFile(controller.IRF_file)
                max_irf = np.argmax(IRF_file_sdt[0][0])
                start_irf = max_irf-100
                stop_irf = max_irf+100
                IRF_pix = IRF_file_sdt.data[0][0][start_point:end_point] 
                plt.plot(IRF_pix)
                plt.show()

            #prepare an empty 2d arr
            img_arr = np.zeros((image_size_y, image_size_x))

            #either fit or X-corr
            if controller.fit_var.get() == 1:
                for i in range(image_size_y):
                    for j in range(image_size_x):
                        try:
                            centre, scale = fit_exp_mod_gauss(times[start_point:end_point], sdt_file.data[0][i][j][start_point:end_point], dt)
                            img_arr[i,j]=centre
                        except TypeError:
                            img_arr[i,j]=float('nan')
                        except RuntimeError:
                            img_arr[i,j]=float('nan')                
            else:
                #cross correlates the data#
                for i in range(image_size_y):
                    for j in range(image_size_x):
                        corr = cross_correlate(sdt_file.data[0][i][j][start_point:end_point] , IRF_pix)
                        max_val = np.argmax(corr)
                        img_arr[i,j]=max_val       
            
            #Scale and convert to mm
            img_arr= -img_arr - np.nanmean(-img_arr)
            img_arr=img_arr*1e-6*3e8*0.5*0.5

            ####PLOTS####
            #3d#
            fig1=plt.figure()
            f1_ax=fig1.gca(projection='3d')
            x_data=np.arange(image_size_x)
            y_data=np.arange(image_size_y)
            x_data, y_data=np.meshgrid(x_data, y_data)
            surf=f1_ax.plot_surface(x_data, y_data, img_arr, cmap=cm.jet)
            fig1.colorbar(surf, shrink=0.5, aspect=5)
            fig1.suptitle('3D')
            
            #2d#
            fig2=plt.figure()
            flat_plot=plt.imshow(img_arr, cmap=cm.jet)
            fig2.colorbar(flat_plot, shrink=0.5, aspect=5)
            fig2.suptitle('2D')
            
            #counts#
            fig3=plt.figure()            
            cnt_map=plt.imshow(max_arr, cmap=cm.jet)
            fig3.colorbar(cnt_map, shrink=0.5, aspect=5)
            fig3.suptitle('Counts Map')
            
            #plots one histogram from a bright pixel to check for gating errors etc and check fit
            #fig4=plt.figure()
            #plt.plot(times[start_point:end_point], sdt_file.data[0][pixel[0]][pixel[1]][start_point:end_point],'o')
            #fig4.suptitle('Brightest histogram')
            
            plt.show()


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
    app.geometry("250x200")
    app.mainloop()

if __name__ == '__main__':
    main()

