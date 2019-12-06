#############################
#G Taylor
#SDT file reader
#requires sdtfile library installed - 'pip install sdtfile'
#
#Don't use gating and noise subtraction simultaneously as it wont work
#############################
import numpy as np
import tkinter as tk
from tkinter.filedialog import askopenfilename, askdirectory
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from sdtfile import SdtFile

#X-corr function
def cross_correlate(data_1, data_2):
    corr = np.fft.ifft(np.fft.fft(data_1)*np.conj(np.fft.fft(data_2)))
    return corr

####SETUP####
#If gating required. If no gating required leave end_point at 0:
start_point=0
end_point=0
#Select if not automatically getting IRF from brightest pixel
external_IRF = False
#Select to remove back reflections from data and avoid gating - NOT TESTED
remove_BR=False

####MAIN####
root = tk.Tk()
filename = askopenfilename(initialdir="Z:\\", title="Choose sdt file")
root.withdraw()

#Takes main file with all pixels#
sdt_file=SdtFile(filename)
image_size_x=sdt_file.data[0].shape[0]
image_size_y=sdt_file.data[0].shape[1]
#sets end point to end if gating not required
if end_point==0:
    end_point=len(sdt_file.data[0][0][0])

#If removing the BR/Noise this takes a seperate sdt file#
if remove_BR == True:
    noise_file = askopenfilename(title="Choose a noise file")
    noise_file_sdt = SdtFile(noise_file)
    noise_data=noise_file_sdt.data[0][0][start_point:end_point]

#Processes the max counts - also finds brightest pixel for IRF if desired#
max_arr = np.zeros((image_size_y, image_size_x))
max_count_all=0
pixel=(0,0)
for i in range(image_size_y):
    for j in range(image_size_x):
        if remove_BR==True:
            sdt_file.data[0][i][j]=sdt_file.data[0][i][j]-noise_data
        max_count=np.amax(sdt_file.data[0][i][j][start_point:end_point])
        max_arr[i][j]=max_count
        if max_count > max_count_all:
            max_count_all=max_count
            pixel=(i,j)
#Takes brightest pixel as IRF or takes external sdt file#
if external_IRF == False:
    IRF_pix = sdt_file.data[0][pixel[0]][pixel[1]][start_point:end_point]
else:
    IRF_file = askopenfilename(title="Choose IRF sdt file")
    IRF_file_sdt = SdtFile(IRF_file)
    IRF_pix = IRF_file_sdt.data[0][0][start_point:end_point] 

#cross correlates the data#
img_arr = np.zeros((image_size_y, image_size_x))
for i in range(image_size_y):
    for j in range(image_size_x):
        corr = cross_correlate(sdt_file.data[0][i][j][start_point:end_point] , IRF_pix)
        max_val = np.argmax(corr)
        img_arr[i,j]=max_val

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

#plots one histogram from a bright pixel to check for gating errors etc
fig4=plt.figure()
plt.plot(sdt_file.data[0][pixel[0]][pixel[1]][start_point:end_point] )
fig4.suptitle('Brightest histogram')

plt.show()

