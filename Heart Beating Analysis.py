# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 16:46:04 2021

@author: natha
"""
from scipy.fft import fft, ifft
from scipy.signal import find_peaks
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

timepoints = 4 #input('How many timepoints: ')
tps = []
treatment = 'ctrl' #string.empty
con = "collagenase" #input('Enter Conditon: ')

# if(input('Treatment or Control (t/c):')== 't'):
#      treatment = 'trtd'
# else:
#      treatment = 'ctrl'



hrt = input('Enter Sample #: ')

if(timepoints == 4):
    tps = [0, 15, 30, 45]

elif(timepoints == 3):
    tps = [0, 57600, 144000]

else:
    tps = [0, 45, 90, 180, 270]

#Set working directory
#if gell run first
if 'n' == 'y':
    os.chdir( r'C:\Users\natha\Desktop\CEMB Summer Program 2021\Discher Lab\Ex-vivo chick heart experiments\Experiments\Beating data_20210728\Data analysis\gels\{c}\{t}{h}'.format(c = con, h = hrt, t = treatment))
else:
    os.chdir( r'C:\Users\natha\Desktop\CEMB Summer Program 2021\Discher Lab\Ex-vivo chick heart experiments\Experiments\Beating data_20210728\Data analysis\{c}\{t}{h}'.format(c = con, h = hrt, t = treatment))



#%%
data_tp0 = pd.read_excel('{c}_{t}{h}_t{tp}_analysis.xlsx'.format(c = con, h = hrt, t = treatment, tp = tps[0]), usecols = "A:B")
Time_tp0 = data_tp0['Time (s)']

signal_tp0 = (data_tp0['delta AR/ARref'])
# signal = sin(Time);   %%  Raw signal with time period of ~6.5 sec i.e. Freq = 1/6.5 = 0.15 Hz
L = len(Time_tp0)
Fs= L/max(Time_tp0)     #% Sampling frequency                    
T = 1/Fs #             % Sampling period  


Y_tp0 = fft(signal_tp0.to_numpy())         #Compute DFT of Signal and get real and imaginary components
P2_tp0 = abs(Y_tp0/L)
                                           #Whole spectrum
P1_tp0 = 2*P2_tp0[:int(L/2)]
P1_tp0[0] = P1_tp0[0]/2
                                           #      % One-sided spectrum

#P1_tp0(2:end-1) = 2*P1_tp0(2:end-1);              % Even-valued signal length L.

#Filter heights less than this when doing peak to peak analysis
filterfraction = 0.66;

# Obtain peaks
minimum_peakheight_tp0 = filterfraction*max(signal_tp0)
peaks_tp0 = find_peaks(signal_tp0, height = minimum_peakheight_tp0)
time_peaks_tp0 = Time_tp0[peaks_tp0[0]]

# %%To get valleys, you invert the signal
invertedsignal_tp0 = max(signal_tp0) - signal_tp0
minimum_valleyheight_tp0 = filterfraction*max(invertedsignal_tp0);
valleys_tp0 = find_peaks(invertedsignal_tp0, height = minimum_valleyheight_tp0)
time_valleys_tp0 = Time_tp0[valleys_tp0[0]]
Values_valleys_tp0 = max(signal_tp0)-valleys_tp0[1]['peak_heights']  #%%%  True valley values
amp_tp0 = peaks_tp0[1]['peak_heights'].mean()-Values_valleys_tp0.mean()

#%%
data_tp1 = pd.read_excel('{c}_{t}{h}_t{tp}_analysis.xlsx'.format(c = con, h = hrt, t = treatment, tp = tps[1]), usecols = "A:B")
Time_tp1 = data_tp1['Time (s)']

signal_tp1 = (data_tp1['delta AR/ARref'])
# signal = sin(Time);   %%  Raw signal with time period of ~6.5 sec i.e. Freq = 1/6.5 = 0.15 Hz
L = len(Time_tp1)
Fs= L/max(Time_tp1)     #% Sampling frequency                    
T = 1/Fs #             % Sampling period  


Y_tp1 = fft(signal_tp1.to_numpy())         #Compute DFT of Signal and get real and imaginary components
P2_tp1 = abs(Y_tp1/L)
                                           #Whole spectrum
P1_tp1 = 2*P2_tp1[:int(L/2)]
P1_tp1[0] = P1_tp1[0]/2
                                           #      % One-sided spectrum

#P1_tp1(2:end-1) = 2*P1_tp1(2:end-1);              % Even-valued signal length L.

#Filter heights less than this when doing peak to peak analysis
filterfraction = 0.66;

# Obtain peaks
minimum_peakheight_tp1 = filterfraction*max(signal_tp1)
peaks_tp1 = find_peaks(signal_tp1, height = minimum_peakheight_tp1)
time_peaks_tp1 = Time_tp1[peaks_tp1[0]]

# %%To get valleys, you invert the signal
invertedsignal_tp1 = max(signal_tp1) - signal_tp1
minimum_valleyheight_tp1 = filterfraction*max(invertedsignal_tp1);
valleys_tp1 = find_peaks(invertedsignal_tp1, height = minimum_valleyheight_tp1)
time_valleys_tp1 = Time_tp1[valleys_tp1[0]]
Values_valleys_tp1 = max(signal_tp1)-valleys_tp1[1]['peak_heights']  #%%%  True valley values
amp_tp1 = peaks_tp1[1]['peak_heights'].mean()-Values_valleys_tp1.mean()





#%%
data_tp2 = pd.read_excel('{c}_{t}{h}_t{tp}_analysis.xlsx'.format(c = con, h = hrt, t = treatment, tp = tps[2]), usecols = "A:B")
Time_tp2 = data_tp2['Time (s)']

signal_tp2 = (data_tp2['delta AR/ARref'])
# signal = sin(Time);   %%  Raw signal with time period of ~6.5 sec i.e. Freq = 1/6.5 = 0.15 Hz
L = len(Time_tp2)
Fs= L/max(Time_tp2)     #% Sampling frequency                    
T = 1/Fs #             % Sampling period  


Y_tp2 = fft(signal_tp2.to_numpy())         #Compute DFT of Signal and get real and imaginary components
P2_tp2 = abs(Y_tp2/L)
                                           #Whole spectrum
P1_tp2 = 2*P2_tp2[:int(L/2)]
P1_tp2[0] = P1_tp2[0]/2
                                           #      % One-sided spectrum

#P1_tp2(2:end-1) = 2*P1_tp2(2:end-1);              % Even-valued signal length L.

#Filter heights less than this when doing peak to peak analysis
filterfraction = 0.66;

# Obtain peaks
minimum_peakheight_tp2 = filterfraction*max(signal_tp2)
peaks_tp2 = find_peaks(signal_tp2, height = minimum_peakheight_tp2)
time_peaks_tp2 = Time_tp2[peaks_tp2[0]]

# %%To get valleys, you invert the signal
invertedsignal_tp2 = max(signal_tp2) - signal_tp2
minimum_valleyheight_tp2 = filterfraction*max(invertedsignal_tp2);
valleys_tp2 = find_peaks(invertedsignal_tp2, height = minimum_valleyheight_tp2)
time_valleys_tp2 = Time_tp2[valleys_tp2[0]]
Values_valleys_tp2 = max(signal_tp2)-valleys_tp2[1]['peak_heights']  #%%%  True valley values
amp_tp2 = peaks_tp2[1]['peak_heights'].mean()-Values_valleys_tp2.mean()





#%%
if(timepoints == 4):
    data_tp3 = pd.read_excel('{c}_{t}{h}_t{tp}_analysis.xlsx'.format(c = con, h = hrt, t = treatment, tp = tps[3]), usecols = "A:B")
    Time_tp3 = data_tp3['Time (s)']
    
    signal_tp3 = (data_tp3['delta AR/ARref'])
    # signal = sin(Time);   %%  Raw signal with time period of ~6.5 sec i.e. Freq = 1/6.5 = 0.15 Hz
    L = len(Time_tp3)
    Fs= L/max(Time_tp3)     #% Sampling frequency                    
    T = 1/Fs #             % Sampling period  
    
    
    Y_tp3 = fft(signal_tp3.to_numpy())         #Compute DFT of Signal and get real and imaginary components
    P2_tp3 = abs(Y_tp3/L)
                                               #Whole spectrum
    P1_tp3 = 2*P2_tp3[:int(L/2)]
    P1_tp3[0] = P1_tp3[0]/2
                                               #      % One-sided spectrum
    
    #P1_tp3(2:end-1) = 2*P1_tp3(2:end-1);              % Even-valued signal length L.
    
    #Filter heights less than this when doing peak to peak analysis
    filterfraction = 0.66;
    
    # Obtain peaks
    minimum_peakheight_tp3 = filterfraction*max(signal_tp3)
    peaks_tp3 = find_peaks(signal_tp3, height = minimum_peakheight_tp3)
    time_peaks_tp3 = Time_tp3[peaks_tp3[0]]
    
    # %%To get valleys, you invert the signal
    invertedsignal_tp3 = max(signal_tp3) - signal_tp3
    minimum_valleyheight_tp3 = filterfraction*max(invertedsignal_tp3);
    valleys_tp3 = find_peaks(invertedsignal_tp3, height = minimum_valleyheight_tp3)
    time_valleys_tp3 = Time_tp3[valleys_tp3[0]]
    Values_valleys_tp3 = max(signal_tp3)-valleys_tp3[1]['peak_heights']  #%%%  True valley values
    amp_tp3 = peaks_tp3[1]['peak_heights'].mean()-Values_valleys_tp3.mean()

#%%
if(timepoints == 5):
    data_tp4 = pd.read_excel('{c}_{t}{h}_t{tp}_analysis.xlsx'.format(c = con, h = hrt, t = treatment, tp = tps[4]), usecols = "A:B")
    Time_tp4 = data_tp4['Time (s)']
    
    signal_tp4 = (data_tp4['delta AR/ARref'])
    # signal = sin(Time);   %%  Raw signal with time period of ~6.5 sec i.e. Freq = 1/6.5 = 0.15 Hz
    L = len(Time_tp4)
    Fs= L/max(Time_tp4)     #% Sampling frequency                    
    T = 1/Fs #             % Sampling period  
    
    
    Y_tp4 = fft(signal_tp4.to_numpy())         #Compute DFT of Signal and get real and imaginary components
    P2_tp4 = abs(Y_tp4/L)
                                               #Whole spectrum
    P1_tp4 = 2*P2_tp4[:int(L/2)]
    P1_tp4[0] = P1_tp4[0]/2
                                               #      % One-sided spectrum
    
    #P1_tp4(2:end-1) = 2*P1_tp4(2:end-1);              % Even-valued signal length L.
    
    #Filter heights less than this when doing peak to peak analysis
    filterfraction = 0.66;
    
    # Obtain peaks
    minimum_peakheight_tp4 = filterfraction*max(signal_tp4)
    peaks_tp4 = find_peaks(signal_tp4, height = minimum_peakheight_tp4)
    time_peaks_tp4 = Time_tp4[peaks_tp4[0]]
    
    # %%To get valleys, you invert the signal
    invertedsignal_tp4 = max(signal_tp4) - signal_tp4
    minimum_valleyheight_tp4 = filterfraction*max(invertedsignal_tp4);
    valleys_tp4 = find_peaks(invertedsignal_tp4, height = minimum_valleyheight_tp4)
    time_valleys_tp4 = Time_tp4[valleys_tp4[0]]
    Values_valleys_tp4 = max(signal_tp4)-valleys_tp4[1]['peak_heights']  #%%%  True valley values
    amp_tp4 = peaks_tp4[1]['peak_heights'].mean()-Values_valleys_tp4.mean()

# %%% Write to excel
f = Fs*np.array(range(0,int(L/2)))/L

if len(tps) == 3:
   freq = [f[np.argmax(P1_tp0)], f[np.argmax(P1_tp1)], f[np.argmax(P1_tp2)]]
   strains = [amp_tp0, amp_tp1, amp_tp2]
   frequency_distribution = pd.DataFrame([f, P1_tp0, P1_tp1, P1_tp2], index = ['frequncies', 'magnitude_tp0','magnitude_tp1','magnitude_tp2'])

elif len(tps) == 4:
    freq = [f[np.argmax(P1_tp0)], f[np.argmax(P1_tp1)], f[np.argmax(P1_tp2)], f[np.argmax(P1_tp3)]]
    strains = [amp_tp0, amp_tp1, amp_tp2,amp_tp3]
    frequency_distribution = pd.DataFrame([f, P1_tp0, P1_tp1, P1_tp2, P1_tp3], index = ['frequncies', 'magnitude_tp0','magnitude_tp1','magnitude_tp2','magnitude_tp3'])
    
elif len(tps) == 5:
    freq = [f[np.argmax(P1_tp0)], f[np.argmax(P1_tp1)], f[np.argmax(P1_tp2)], f[np.argmax(P1_tp3)],f[np.argmax(P1_tp4)] ]
    strains = [amp_tp0, amp_tp1, amp_tp2,amp_tp3, amp_tp4]
    frequency_distribution = pd.DataFrame([f, P1_tp0, P1_tp1, P1_tp2, P1_tp3, P1_tp4], index = ['frequncies', 'magnitude_tp0','magnitude_tp1','magnitude_tp2','magnitude_tp3','magnitude_tp4'])
    
analysis = {'Time': tps , 'Frequency': freq, 'Beating Strain': strains}
Analysis_Summary = pd.DataFrame(analysis)


 

#make Excel writer with user input for file name and workbook object
writer = pd.ExcelWriter('{c}_{t}{h}_analysis.xlsx'.format(h = hrt, c = con, t = treatment))
workbook = writer.book

##########################################################################
Analysis_Summary.to_excel(excel_writer = writer, sheet_name = 'Summary' , index = True)
frequency_distribution.to_excel(excel_writer = writer, sheet_name = 'Freq_Distribution' , index = True)

summary =writer.sheets['Summary']
freqDisSheet = writer.sheets['Freq_Distribution']

#create chart object for strain over time
strainOt = workbook.add_chart({'type': 'scatter','subtype': 'straight_with_markers'})
strainOt.add_series({
    'categories': ['Summary', 1, 1, len(tps)+1, 1],
    'values':     ['Summary', 1, 2, len(tps)+1, 2],
})

#create chart object for frequency over time
freqOt = workbook.add_chart({'type': 'scatter','subtype': 'straight_with_markers'})
freqOt.add_series({
    'categories': ['Summary', 1, 1, len(tps)+1, 1],
    'values':     ['Summary', 1, 3, len(tps)+1, 3],
})


#create chart object for Frequency distribution
freqDis = workbook.add_chart({'type': 'scatter','subtype': 'straight'})

#Add frequency distribution for tp0
freqDis.add_series({
    'name':       'tp0',
    'categories': ['Freq_Distribution', 1, 1, 1, len(f)+1],
    'values':     ['Freq_Distribution', 2, 1, 2, len(P1_tp0)+1],
})

#Add frequency distribution for tp1
freqDis.add_series({
    'name':       'tp1',
    'categories': ['Freq_Distribution', 1, 1, 1, len(f)+1],
    'values':     ['Freq_Distribution', 3, 1, 3, len(P1_tp1)+1],
})

#Add frequency distribution for tp2
freqDis.add_series({
    'name':       'tp2',
    'categories': ['Freq_Distribution', 1, 1, 1, len(f)+1],
    'values':     ['Freq_Distribution', 4, 1, 4, len(P1_tp2)+1],
})

if len(tps) ==4:
#Add frequency distribution for tp3
    freqDis.add_series({
        'name':       'tp3',
        'categories': ['Freq_Distribution', 1, 1, 1, len(f)+1],
        'values':     ['Freq_Distribution', 5, 1, 5, len(P1_tp3)+1],
    })


# Set name on axis of colChartBkwd_Fwd object and insert to 920vs860nm sheet
freqDis.set_title({'name': 'Frequency Distribution of {c} {t}{h}'.format(c = con, h = hrt, t =treatment)})
freqDis.set_x_axis({'name': 'Frequency (Hz)', 'min': 0, 'max': 5})
freqDis.set_y_axis({'name': 'Magnitude','min': 0, 'max': 2,
                  'major_gridlines': {'visible': True},
})

# Set name on axis of colChartBkwd_Fwd object and insert to 920vs860nm sheet
strainOt.set_title({'name': 'Strain over time {c} {t}{h}'.format(c = con, h = hrt, t =treatment)})
strainOt.set_x_axis({'name': 'Time (min)', 'min': 0, 'max': tps[-1]+15})
strainOt.set_y_axis({'name': 'AR/AR_ref','min': 0, 'max': 5,
                  'major_gridlines': {'visible': True},
})
strainOt.set_legend({'none': True})


# Set name on axis of colChartBkwd_Fwd object and insert to 920vs860nm sheet
freqOt.set_title({'name': 'Frequency over time {c} {t}{h}'.format(c = con, h = hrt, t =treatment)})
freqOt.set_x_axis({'name': 'Time (min)', 'min': 0, 'max': tps[-1]+15})
freqOt.set_y_axis({'name': 'Frequency (Hz)','min': 0, 'max': 6,
                  'major_gridlines': {'visible': True},
})
freqOt.set_legend({'none': True})

#Insert chartsheet
freqDisSheet.insert_chart('G8', freqDis)
summary.insert_chart('F1',strainOt)
summary.insert_chart('F17',freqOt)

#close excel file       
writer.close()
