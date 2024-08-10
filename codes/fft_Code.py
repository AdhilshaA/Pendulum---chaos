import numpy as np
from scipy.fft import fft
from scipy.fft import fftfreq
import matplotlib.pyplot as plt
from library_v2 import parse_trajectory_data
from scipy.signal.windows import blackman, hann

base_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\endsem" #where your trial.txt files are.
graphs_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\endsem" #where you want to save the graphs
file_name = 'trial2' # name of the datafile (no .txt)
offset = 'auto' # no. of lines to skip in the datafile (including the header)
                # if using auto, first go open the txt file see how many lines to skip and add it after '#multi:' in same first line of data
                # for e.g. '#multi: 3' means 3 lines to skip, auto will recognize this and skip lines
                # otherwise specify value of offset manually here

save_fig = True # save the figures or not (True/False)
freq_limit = 10 # in Hz for limit on graph   
window = True # apply window or not (True/False)
# which mode we are processing

# Parse the datafile into a matrix
data1 = parse_trajectory_data(f'{base_dir}\{file_name}.txt', offset = offset, updates=True)
print(f"Data parsed !")

# creating t, theta1, theta2 matrix
print(f"Creating t, theta1, theta2 matrix !")

t = data1[:,0] #adding time column
theta1 =  np.degrees(np.arctan2(data1[:,1] , (- data1[:,2]))) # theta 1 calculation
theta2 =  np.degrees(np.arctan2((data1[:,3]- data1[:,1]) , (data1[:,2]- data1[:,4]))) # theta 2 calculation
print(f"t, theta1, theta2 matrix created !")

omega1 = np.gradient(theta1, t) # omega 1 calculation
omega2 = np.gradient(theta2, t) # omega 2 calculation

# setting the initial values of omega1 and omega2 to zero
omega1[0] = 0
omega2[0] = 0

# plotting t vs theta1 and t vs theta2
print(f"Plotting t vs theta1 and t vs theta2 !")
plt.plot(t, theta1, label = r'$\theta_1$')
plt.plot(t, theta2, label = r'$\theta_2$')
plt.xlabel('t (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta_1$ and $\theta_2$ vs t -  {}'.format(file_name))
plt.legend()
if save_fig: plt.savefig(f'{graphs_dir}\{file_name}_trajectory.png')
plt.show()
print(f"t vs theta1 and t vs theta2 plotted !")

# plotting t vs omega1 and t vs omega2
print(f"Plotting t vs omega1 and t vs omega2 !")
plt.plot(t, omega1, label = r'$\omega_1$')
plt.plot(t, omega2, label = r'$\omega_2$')
plt.xlabel('t (s)')
plt.ylabel(r'$\omega$ (degrees/s)')
plt.title(r'$\omega_1$ and $\omega_2$ vs t - {}'.format(file_name))
plt.legend()
if save_fig: plt.savefig(f'{graphs_dir}\{file_name}_omega.png')
plt.show()

if window:
    N = len(theta1) # number of samples
    window = blackman(N) # window function
    
    print(f"Applying FFT on the both angles with window !")
    fft_matrix1 = fft(theta1 * window)
    fft_matrix2 = fft(theta2 * window)
    print(f"FFT applied on the both angles with window !")
    
    freq = fftfreq(t.shape[0], 1/60) # getting related frequencies for FFT
    
    # # plotting t vs theta1 and t vs theta2 windowed
    # plt.plot(theta1 * window, label = r'$\theta_1$')
    # plt.plot(theta2 * window, label = r'$\theta_2$')
    # plt.xlabel('t (s)')
    # plt.ylabel(r'$\theta$ (degrees)')
    # plt.title(r'$\theta_1$ and $\theta_2$ vs t - {} window applied'.format(file_name))
    # plt.legend()
    # plt.show()
else:
    print(f"Applying FFT on the both angles !")
    # Apply FFT on the matrix
    fft_matrix1 = fft(theta1)
    fft_matrix2 = fft(theta2)
    print(f"FFT applied on the both angles !")


    freq = fftfreq(t.shape[0], 1/60) # getting related frequencies for FFT
    

normal_freq_theta1 = freq[np.argmax(np.abs(fft_matrix1))] # normal frequency for theta1
normal_freq_theta2 = freq[np.argmax(np.abs(fft_matrix2))] # normal frequency for theta2
max_freq_idx = np.argmax(freq) # index of max frequency needed on graph

# check if the maximum frequency is zero
if normal_freq_theta1 == 0:
    # find the second biggest peak
    normal_freq_theta1 = freq[np.argmax(np.abs(fft_matrix1[1:]))]#normal frequency for theta1

if normal_freq_theta2 == 0:
    # find the second biggest peak
    normal_freq_theta2 = freq[np.argmax(np.abs(fft_matrix2[1:]))]#normal frequency for theta2


# getting max frequency till graph has to be made
if freq_limit != None:
    for i in range(len(freq)):
        if freq[i] > freq_limit:
            max_freq_idx = i
            break
max_freq = freq[max_freq_idx]

print(f"Plotting FFT of the both angles !")

# plotting FFT of theta1
plt.plot(freq[: max_freq_idx], np.abs(fft_matrix1)[: max_freq_idx])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title(r'FFT of the $\theta_1$ vs frequency (max at {:.4f} Hz) - {}'.format(normal_freq_theta1, file_name))
if save_fig: plt.savefig(f'{graphs_dir}\{file_name}_FFT_theta1.png')
plt.show()

# plotting FFT of theta2
plt.plot(freq[: max_freq_idx], np.abs(fft_matrix2)[: max_freq_idx])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title(r'FFT of the $\theta_2$ vs frequency (max at {:.4f} Hz) - {}'.format(normal_freq_theta2, file_name))
if save_fig: plt.savefig(f'{graphs_dir}\{file_name}_FFT_theta2.png')
plt.show()
print(f"FFT of the both angles plotted !")

if normal_freq_theta1 == normal_freq_theta2:
    print(f"Both angles have same normal frequency !")
else:
    print(f"differece in normal frequency of both angles is {np.abs(normal_freq_theta1 - normal_freq_theta2):.4f} Hz !")