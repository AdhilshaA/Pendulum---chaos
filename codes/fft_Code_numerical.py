import numpy as np
from scipy.fft import fft
from scipy.fft import fftfreq
import matplotlib.pyplot as plt
from library_v2 import parse_trajectory_data
from scipy.signal.windows import blackman, hann
import library_v2 as lib
from math import radians,sin,cos,sqrt


graphs_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\graphs\sim_sol" # where you want to save the graphs
save_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\savefiles"
file_name = "num2"

g = np.float64(9.80665) # m/s^2 acceleration due to gravity
m1 = np.float64(2 * 0.3048 * 0.0381 * 0.009525 * 2710) # kg - mass of upper rod
m2 = np.float64(0.2286 * 0.0381 * 0.009525 * 2710) # kg - mass of lower rod
L1 = np.float64(0.3048) # m - length of upper rod
L2 = np.float64(0.2286) # m - length of lower rod
B1 = np.float64(0.0381)
B2 = np.float64(0.0381)
t1i = np.float64(radians(2))

NM1_ratio = 3.5 / 2
NM2_ratio = -6.8 / 2

t2i = t1i * NM1_ratio # normal mode1 by trial and error
mode = "NM1"
# t2i = t1i * NM2_ratio # normal mode2 by trial and error
# mode = "NM2"

#scaled NM modes



# t2i = sqrt(2)*t1i # normal mode1 for equal length and stuff
# t2i = -sqrt(2)*t1i # normal mode2 for equal length and stuff

w1i = np.float64(0)
w2i = np.float64(0)
ti = np.float64(0)
tf = np.float64(14.5)
h = np.float64(0.0005)
tolerance = np.float64(1e-6)

mu = 1 + m1/m2

I1 = (1/12) * m1 * (L1 ** 2) # moment of inertia of both arms in upper half
I2 = (1/12) * m2 * (L2 ** 2) # moment of inertia of only arm in lower half

# I1 = (m1 * ((L1 ** 2) + (B2 ** 2)) / 12)
# I2 = (m2 * ((L2 ** 2) + (B2 ** 2)) / 12)

Ja = ((1/4) * m1 * (L1 ** 2)) + (m2 * (L1 ** 2) + I1)
Jb = ((1/4) * m2 * (L2 ** 2) + I2)
Jx = ((1/2) * m2 * L1 * L2)
mu1 = (((1/2) * m1) + m2) * g * L1
mu2 = (1/2) * m2 * g * L2

save_resultnpy = False
print_updates = False

save_fig = False # save the figures or not (True/False)
freq_limit = 10 # in Hz for limit on graph   
window = True # apply window or not (True/False)
# which mode we are processing

# twlist is in format [t1, w1, t2, w2], these are the first order ODEs to be solved
def dt1(twlist,t):
    return twlist[1]

def dt2(twlist,t):
    return twlist[3]

#mphysics formulae for the second order ODEs
def dw1(twlist,t):
    return ((-1*g*((2*m1) + m2)*sin(twlist[0])) - (m2*g*sin(twlist[0] - (2*twlist[2]))) - (2*sin(twlist[0]-twlist[2])*m2*(((twlist[3]**2)*L2) + ((twlist[1]**2)*L1*cos(twlist[0]-twlist[2])))))/(L1*((2*m1)+m2-(m2*cos((2*twlist[0]) - (2*twlist[2])))))

def dw2(twlist,t):
    return (2*sin(twlist[0]-twlist[2])* (((twlist[1]**2)*L1*(m1+m2)) + (g*(m1 + m2)*cos(twlist[0])) + ((twlist[3]**2)*L2*m2*cos(twlist[0]-twlist[2]))))/(L2*((2*m1)+m2-(m2*cos((2*twlist[0]) - (2*twlist[2])))))

# alterate equations from paper
# def dw1(twlist,t):
#     return ( ( g * ( ( sin(twlist[2]) * cos(twlist[0] - twlist[2])) - ( mu * sin(twlist[0]) ) ) ) - ( ( ( L2 * (twlist[3] ** 2) ) + ( L1 * (twlist[1] ** 2) * cos(twlist[0] - twlist[2]) ) ) * ( sin(twlist[0] - twlist[2]) ) ) ) / ( L1 * (mu - ((cos(twlist[0] - twlist[2])) ** 2) ) )

# def dw2(twlist,t):
#     return ( ( g * mu * ( ( sin(twlist[0]) * cos(twlist[0] - twlist[2])) - (  sin(twlist[2]) ) ) ) - ( ( ( mu * L1 * (twlist[1] ** 2) ) + ( L2 * (twlist[3] ** 2) * cos(twlist[0] - twlist[2]) ) ) * ( sin(twlist[0] - twlist[2]) ) ) ) / ( L2 * (mu - ((cos(twlist[0] - twlist[2])) ** 2) ) )

# my own equations
def dw1(twlist,t):
    return ( (((Jx ** 2) / Jb) * cos(twlist[2] - twlist[0]) * sin(twlist[2] - twlist[0]) * (twlist[1] ** 2)) + (((Jx * mu2) / (Jb)) * cos(twlist[2] - twlist[0]) * sin(twlist[2])) + (Jx * sin(twlist[2] - twlist[0]) * (twlist[3] ** 2)) - (mu1 * sin(twlist[0])) ) / (Ja - (((Jx ** 2) / Jb) * (cos(twlist[2] - twlist[0]) ** 2)))


def dw2(twlist,t):
    return ( (-1 * ((Jx ** 2) / Ja) * cos(twlist[2] - twlist[0]) * sin(twlist[2] - twlist[0]) * (twlist[3] ** 2)) + (((Jx * mu1) / (Ja)) * cos(twlist[2] - twlist[0]) * sin(twlist[0])) - (Jx * sin(twlist[2] - twlist[0]) * (twlist[1] ** 2)) - (mu2 * sin(twlist[2])) ) / (Jb - (((Jx ** 2) / Ja) * (cos(twlist[2] - twlist[0]) ** 2)))


print("Running simple RK4")
t,Y1 = lib.RK4([dt1,dw1,dt2,dw2],ti,[t1i,w1i,t2i,w2i],tf,h,print_updates)

t = np.array(t)
Y1 = np.array(Y1).T

theta1 =  np.degrees(Y1[:,0]) # theta 1 calculation
theta2 =  np.degrees(Y1[:,2]) # theta 2 calculation
print(f"t, theta1, theta2 matrix created !")

# plotting t vs theta1 and t vs theta2
print(f"Plotting t vs theta1 and t vs theta2 !")
plt.plot(t, theta1, label = r'$\theta_1$')
plt.plot(t, theta2, label = r'$\theta_2$')
plt.xlabel('t (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta_1$ and $\theta_2$ vs t - {}'.format(file_name+'_'+mode))
plt.legend()
if save_fig: plt.savefig(f'{graphs_dir}\{file_name}_{mode}_trajectory.png')
plt.show()
print(f"t vs theta1 and t vs theta2 plotted !")


if window:
    N = len(theta1) # number of samples
    window = blackman(N) # window function
    
    print(f"Applying FFT on the both angles with window !")
    fft_matrix1 = fft(theta1 * window)
    fft_matrix2 = fft(theta2 * window)
    print(f"FFT applied on the both angles with window !")
    
    freq = fftfreq(t.shape[0], h) # getting related frequencies for FFT
    
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


    freq = fftfreq(t.shape[0], h) # getting related frequencies for FFT
    

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
plt.title(r'FFT of the $\theta_1$ (max at {:.4f} Hz) - {}'.format(normal_freq_theta1, file_name+'_'+mode))
if save_fig: plt.savefig(f'{graphs_dir}\{file_name}_{mode}_FFT_theta1.png')
plt.show()

# plotting FFT of theta2
plt.plot(freq[: max_freq_idx], np.abs(fft_matrix2)[: max_freq_idx])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title(r'FFT of the $\theta_2$ (max at {:.4f} Hz) - {}'.format(normal_freq_theta2, file_name+'_'+mode))
if save_fig: plt.savefig(f'{graphs_dir}\{file_name}_{mode}_FFT_theta2.png')
plt.show()
print(f"FFT of the both angles plotted !")

if normal_freq_theta1 == normal_freq_theta2:
    print(f"Both angles have same normal frequency !")
else:
    print(f"differece in normal frequency of both angles is {np.abs(normal_freq_theta1 - normal_freq_theta2):.4f} Hz !")