import numpy as np
import matplotlib.pyplot as plt
from library_v2 import parse_trajectory_data, RK4
from tqdm import tqdm
from math import radians,sin,cos,sqrt


base_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\Tracker_files" #where your trial.txt files are.

graphs_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\graphs" #where you want to save the graphs

mode = '90_0' # '90_180' or '90_0' or '170_0'
offset = 'auto' # no. of lines to23 skip in the datafile (including the header)
                # if using auto, first go open the txt file see how many lines to skip and add it after '#multi:' in same first line of data
                # for e.g. '#multi: 3' means 3 lines to skip, auto will recognize this and skip lines
                # otherwise specify value of offset manually here


save_fig = True # save the figures or not (True/False)

check_best_fit = False # check best fit or not (True/False)
manual_pts = 20 # no. of points to fit for best fit if check_best_fit is False

g = np.float64(9.80665)
m1 = np.float64(2 * 0.3048 * 0.0381 * 0.009525 * 2710)
m2 = np.float64(0.2286 * 0.0381 * 0.009525 * 2710)
L1 = np.float64(0.3048)
L2 = np.float64(0.2286)


initials = [(179.8,179.8),(179.9,179.9),(180,180),(180.1,180.1),(180.2,180.2)]

w1i = np.float64(0)
w2i = np.float64(0)
ti = np.float64(0)
tf = np.float64(10)
h = np.float64(0.0005)
tolerance = np.float64(1e-6)


n_points = int(1 / h) # no. of data points under analysis
le_guess = 15 # guess for le

Ja = ((1/3) * m1 * (L1 ** 2)) + (m2 * (L1 ** 2))
Jb = ((1/3) * m2 * (L2 ** 2))
Jx = ((1/2) * m2 * L1 * L2)
mu1 = (((1/2) * m1) + m2) * g * L1
mu2 = (1/2) * m2 * g * L2

def dt1(twlist,t):
    return twlist[1]

def dt2(twlist,t):
    return twlist[3]

def dw1(twlist,t):
    return ( (((Jx ** 2) / Jb) * cos(twlist[2] - twlist[0]) * sin(twlist[2] - twlist[0]) * (twlist[1] ** 2)) + (((Jx * mu2) / (Jb)) * cos(twlist[2] - twlist[0]) * sin(twlist[2])) + (Jx * sin(twlist[2] - twlist[0]) * (twlist[3] ** 2)) - (mu1 * sin(twlist[0])) ) / (Ja - (((Jx ** 2) / Jb) * (cos(twlist[2] - twlist[0]) ** 2)))


def dw2(twlist,t):
    return ( (-1 * ((Jx ** 2) / Ja) * cos(twlist[2] - twlist[0]) * sin(twlist[2] - twlist[0]) * (twlist[3] ** 2)) + (((Jx * mu1) / (Ja)) * cos(twlist[2] - twlist[0]) * sin(twlist[0])) - (Jx * sin(twlist[2] - twlist[0]) * (twlist[1] ** 2)) - (mu2 * sin(twlist[2])) ) / (Jb - (((Jx ** 2) / Ja) * (cos(twlist[2] - twlist[0]) ** 2)))

for i in range(len(initials)):
    
    t1i,t2i = initials[i]
    
    t1,Y1 = RK4([dt1,dw1,dt2,dw2],ti,[t1i,w1i,t2i,w2i],tf,h,print_updates = False)

    t1 = np.array(t1)
    Y1 = np.array(Y1).T


    theta1 =  np.degrees(Y1[:,0]) # theta 1 calculation
    theta2 =  np.degrees(Y1[:,2]) # theta 2 calculation
    omega1 = np.degrees(Y1[:,1])
    omega2 = np.degrees(Y1[:,3])

    theta1 = theta1[:n_points] # taking only n_points
    theta2 = theta2[:n_points]
    omega1 = omega1[:n_points]
    omega2 = omega2[:n_points]
    t1 = t1[:n_points]

    plt.plot(t1, theta1, '--o', markersize = 1, label = r'({},{})-$\theta_1$'.format(t1i,t2i))
    plt.plot(t1, theta2, '--o', markersize = 1, label = r'({},{})-$\theta_2$'.format(t1i,t2i))
    plt.xlabel('time (s)')
    plt.ylabel(r'angle, $\theta$ (deg)')
plt.title('Trajectory vs t - from ({})'.format(mode))
plt.legend()
if save_fig:
    plt.savefig(f'{graphs_dir}\chaos_traj_{mode}.png')
plt.show()