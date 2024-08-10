import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from library_v2 import parse_trajectory_data,RK4
from math import radians,sin,cos,sqrt

base_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\Tracker_files" #where your trial.txt files are.

graphs_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\graphs" #where you want to save the graphs

mode = '90_0' # '90_180' or '90_0' or '170_0'
offset = 'auto' # no. of lines to23 skip in the datafile (including the header)
                # if using auto, first go open the txt file see how many lines to skip and add it after '#multi:' in same first line of data
                # for e.g. '#multi: 3' means 3 lines to skip, auto will recognize this and skip lines
                # otherwise specify value of offset manually here



g = np.float64(9.80665)
m1 = np.float64(2 * 0.3048 * 0.0381 * 0.009525 * 2710)
m2 = np.float64(0.2286 * 0.0381 * 0.009525 * 2710)
L1 = np.float64(0.3048)
L2 = np.float64(0.2286)


t1i_1 = np.float64(radians(90))
t2i_1 = np.float64(radians(179))

# using 90, 90.5, 
t1i_2 = np.float64(radians(90))
t2i_2 = np.float64(radians(181))

w1i = np.float64(0)
w2i = np.float64(0)
ti = np.float64(0)
tf = np.float64(14.5)
h = np.float64(0.0005)
tolerance = np.float64(1e-6)


n_points = int(0.2 / h) # no. of data points under analysis
le_guess = 15 # guess for le


save_fig = True # save the figures or not (True/False)
check_best_fit = False # check best fit or not (True/False)
manual_pts = n_points # no. of points to fit for best fit if check_best_fit is False

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


def exponential_func(x, le, d0):
    return d0 * np.exp(le * x)

def fit_exponential(t, diff, d0, le0):
    
    # Fit the exponential function to the data
    popt, pcov = curve_fit(exponential_func, t, diff, bounds=([-np.inf, d0], [np.inf, d0 + 1e-7])) # fit with d0 fixed
    
    # popt, pcov = curve_fit(exponential_func, t, diff, p0=[le0, d0]) # fit with d0 and le0 guesses
    

    # Calculate the chi-squared value
    residuals = diff - exponential_func(t, *popt)
    chi_squared = np.sum(residuals**2 / exponential_func(t, *popt))

    # Calculate the errors on the fitted parameters
    perr = np.sqrt(np.diag(pcov))

    # Calculate the reduced chi-squared value
    n = len(diff) # number of data points
    k = len(popt) # number of fitted parameters
    red_chi_squared = chi_squared / (n - k)

    return popt, perr, chi_squared, red_chi_squared

t1,Y1 = RK4([dt1,dw1,dt2,dw2],ti,[t1i_1,w1i,t2i_1,w2i],tf,h,print_updates = False)
t2,Y2 = RK4([dt1,dw1,dt2,dw2],ti,[t1i_2,w1i,t2i_2,w2i],tf,h,print_updates = False)

t1 = np.array(t1)
Y1 = np.array(Y1).T


theta1 =  np.degrees(Y1[:,0]) # theta 1 calculation
theta2 =  np.degrees(Y1[:,2]) # theta 2 calculation
omega1 = np.degrees(Y1[:,1])
omega2 = np.degrees(Y1[:,3])

# print(f"shape of t, theta1, theta2 : {t1.shape}, {theta1.shape}, {theta2.shape}")

theta1 = theta1[:n_points] # taking only n_points
theta2 = theta2[:n_points]
omega1 = omega1[:n_points]
omega2 = omega2[:n_points]
t1 = t1[:n_points]

# print(f'{file_name1}: theta1 = {theta1[0]}, theta2 = {theta2[0]}')

phase_space1 = np.sqrt(omega1**2 + theta1**2 + omega2**2 + theta2**2) # phase space calculation
# phase_space1 = np.sqrt(theta1**2 + theta2**2) # phase space calculation

# print(f"Calculating theta1, theta2, omega1, omega2, phase_space2 !")

t2 = np.array(t2)
Y2 = np.array(Y2).T

theta1 =  np.degrees(Y2[:,0]) # theta 1 calculation
theta2 =  np.degrees(Y2[:,2]) # theta 2 calculation
omega1 = np.degrees(Y2[:,1])
omega2 = np.degrees(Y2[:,3])

theta1 = theta1[:n_points] # taking only n_points
theta2 = theta2[:n_points]
omega1 = omega1[:n_points]
omega2 = omega2[:n_points]
t2 = t2[:n_points]

# print(f'{file_name2}: theta1 = {theta1[0]}, theta2 = {theta2[0]}')

phase_space2 = np.sqrt(omega1**2 + theta1**2 + omega2**2 + theta2**2) # phase space calculation
# phase_space2 = np.sqrt(theta1**2 + theta2**2) # phase space calculation


# print("Calculating trajectory difference !")
trajectory_diff = np.abs(phase_space1 - phase_space2)

if not np.allclose(t1, t2, atol=1e-3):
    print("t1 and t2 are not the same!")
t = t1

plt.plot(t, trajectory_diff, 'ko' ,label = 'Data points', markersize = 2)
plt.xlabel('t (s)')
plt.ylabel('Phase space difference')
plt.title('Phase space difference vs t - chaos{} - {:.2f}'.format(mode, trajectory_diff[0]))
plt.legend()
plt.yscale('log')
plt.show()
plt.close()

if check_best_fit:
    
    for k in range(4, n_points + 1):
        t_now = t[:k]
        trajectory_diff_now = trajectory_diff[:k]
        
        popt, perr, chi_squared, red_chi_squared = fit_exponential(t_now, trajectory_diff_now, trajectory_diff_now[0], le_guess)
        
        print(f'({k} points): red_chi = {red_chi_squared}, le = {popt[0]} +/- {perr[0]}, d0 = {popt[1]} +/- {perr[1]}')
    
    
    best_points = int(input("Enter the best points: "))
    t_now = t[:best_points]
    trajectory_diff_now = trajectory_diff[:best_points]
else:
    t_now = t[:manual_pts]
    trajectory_diff_now = trajectory_diff[:manual_pts]
    best_points = manual_pts

popt, perr, chi_squared, red_chi_squared = fit_exponential(t_now, trajectory_diff_now, trajectory_diff_now[0], le_guess)

# plotting all data

# print(f"Plotting all data !")
plt.figure() # create a new figure each time
plt.plot(t, trajectory_diff, 'ko' ,label = 'Data points', markersize = 2)
plt.plot(t_now, exponential_func(t_now, *popt),'r-', label = f'Best fit - {best_points} pts')
plt.xlabel('t (s)')
plt.ylabel('Phase space difference')
plt.title(r'LE_fit - chaos from {} - $d_0$ = {:.3f}, LE = ${:.3f}\pm{:.3f}$'.format(mode, trajectory_diff_now[0], popt[0], perr[0]))
plt.yscale('log')
plt.legend()
if save_fig: plt.savefig(r'{}\chaos from {}_{:.2f}_le.png'.format(graphs_dir, mode,trajectory_diff_now[0]))
plt.show()
plt.close()
