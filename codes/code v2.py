# DIY code and output

import library_v2 as lib
from math import radians,sin,cos,sqrt
import matplotlib.pyplot as plt
import numpy as np
import time
# your code here



#------INPUTS AND MEANINGS----------_#
# g   : acceleration due to gravity (m/s)
# m1  : mass (kg) of first bob from top
# m2  : mass (kg) of second bob from top 
# L1  : Length (m) of massless rigid rod connecting from the fixed point to first bob from top
# L2  : Length (m) of massless rigid rod connecting from the first bob to second bob
# ti  : initial time (s)
# t1i : initial theta1 (rad) at initial time
# t2i : initial theta2 (rad) at initial time
# w1i : initial omega1 (rad/s) at initial time
# w2i : initial omega2 (rad/s) at initial time
# tf  : time (s) till the solution is needed
# hi  : initial step size for the Adaptive step-size RK4
# h   : fixed step-size for simple RK4
# tolerance : tolerance level for adaptive step-size RK4

#---------VALUES USED------------------#
# Uncomment and alter values if needed
g = np.float64(9.80665)
m1 = np.float64(2 * 0.3048 * 0.0381 * 0.009525 * 2710)
m2 = np.float64(0.2286 * 0.0381 * 0.009525 * 2710)
L1 = np.float64(0.3048)
L2 = np.float64(0.2286)
t1i = np.float64(radians(2))

NM1_ratio = 3.5 / 2
NM2_ratio = -6.8 / 2

t2i = t1i * NM1_ratio # normal mode1 by trial and error
t2i = t1i * NM2_ratio # normal mode2 by trial and error

#scaled NM modes



# t2i = sqrt(2)*t1i # normal mode1 for equal length and stuff
# t2i = -sqrt(2)*t1i # normal mode2 for equal length and stuff

w1i = np.float64(0)
w2i = np.float64(0)
ti = np.float64(0)
tf = np.float64(24)
h = np.float64(0.0005)
tolerance = np.float64(1e-6)

mu = 1 + m1/m2

Ja = ((1/3) * m1 * (L1 ** 2)) + (m2 * (L1 ** 2))
Jb = ((1/3) * m2 * (L2 ** 2))
Jx = ((1/2) * m2 * L1 * L2)
mu1 = (((1/2) * m1) + m2) * g * L1
mu2 = (1/2) * m2 * g * L2

# figure savefile directory
save_resultnpy = False
plot_img = True
save_img = False
print_updates = True
img_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\graphs\sim_sol"
file_name = "graph_trialoldeqn"

# twlist is in format [t1, w1, t2, w2], these are the first order ODEs to be solved
def dt1(twlist,t):
    return twlist[1]

def dt2(twlist,t):
    return twlist[3]

#mphysics formulae for the second order ODEs
# def dw1(twlist,t):
#     return ((-1*g*((2*m1) + m2)*sin(twlist[0])) - (m2*g*sin(twlist[0] - (2*twlist[2]))) - (2*sin(twlist[0]-twlist[2])*m2*(((twlist[3]**2)*L2) + ((twlist[1]**2)*L1*cos(twlist[0]-twlist[2])))))/(L1*((2*m1)+m2-(m2*cos((2*twlist[0]) - (2*twlist[2])))))

# def dw2(twlist,t):
#     return (2*sin(twlist[0]-twlist[2])* (((twlist[1]**2)*L1*(m1+m2)) + (g*(m1 + m2)*cos(twlist[0])) + ((twlist[3]**2)*L2*m2*cos(twlist[0]-twlist[2]))))/(L2*((2*m1)+m2-(m2*cos((2*twlist[0]) - (2*twlist[2])))))

# my own equations
def dw1(t1,t2,w1,w2,t):
    return ( (((Jx ** 2) / Jb) * cos(t2 - t1) * sin(t2 - t1) * (w1 ** 2)) + (((Jx * mu2) / (Jb)) * cos(t2 - t1) * sin(t1)) + (Jx * sin(t2 - t1) * (w2 ** 2)) - (mu1 * sin(t1)) ) / (Ja - (((Jx ** 2) / Jb) * (cos(t2 - t1) ** 2)))


def dw2(t1,t2,w1,w2,t):
    return ( (-1 * ((Jx ** 2) / Ja) * cos(t2 - t1) * sin(t2 - t1) * (w2 ** 2)) + (((Jx * mu1) / (Ja)) * cos(t2 - t1) * sin(t1)) - (Jx * sin(t2 - t1) * (w1 ** 2)) - (mu2 * sin(t2)) ) / (Jb - (((Jx ** 2) / Ja) * (cos(t2 - t1) ** 2)))

#reference papers formula
# def dw1(twlist,t):
#     return ( ( g * ( ( sin(twlist[2]) * cos(twlist[0] - twlist[2])) - ( mu * sin(twlist[0]) ) ) ) - ( ( ( L2 * (twlist[3] ** 2) ) + ( L1 * (twlist[1] ** 2) * cos(twlist[0] - twlist[2]) ) ) * ( sin(twlist[0] - twlist[2]) ) ) ) / ( L1 * (mu - ((cos(twlist[0] - twlist[2])) ** 2) ) )

# def dw2(twlist,t):
#     return ( ( g * mu * ( ( sin(twlist[0]) * cos(twlist[0] - twlist[2])) - (  sin(twlist[2]) ) ) ) - ( ( ( mu * L1 * (twlist[1] ** 2) ) + ( L2 * (twlist[3] ** 2) * cos(twlist[0] - twlist[2]) ) ) * ( sin(twlist[0] - twlist[2]) ) ) ) / ( L2 * (mu - ((cos(twlist[0] - twlist[2])) ** 2) ) )


start_time = time.time()
### ---------------simple RK4 application ----------------####
print("Running simple RK4")
X1,Y1 = lib.RK4([dt1,dw1,dt2,dw2],ti,[t1i,w1i,t2i,w2i],tf,h,print_updates)

# uncomment to see the data for simple RK4
# lib.print_coltable({"Sl. NO.":[i for i in range(len(X1))],"t":X1,"t1":Y1[0],"w1":Y1[1],"t2":Y1[2],"w2":Y1[3],})
run_time1 = time.time()

if plot_img == True:
    ### ---------------plotting simple RK4 ----------------####
    # plt.vlines(X1,min(min(Y1[0]),min(Y1[2])),max(max(Y1[0]),max(Y1[2])),color = 'gray',lw = 1)
    plt.plot(X1,Y1[0],'-',label=r"RK4 $\theta_1(t)$")
    # plt.plot(X1,Y1[1],'-',label="RK4 w1")
    # plt.plot(X1,Y1[3],'-',label="RK4 w2")
    plt.plot(X1,Y1[2],'-',label=r"RK4 $\theta_2(t)$")
plot_time1 = time.time()
del X1,Y1
save_time1 = time.time()


### ---------------adaptive step size RK4 application----------------####
print("Running adaptive step size RK4")
X2,Y2,counts = lib.ASRK4([dt1,dw1,dt2,dw2],ti,[t1i,w1i,t2i,w2i],tf,h,tolerance,errprioirity=[0,1,2,3],print_updates=print_updates)
# uncomment the below line to see the data for ASRK4
# lib.print_coltable({"Sl. NO.":[i for i in range(len(X2))],"t":X2,"t1":Y2[0],"w1":Y2[1],"t2":Y2[2],"w2":Y2[3],})
run_time2 = time.time()


if plot_img == True:
    ### ---------------plotting adaptive step size RK4 ----------------####
    # plt.vlines(X2,min(min(Y2[0]),min(Y2[2])),max(max(Y2[0]),max(Y2[2])),color = 'gray',lw = 1)
    plt.plot(X2,Y2[0],'-',ms=2,label=r"ASRK4 $\theta_1(t)$")
    # plt.plot(X2,Y2[1],'-',ms=2,label="ASRK4 w1(t)")
    # plt.plot(X2,Y2[3],'-',ms=2,label="ASRK4 w2(t)")
    plt.plot(X2,Y2[2],'-',ms=2,label=r"ASRK4 $\theta_2(t)$")

    #-------------------- labels and showing graph----------------####
    plt.xlabel(r"time, $t(s)$")
    plt.ylabel(r"angles, $\theta_1$ and $\theta_2$ (rad)")
    plt.legend(fontsize = 8)
    plt.title(r"$\theta_1,\theta_2$ vs $t$ graph - (${}$ stepsize)".format(h))
    if save_img == True:
        plt.savefig(f"{img_dir}\{file_name}.png")
    plt.show()
    
plot_time2 = time.time()
print("\n\tCalculations complete !")


del X2,Y2

save_time2 = time.time()

print("{} stepsize RK4 completed".format(h))
print("Simple RK4 run time: {}s".format(run_time1 - start_time))
print("Simple RK4 plot time: {}s".format(plot_time1 - run_time1))
print("Simple RK4 save time: {}s".format(save_time1 - plot_time1))
print("ASRK4 run time: {}s".format(run_time2 - save_time1))
print("ASRK4 plot time: {}s".format(plot_time2 - run_time2))
print("ASRK4 save time: {}s".format(save_time2 - plot_time2))
print("Total run time: {}s".format(save_time2 - start_time))


# print(f"shape of X1: (1,{len(X1)})")
# print(f"shape of Y1: ({len(Y1)},{len(Y1[0])})")
# print(f"shape of X2: (1,{len(X2)})")
# print(f"shape of Y2: ({len(Y2)},{len(Y2[0])})")

# print(f"type of X1: {type(X1[0])}")
# print(f"type of Y1s: {[type(y[0]) for y in Y1]}")
# print(f"type of X2: {type(X2[0])}")
# print(f"type of Y2s: {[type(y[0]) for y in Y2]}")

