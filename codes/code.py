# DIY code and output

import library as lib
import math as m
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#------INPUTS AND MEANINGS----------_#
# g : acceleration due to gravity (m/s)
# m1 : mass (kg) of first bob from top
# m2 : mass (kg) of second bob from top 
# L1 : Length (m) of massless rigid rod connecting from the fixed point to first bob from top
# L2 : Length (m) of massless rigid rod connecting from the first bob to second bob
# ti : initial time (s)
# t1i : initial theta1 (rad) at initial time
# t2i : initial theta2 (rad) at initial time
# w1i : initial omega1 (rad/s) at initial time
# w2i : initial omega2 (rad/s) at initial time
# tf : time (s) till the solution is needed
# hi  : initial step size for the Adaptive step-size RK4
# h : fixed step-size for simple RK4
# tolerance : tolerance level for adaptive step-size RK4

#---------VALUES USED------------------#
# uncomment and alter values if needed
g = 9.8
m1 = 2
m2 = 2
L1 = 1
L2 = 1
t1i = m.pi-0.01
t2i = 0
w1i = 0
w2i = 0
ti = 0
tf = 10
hi  = 0.001
h = 0.025
tolerance = 1e-6


# twlist is in format [t1, w1, t2, w2], these are the first order ODEs to be solved
def dt1(twlist,t):
    return twlist[1]

def dt2(twlist,t):
    return twlist[3]

def dw1(twlist,t):
    return ((-1*g*((2*m1) + m2)*m.sin(twlist[0])) - (m2*g*m.sin(twlist[0] - (2*twlist[2]))) - (2*m.sin(twlist[0]-twlist[2])*m2*(((twlist[3]**2)*L2) + ((twlist[1]**2)*L1*m.cos(twlist[0]-twlist[2])))))/(L1*((2*m1)+m2-(m2*m.cos((2*twlist[0]) - (2*twlist[2])))))

def dw2(twlist,t):
    return (2*m.sin(twlist[0]-twlist[2])* (((twlist[1]**2)*L1*(m1+m2)) + (g*(m1 + m2)*m.cos(twlist[0])) + ((twlist[3]**2)*L2*m2*m.cos(twlist[0]-twlist[2]))))/(L2*((2*m1)+m2-(m2*m.cos((2*twlist[0]) - (2*twlist[2])))))


### ---------------simple RK4 application ----------------####
X1,Y1 = lib.RK4([dt1,dw1,dt2,dw2],ti,[t1i,w1i,t2i,w2i],tf,h)
# uncomment to see the data for simple RK4
# lib.print_coltable({"Sl. NO.":[i for i in range(len(X1))],"t":X1,"t1":Y1[0],"w1":Y1[1],"t2":Y1[2],"w2":Y1[3],})

### ---------------adaptive step size RK4 application----------------####
X2,Y2,counts = lib.ASRK4([dt1,dw1,dt2,dw2],ti,[t1i,w1i,t2i,w2i],tf,hi,tolerance,errprioirity=[0,2])
# uncomment the below line to see the data for ASRK4
# lib.print_coltable({"Sl. NO.":[i for i in range(len(X2))],"t":X2,"t1":Y2[0],"w1":Y2[1],"t2":Y2[2],"w2":Y2[3],})

# ### ---------------plotting simple RK4 ----------------####
# # plt.vlines(X1,min(min(Y1[0]),min(Y1[2])),max(max(Y1[0]),max(Y1[2])),color = 'gray',lw = 1)
# plt.plot(X1,Y1[0],'-',label="RK4 t1(t)")
# # plt.plot(X1,Y1[1],'-',label="RK4 w1")
# # plt.plot(X1,Y1[3],'-',label="RK4 w2")
# plt.plot(X1,Y1[2],'-',label="RK4 t2(t)")

# ### ---------------plotting adaptive step size RK4 ----------------####
# # plt.vlines(X2,min(min(Y2[0]),min(Y2[2])),max(max(Y2[0]),max(Y2[2])),color = 'gray',lw = 1)
# plt.plot(X2,Y2[0],'-',ms=2,label="ASRK4 t1(t)")
# # plt.plot(X2,Y2[1],'-',ms=2,label="ASRK4 w1(t)")
# # plt.plot(X2,Y2[3],'-',ms=2,label="ASRK4 w2(t)")
# plt.plot(X2,Y2[2],'-',ms=2,label="ASRK4 t2(t)")

# #-------------------- labels and showing graph----------------####
# plt.xlabel("time, t(s)")
# plt.ylabel("angles, t1 and t2 (rad)")
# plt.legend(fontsize = 8)
# plt.title("t1,t2 vs t graph (10 sec)")
# plt.show()

# print("The different graphs attained by using different plotting parameter are given in the folder starting with names 'graph_<name>.png'")
# print("Check report for more specific details, the current graph is set to show theta1 and theta2 for 10 seconds for both programs.")

# ######-----  OUTPUT  -----######
# '''
# (The different graphs attained by using different plotting parameter are given in the folder starting with names "graph_<name>.png")
# (Check report for more specific details, the current graph is set to show theta1 and theta2 for 10 seconds for both programs.)
# '''

print(f"shape of X1: (1,{len(X1)})")
print(f"shape of Y1: ({len(Y1)},{len(Y1[0])})")
print(f"shape of X2: (1,{len(X2)})")
print(f"shape of Y2: ({len(Y2)},{len(Y2[0])})")

print(f"type of X1: {type(X1[0])}")
print(f"type of Y1s: {[type(y[0]) for y in Y1]}")
print(f"type of X2: {type(X2[0])}")
print(f"type of Y2s: {[type(y[0]) for y in Y2]}")