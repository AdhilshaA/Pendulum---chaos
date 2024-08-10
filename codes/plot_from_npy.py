import numpy as np
import matplotlib.pyplot as plt

save_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\savefiles"
filename = "graph_trial1"
img_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\graphs"
save_fig = True

X1 = np.load(f"{save_dir}\{filename}_X1.npy")
Y1 = np.load(f"{save_dir}\{filename}_Y1.npy")

# print("plotting RK4")
plt.plot(X1,Y1[0],'-',label=r"RK4 $\theta_1(t)$")
plt.plot(X1,Y1[2],'-',label=r"RK4 $\theta_2(t)$")

X2 = np.load(f"{save_dir}\{filename}_X2.npy")
Y2 = np.load(f"{save_dir}\{filename}_Y2.npy")

# print("plotting Adaptive RK4")
plt.plot(X2,Y2[0],'-',label=r"Adaptive RK4 $\theta_1(t)$")
plt.plot(X2,Y2[2],'-',label=r"Adaptive RK4 $\theta_2(t)$")

plt.xlabel("t")
plt.ylabel(r"$\theta$ (rad )")
plt.legend()
if save_fig == True:
    plt.savefig(f"{img_dir}\{filename}.png") 
plt.show()
