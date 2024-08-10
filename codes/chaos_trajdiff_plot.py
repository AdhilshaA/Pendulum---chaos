import numpy as np
import matplotlib.pyplot as plt
from library_v2 import parse_trajectory_data
from tqdm import tqdm

base_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\Tracker_files" #where your trial.txt files are.
late = []
used = []

# graphs_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\graphs\chaos_90_180_trials" #where you want to save the graphs
# start_trial = 31
# end_trial = 41
# used = list(range(start_trial, end_trial + 1))
# used = [31,38,39,40,41]
# late = [31,32,37,41]
# mode = '90_180'

# graphs_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\graphs\chaos_150_0_trials" #where you want to save the graphs
# start_trial = 1
# end_trial = 7
# used = [1, 2, 3, 4, 5, 6, 7]
# late = [1]
# mode = '150_0'

# graphs_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\graphs\chaos_90_0_trials" #where you want to save the graphs
# start_trial = 11
# end_trial = 15
# used = [11, 13, 14, 15]
# late = [12,13,14,15]
# mode = '90_0'

graphs_dir = "H:\PERSONAL\PROJECTS\Programming\Openlab_7thsem\graphs\chaos_170_0_trials" #where you want to save the graphs
start_trial = 51
end_trial = 58
used = [52, 54, 55, 57, 58]
late = [51,53]
mode = '170_0'


file_name = 'chaos' # name of the datafile(no .txt)
offset = 'auto' # no. of lines to23 skip in the datafile (including the header)
                # if using auto, first go open the txt file see how many lines to skip and add it after '#multi:' in same first line of data
                # for e.g. '#multi: 3' means 3 lines to skip, auto will recognize this and skip lines
                # otherwise specify value of offset manually here

n_points = 30 # no. of data points under analysis

save_fig = True # save the figures or not (True/False)


for i in tqdm(range(start_trial, end_trial + 1)):
    if i not in used:
        continue
    file_name1 = f'{file_name}{i}'
    for j in range(i + 1, end_trial + 1):
        
        file_name2 = f'{file_name}{j}'
        if j not in used or int(i) == int(j):
            continue
        
        # Parse the datafile into a matrix
        data1 = parse_trajectory_data(f'{base_dir}\{file_name1}.txt', offset = offset, updates = False)
        data2 = parse_trajectory_data(f'{base_dir}\{file_name2}.txt', offset = offset, updates = False)
        # print(f"Data parsed !")


        # # creating t, theta1, theta2 matrix
        # print(f"Calculating theta1, theta2, omega1, omega2, phase_space1 !")

        t1 = data1[:,0] #adding time column
        t1 = t1 - t1[0]
        theta1 =  np.degrees(np.arctan2(data1[:,1] , (- data1[:,2]))) # theta 1 calculation
        theta2 =  np.degrees(np.arctan2((data1[:,3]- data1[:,1]) , (data1[:,2]- data1[:,4]))) # theta 2 calculation
        # N.B. These theta1 and theta2 are in (-180, 180) range

        # converting theta1 and theta2 to (-inf, inf) range
        if mode == '90_180':
            for k in range(len(theta1)):
                if theta1[k] < 0:
                    theta1[k] += 360
                if theta2[k] < 0:
                    theta2[k] += 360
                    
        add_temp = 0
        temp_next = theta1[0]
        for k in range(len(theta1) - 1):
            if theta1[k+1] - temp_next > 170:
                add_temp -= 360
            elif theta1[k+1] - temp_next < -170:
                add_temp += 360
            temp_next = theta1[k+1]
            theta1[k+1] += add_temp

        add_temp = 0
        temp_next = theta2[0]
        for k in range(len(theta2) - 1):
            if theta2[k+1] - temp_next > 170:
                add_temp -= 360
            elif theta2[k+1] - temp_next < -170:
                add_temp += 360
            temp_next = theta2[k+1]
            theta2[k+1] += add_temp

        omega1 = np.gradient(theta1, t1) # angular velocity calculation
        omega2 = np.gradient(theta2, t1)

        omega1[0] = 0 # setting initial velocity to be 0
        omega2[0] = 0
        
        if i in late:
            theta1 = theta1[1:n_points + 1] # taking only n_points
            theta2 = theta2[1:n_points + 1]
            omega1 = omega1[1:n_points + 1]
            omega2 = omega2[1:n_points + 1]
            t1 = t1[1:n_points + 1] - t1[1]
        else:
            theta1 = theta1[:n_points] # taking only n_points
            theta2 = theta2[:n_points]
            omega1 = omega1[:n_points]
            omega2 = omega2[:n_points]
            t1 = t1[:n_points]

        # print(f'{file_name1}: theta1 = {theta1[0]}, theta2 = {theta2[0]}')

        phase_space1 = np.sqrt(omega1**2 + theta1**2 + omega2**2 + theta2**2) # phase space calculation
        # phase_space1 = np.sqrt(theta1**2 + theta2**2) # phase space calculation

        # print(f"Calculating theta1, theta2, omega1, omega2, phase_space2 !")

        t2 = data2[:,0] #adding time column
        t2 = t2 - t2[0]
        theta1 =  np.degrees(np.arctan2(data2[:,1] , (- data2[:,2]))) # theta 1 calculation
        theta2 =  np.degrees(np.arctan2((data2[:,3]- data2[:,1]) , (data2[:,2]- data2[:,4]))) # theta 2 calculation

        # converting theta1 and theta2 to (-inf, inf) range
        if mode == '90_180':
            for k in range(len(theta1)):
                if theta1[k] < 0:
                    theta1[k] += 360
                if theta2[k] < 0:
                    theta2[k] += 360
        add_temp = 0
        temp_next = theta1[0]
        for k in range(len(theta1) - 1):
            if theta1[k+1] - temp_next > 170:
                add_temp -= 360
            elif theta1[k+1] - temp_next < -170:
                add_temp += 360
            temp_next = theta1[k+1]
            theta1[k+1] += add_temp
            
        add_temp = 0
        temp_next = theta2[0]
        for k in range(len(theta2) - 1):
            if theta2[k+1] - temp_next > 170:
                add_temp -= 360
            elif theta2[k+1] - temp_next < -170:
                add_temp += 360
            temp_next = theta2[k+1]
            theta2[k+1] += add_temp

        omega1 = np.gradient(theta1, t2) # angular velocity calculation
        omega2 = np.gradient(theta2, t2)

        omega1[0] = 0 # setting initial velocity to be 0
        omega2[0] = 0

        if j in late:
            theta1 = theta1[1:n_points + 1]
            theta2 = theta2[1:n_points + 1]
            omega1 = omega1[1:n_points + 1]
            omega2 = omega2[1:n_points + 1]
            t2 = t2[1:n_points + 1] - t2[1]

        else:
            theta1 = theta1[:n_points] # taking only n_points
            theta2 = theta2[:n_points]
            omega1 = omega1[:n_points]
            omega2 = omega2[:n_points]
            t2 = t2[:n_points]

        # print(f'{file_name2}: theta1 = {theta1[0]}, theta2 = {theta2[0]}')

        phase_space2 = np.sqrt(omega1**2 + theta1**2 + omega2**2 + theta2**2) # phase space calculation
        # phase_space2 = np.sqrt(theta1**2 + theta2**2) # phase space calculation


        # print("Calculating trajectory difference !")
        trajectory_diff = np.abs(phase_space2 - phase_space1)

        if not np.allclose(t1, t2, atol=1e-3):
            print("t1 and t2 are not the same!")
        t = t1

        plt.figure() # create a new figure each time
        plt.plot(t, trajectory_diff, 'ko' ,label = 'Data points')
        plt.xlabel('t (s)')
        plt.ylabel('Phase space difference')
        plt.title('Phase space difference vs t - chaos{}_{}'.format(file_name1[5:],file_name2[5:]))
        plt.legend()
        plt.yscale('log')
        if save_fig:
            plt.savefig(f'{graphs_dir}\chaos_trajdiff_{file_name1[5:]}_{file_name2[5:]}.png')
        plt.close()