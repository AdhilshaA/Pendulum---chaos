import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


# Define the exponential function to fit
def exponential_func(x, le, d0):
    return d0 * np.exp(le * x)
# def exponential_func(x, le, d0, c):
#     return d0 * np.exp(le * x) + c

def fit_exponential(t, diff, d0, le0):
    
    # Fit the exponential function to the data
    # popt, pcov = curve_fit(exponential_func, t, diff) # complete fit
    popt, pcov = curve_fit(exponential_func, t, diff, bounds=([-np.inf, d0], [np.inf, d0 + 1e-7])) # fit with d0 fixed
    # popt, pcov = curve_fit(exponential_func, t, diff, p0=[le0, d0]) # fit with d0 and le0 guesses
    # popt, pcov = curve_fit(exponential_func, t, diff, bounds=([-np.inf, d0,-np.inf ], [np.inf, d0 + 1e-7,np.inf]))
    

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

# trial data
t = (np.arange(1,15) - 1) / 25
print(t)
diff1 = [1,14.5,15.3,15.8,16.9,23.9,10.7,17.3,19,20.1,16.2,9.9,8.4,4.5]
diff2 = [2.4,9.6,10.7,16.3,20.7,26.5,36.3,51.1,112.6,44,32.6,13.9,14.2,18]
diff3 = [2.5,6.2,10.1,6.9,16.3,21.2,29.4,41.2,104.3,63.7,47.5,14.9,21.8,20.9]
diff4 = [7,11,4.5,8.9,24.6,34.2,3.2,3.9,18.9]
diff5 = [6.9,23.3,14.8,24.1,8.1,13.1,9.7,14.1,24.6]
diff6 = [6.8,15.4,7.1,11.8,27.2,30.9,29.8,41.3,108.1]
diff7 = [5.2,7.1,7.6,7.6,9.1,41.6,89.3,50.8,21.6]

le_found = [3.50,9.73,9.07,7.88,-0.248,6.3,4.38]
len_found = [11,9,8,6,8,8,4]

diffs = [diff1, diff2, diff3, diff4, diff5, diff6, diff7]
d0s = [diff[0] for diff in diffs]
# d0s = [1 for i in range(len(d0s))]
check_best = False


labels = ['trial1_2', 'trial2_3', 'trial1_3', 'trial1_4', 'trial2_4', 'trail3_4', 'trial1_5']

for i, diff in enumerate(diffs):
    
    if i != 2:
        continue

    if check_best == True:
        # Define the range of points to consider
        max_points = len(diffs[i])
        min_points = 4
        # Initialize variables to store the best fit
        best_fit = None
        best_points = None
        best_chi_squared = np.inf

        # Loop over all possible combinations of points
        for n_points in range(min_points, max_points+1):
            t_now = t[:n_points]
            diff_now = diffs[i][:n_points]
            
            # resetting time
            t_now = t_now - t_now[0]

            # Fit the exponential function to the data
            popt, perr, chi_squared, red_chi_squared = fit_exponential(t_now, diff_now, diff_now[0], le_found[i])
            # popt, perr, chi_squared, red_chi_squared = fit_exponential(t_now, diff_now, d0s[i], le_found[i])
            
            
            print(f"n_points: {n_points}, red_chi_squared: {red_chi_squared}")
            
            # Check if this fit is the best so far
            if red_chi_squared < best_chi_squared:
                best_fit = popt
                best_points = (0, n_points)
                best_chi_squared = red_chi_squared
        
        j_start, j_end = best_points
        t_now = t[j_start:j_end]
        diff_now = diffs[i][j_start:j_end]
        popt = best_fit
    else:
        start = 0
        j_start, j_end = (start,start + len_found[i])
        t_now = t[j_start:j_end]
        diff_now = diffs[i][j_start:j_end]
        
        t_now = t[j_start:j_end] - t[j_start] # resetting time
        
        # popt, perr, chi_squared, red_chi_squared = fit_exponential(t_now, diff_now, diff_now[0], le_found[i])
        popt, perr, chi_squared, red_chi_squared = fit_exponential(t_now, diff_now, d0s[i], le_found[i])
        
        
        print(f"n_points: {len_found[i]}, red_chi_squared: {red_chi_squared}")
        best_chi_squared = red_chi_squared

    # print the results
    print("For diff{}:".format(i+1))
    print("Fitted parameters: le = {:.3f} +/- {:.3f}, d0 = {:.3f} +/- {:.3f}".format(popt[0], perr[0], popt[1], perr[1]))
    # print("Fitted parameters: le = {:.3f} +/- {:.3f}, d0 = {:.3f} +/- {:.3f}, c = {:.3f} +/- {:.3f}".format(popt[0], perr[0], popt[1], perr[1], popt[2], perr[2]))
    
    print("Chi-squared value: {:.3f}".format(best_chi_squared))
    print("Points used: {}-{}".format(j_start+1, j_end))

    # plot the best fit
    plt.plot(t[:len(diffs[i])], diffs[i], 'o', label=f'data{i}')
    
    # plt.plot(t_now, exponential_func(t_now, *popt), '-', label=f'fit{i}')
    # plt.plot(t_now, exponential_func(t_now, le_found[i], diff_now[0]), '-', label=f'papers{i}')
    
    plt.plot(t_now + t[j_start], exponential_func(t_now, *popt), '-', label=f'fit{i}')
    plt.plot(t_now + t[j_start], exponential_func(t_now, le_found[i], diff_now[0]), '-', label=f'papers{i}')
    # plt.plot(t_now, exponential_func(t_now, le_found[i], diff_now[0], 0), '-', label=f'papers{i}')
    plt.title(labels[i])
    plt.xlabel('t (s)')
    plt.ylabel('diff (deg)')
    plt.yscale('log') # set y-axis to semi-log scale
    plt.legend()
plt.show()


