import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


# Define the exponential function to fit
def exponential_func(x, le, d0 ):
    return d0 * np.exp(le * x) 

def fit_exponential(t, diff, d0, le0):
    # Fit the exponential function to the data
    # popt, pcov = curve_fit(exponential_func, t, diff)
    # popt, pcov = curve_fit(exponential_func, t, diff, bounds=([-np.inf, d0], [np.inf, d0 + 1e-7]))
    popt, pcov = curve_fit(exponential_func, t, diff, p0=[le0, d0])

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
# d0s = [diff[0] for diff in diffs]

for i, diff in enumerate(diffs):
    if i != 1:
        continue
    start = 0
    t_now = np.copy(t[start:len_found[i]+start]-t[0])
    diff = diff[start:len_found[i]+start]
    popt, perr, chi_squared, red_chi_squared = fit_exponential(t_now, diff, diff[0], le_found[i])

    # print the results
    print("For diff{}:".format(i+1))
    # print("Fitted parameters: le = {:.3f} +/- {:.3f}, d0 = {:.3f} +/- {:.3f}, c = {:.3f} +/- {:.3f}".format(popt[0], perr[0], popt[1], perr[1], popt[2], perr[2]))
    print("Fitted parameters: le = {:.3f} +/- {:.3f}, d0 = {:.3f} +/- {:.3f}".format(popt[0], perr[0], popt[1], perr[1]))
    print("Chi-squared value: {:.3f}, reduced ch-squared value: {:.3f}".format(chi_squared,red_chi_squared))
    print()

    # plot the best fit
    plt.plot(t[:len(diffs[i])], diffs[i], 'o', label=f'data{i}')
    plt.plot(t_now, exponential_func(t_now, *popt), '-', label=f'fit{i}')
    plt.plot(t_now, exponential_func(t_now, le_found[i], diff[0]), '-', label=f'papers{i}')
    # plt.plot(t_now, exponential_func(t_now, le_found[i], diff[i], 0), '-', label=f'papers{i}')
    plt.xlabel('t (s)')
    plt.ylabel('diff (deg)')
    # plt.yscale('log') # set y-axis to semi-log scale
    plt.yscale('log') # set y-axis to semi-log scale
    plt.legend()
plt.show()
