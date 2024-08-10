import numpy as np

# constants

# discarded metric unit values
m1 = 2 * 0.3048 * 0.0381 * 0.009525 * 2710  # mass of both arms in upper half
m2 = 0.2286 * 0.0381 * 0.009525 * 2710 # mass of only arm in lower half
l_upper = 0.3048 # same length of both arms in upper half
l_lower = 0.2286 # length of only arm in lower half
g = 9.80665 # m/s2 acceleration due to gravity

# alternate values for mass and length approximations (in inches)
# m1 = 2 * 12 * 1.5 * 0.375 * 0.04440894
# m2 = 9 * 1.5 * 0.375 * 0.04440894
# l_upper = 12
# l_lower = 9
# g = 386.088583 # in/s2 acceleration due to gravity

# alternate values for length approximations
l1 = l_upper
l2 = l_lower 

I1 = m1 * (l_upper ** 2) / 3
I2 = m2 * (l_lower ** 2) / 3

print(f"m1 = {m1} kg \nm2 = {m2} kg \nl1 = {l1} m \nl2 = {l2} m \nI1 = {I1} kg m^2 \nI2 = {I2} kg m^2")

# constants needed for calculation of normal modes
mu1 = ((m1 / 2) + m2) * g * l1
mu2 = (m2 * g * l2) / 2
Ia = ((m1 * (l1 ** 2)) / 4) + ((l1 ** 2) * m2) + I1
Ib = ((m2 * (l2 ** 2)) / 4) + I2
Iab = (m2 * l1 * l2) / 2

f1 = np.sqrt((-((Ia * mu2) + (Ib * mu2)) + np.sqrt((((Ia * mu2) + (Ib * mu2)) ** 2) + (4 * ((Iab * Iab) - (Ia * Ib)) * mu1 * mu2))) / (2 * ((Iab * Iab) - (Ia * Ib))))
f2 = np.sqrt((-((Ia * mu2) + (Ib * mu2)) - np.sqrt((((Ia * mu2) + (Ib * mu2)) ** 2) + (4 * ((Iab * Iab) - (Ia * Ib)) * mu1 * mu2))) / (2 * ((Iab * Iab) - (Ia * Ib))))

print(f"\nNormal mode frequencies are {f1:.3f} Hz and {f2:.3f} Hz.")
print(f"Normal mode time periods are {1 / f1:.3f} s and {1 / f2:.3f} s.")