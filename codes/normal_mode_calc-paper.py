import numpy as np

# constants

# metric unit values
m1 = 2 * 0.3048 * 0.0381 * 0.009525 * 2710  # mass of both arms in upper half
m2 = 0.2286 * 0.0381 * 0.009525 * 2710 # mass of only arm in lower half
l_upper = 0.3048 # same length of both arms in upper half
l_lower = 0.2286 # length of only arm in lower half
b_upper = 0.0381 # breadth of both arms in upper half
b_lower = 0.0381 # breadth of only arm in lower half
g = 9.80665 # m/s2 acceleration due to gravity

# alternate values for mass and length approximations (in inches)
# m1 = 2 * 12 * 1.5 * 0.375 * 0.04440894
# m2 = 9 * 1.5 * 0.375 * 0.04440894
# l_upper = 12
# l_lower = 9
# g = 386.088583 # in/s2 acceleration due to gravity

# alternate values for mass and length approximations (in metres but physically scaled)
# m1 = 2 * 12 * 1.5 * 0.375 * 2710
# m2 = 9 * 1.5 * 0.375 * 2710
# l_upper = 12
# l_lower = 9
# g = 9.80665 # m/s2 acceleration due to gravity


l3 = l_upper # same length of both arms in upper half
l2 = l_lower / 2 # length of arm till center of mass
l1 = l_upper / 2 # length of arms till center of mass

# alternate values for length approximations
# l3 = l_upper 
# l2 = l_lower 
# l1 = l_upper

I1 = (1/3) * m1 * (l_upper ** 2) # moment of inertia of both arms in upper half
I2 = (1/3) * m2 * (l_lower ** 2) # moment of inertia of only arm in lower half

I1 = (m1 * ((l_upper ** 2) + (b_lower ** 2)) / 12)
I2 = (m2 * ((l_lower ** 2) + (b_lower ** 2)) / 12)

# I1 = (m1 * ((l_upper ** 2) + (b_lower ** 2)) / 12) + (m1 * (l_upper / 2) ** 2)
# I2 = (m2 * ((l_lower ** 2) + (b_lower ** 2)) / 12) + (m2 * (l_lower / 2) ** 2)

print(f"m1 = {m1} kg \nm2 = {m2} kg \nl1 = {l1} m \nl2 = {l2} m \nl3 = {l3} m \nI1 = {I1} kg m^2 \nI2 = {I2} kg m^2")

# constants needed for calculation of normal modes
c1 = ((m2 ** 2) * (l2 ** 2) * g * l3) - (I2 * g * ((m1 * l1) + (m2 * l3)))
c2 = ((m2 ** 2) * (l2 ** 2) * g * l3)
c3 = (I2 * g * ((m1 * l1) + (m2 * l3))) + (m2 * g * l2 * ((m1 * l1 * l3) - I1 - (m2 * l2 * l3)))
c4 = - m2 * g * l2 * ((m2 * l2 * l3) + I1 + (m2 * (l3 ** 2)))
d = (I2 * (I1 + (m2 * (l3 ** 2)))) - ((m2 * l2 * l3) ** 2)

# print(f"c1 = {c1} \nc2 = {c2} \nc3 = {c3} \nc4 = {c4} \nd = {d}")

f1 = np.sqrt((-(c1 + c4) + np.sqrt(((c1 + c4) ** 2) - (4 * ((c1 * c4) - (c2 * c3))))) / (2 * d))
f2 = np.sqrt((-(c1 + c4) - np.sqrt(((c1 + c4) ** 2) - (4 * ((c1 * c4) - (c2 * c3))))) / (2 * d))

print(f"\nNormal mode frequencies are {f1:.3f} Hz and {f2:.3f} Hz.")
print(f"Normal mode time periods are {1 / f1:.3f} s and {1 / f2:.3f} s.")