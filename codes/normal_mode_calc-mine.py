import numpy as np

# constants

# metric unit values
m1 = 2 * 0.3048 * 0.0381 * 0.009525 * 2710  # mass of both arms in upper half
m2 = 0.2286 * 0.0381 * 0.009525 * 2710 # mass of only arm in lower half
l_upper = 0.3048 # same length of both arms in upper half
l_lower = 0.2286 # length of only arm in lower half
g = 9.80665 # m/s2 acceleration due to gravity
b_upper = 0.0381 # breadth of both arms in upper half
b_lower = 0.0381 # breadth of only arm in lower half
d_upper = 0.009525 # depth of both arms in upper half
d_lower = 0.009525 # depth of only arm in lower half


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

l2 = l_lower # length of arm till center of mass
l1 = l_upper # length of arms till center of mass
b1 = b_upper # breadth of arms in upper half
b2 = b_lower # breadth of arm in lower half
d1 = d_upper # depth of arms in upper half
d2 = d_lower # depth of arm in lower half


I1 = (1/12) * m1 * (l1 ** 2) # moment of inertia of both arms in upper half
I2 = (1/12) * m2 * (l2 ** 2) # moment of inertia of only arm in lower half

Ja = ((1/4) * m1 * (l1 ** 2)) + (m2 * (l1 ** 2) + I1)
Jb = ((1/4) * m2 * (l2 ** 2) + I2)
Jx = ((1/2) * m2 * l1 * l2)
mu1 = (((1/2) * m1) + m2) * g * l1
mu2 = (1/2) * m2 * g * l2

print(f"m1 = {m1} kg \nm2 = {m2} kg \nl1 = {l1} m \nl2 = {l2} m \n m \nJa = {Ja} kg m^2 \nJb = {Jb} kg m^2 \nJx = {Jx} kg m^2 \nmu1 = {mu1} kg m \nmu2 = {mu2} kg m")
print()


a = -1 * ((Jx ** 2) - (Ja * Jb))
b = (mu1 * Jb) + (mu2 * Ja)
c = mu1 * mu2

b = -b # to match with the equation in the paper

print(f"a = {a} \nb = {b} \nc = {c}")
print(f"b^2 - 4ac = {(b ** 2) - (4 * a * c)}")
print(f"sqrt(b^2 - 4ac) = {np.sqrt((b ** 2) - (4 * a * c))}")
print(f"2a = {2 * a}")

f1 = np.sqrt(np.abs((-b + np.sqrt((b ** 2) - (4 * a * c))) / (2 * a)))
f2 = np.sqrt(np.abs((-b - np.sqrt((b ** 2) - (4 * a * c))) / (2 * a)))

print(f"\nNormal mode frequencies are {f1:.3f} Hz and {f2:.3f} Hz.")
print(f"Normal mode time periods are {1 / f1:.3f} s and {1 / f2:.3f} s.")