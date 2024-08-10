import numpy as np


# Function to calculate the mean and standard deviation
def calculate_mean_std(data):
    # Extract values and errors from the tuples
    values = [x[0] for x in data]
    errors = [x[1] for x in data]

    # Calculate mean and standard deviation
    mean = np.mean(values)
    std = np.std(values)

    return mean, std

# Function to calculate the weighted average and its error
def calculate_weighted_average(data):
    # Extract values and errors from the tuples
    values = [x[0] for x in data]
    errors = [x[1] for x in data]

    # Calculate weights as inverse of squared errors
    weights = [1 / (error ** 2) for error in errors]

    # Calculate weighted average
    weighted_average = np.average(values, weights=weights)

    # Calculate error of weighted average
    weighted_average_error = np.sqrt(1 / np.sum(weights))

    return weighted_average, weighted_average_error

# Example usage [(value, error), (value, error), ...]
data = [(0.9106,0),(0.8985,0),(0.8985,0)] #NM1
data = [(1.7970,0),(1.8090,0),(1.7970,0)] #NM2
data = [(20.678,1.141),(20.485,0.585),(21.693,1.926),(20.371,1.039),(18.409,2.588),(20.056,0.653)]
data = [(14.498,2.36),(14.365,3.337),(16.943,1.301),(14.785,1.102),(15.559,1.73),(15.729,2.02)]
data = [(18.152,3.111),(17.727,3.086),(17.862,2.499),(17.444,3.347),(18.747,9.298),(18.787,5.941),(18.370,3.895),(18.112,3.633)]

mean, std = calculate_mean_std(data)
print("Mean:", mean)
print("Standard Deviation:", std)

weighted_average, weighted_average_error = calculate_weighted_average(data)
print("Weighted Average:", weighted_average)
print("Weighted Average Error:", weighted_average_error)
