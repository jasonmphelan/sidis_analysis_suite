import numpy as np
import matplotlib.pyplot as plt

# Parameters for the Poisson-distributed yields
lambda_x = 3  # Expected yield for x
lambda_y = 4  # Expected yield for y
n_samples = 1000000  # Number of Monte Carlo samples

# Generate samples of x and y from Poisson distributions
x_samples = np.random.poisson(lambda_x, n_samples)
y_samples = np.random.poisson(lambda_y, n_samples)

# Calculate the yield ratio for each sample
z_samples = x_samples / (x_samples + y_samples)
z_samples = z_samples[~np.isnan(z_samples)]

lambda_z = lambda_x/(lambda_x+lambda_y)

# Calculate the mean and standard deviation of the ratio
z_mean = np.mean(z_samples)
z_std = np.std(z_samples)

prop = 1/(lambda_x + lambda_y)**2 * np.sqrt(lambda_x*lambda_x*lambda_y + lambda_y*lambda_y*lambda_x)

# Display the results
print(f"Mean of the yield ratio: {z_mean:.4f}")
print(f"Standard deviation of the yield ratio: {z_std:.4f}")
print(f"Std from error propagation: {prop:.4f}")

print(f"Sampling: {z_std/lambda_z:.4f}")
print(f"Prop: {prop/lambda_z:.4f}")

# Plot the distribution of the yield ratio
plt.hist(z_samples, bins=50, density=True, alpha=0.7, color='g')
plt.title(f"Distribution of Yield Ratio (x / (x + y))\nMean = {z_mean:.4f}, Std = {z_std:.4f}")
plt.xlabel('Yield Ratio')
plt.ylabel('Density')
plt.show()

