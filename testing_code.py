import numpy as np
import matplotlib.pyplot as plt

# Sample data
data = np.random.random((5, 5))

# Create a heatmap
plt.imshow(data, cmap='viridis', interpolation='nearest')

# Add a colorbar for reference
plt.colorbar()

# Set labels and title
plt.xlabel('X-axis Label')
plt.ylabel('Y-axis Label')
plt.title('Heatmap Example')

# Show the plot
plt.show()
