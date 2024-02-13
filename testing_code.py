import pandas as pd
import matplotlib.pyplot as plt

# Create a dummy DataFrame
data = {
    'Category': ['A', 'B', 'A', 'C', 'B', 'A', 'A', 'B', 'C', 'C'],
    'Value': [10, 20, 15, 25, 30, 12, 18, 22, 28, 35]
}
df = pd.DataFrame(data)

# Count the frequency of each unique entry in the 'Category' column
value_counts = df['Category'].value_counts()

# Plot the pie chart
value_counts.plot.pie(autopct='%1.1f%%', startangle=90)
plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

# Add a title
plt.title('Pie Chart of Unique Entries in "Category" Column')

# Show the plot
plt.show()
