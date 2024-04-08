import pandas as pd
import matplotlib.pyplot as plt



# Example usage:
data = {
    'Category': ['A', 'B', 'A', 'C', 'B', 'A', 'A', 'B', 'C', 'C'],
    'Value': [10, 20, 15, 25, 30, 12, 18, 22, 28, 35]
}
df = pd.DataFrame(data)

plot_pie_with_counts(df, 'Category')
