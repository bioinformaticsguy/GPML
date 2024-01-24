import pandas as pd

# Sample DataFrame
data = {'column1': [1, 2, 3, 4, 5],
        'column2': ['A', 'B', 'A', 'C', 'A']}
df = pd.DataFrame(data)

# Display the original DataFrame
print("Original DataFrame:")
print(df)

# Remove rows where 'column2' has the value 'A'
df = df[df['column2'] != 'A']

# Display the modified DataFrame
print("\nDataFrame after removing rows with 'A' in 'column2':")
print(df)
