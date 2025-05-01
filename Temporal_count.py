import pandas as pd

# Read the Excel file
input_file = 'CVA6.xlsx'
df = pd.read_excel(input_file)

# Parse the date column
df['Date'] = pd.to_datetime(df['Month'])

# Count occurrences for each month of each year
df['Month'] = df['Date'].dt.to_period('M')
result = df.groupby('Month').size().reset_index(name='Count')

# Output the result to a new Excel file
output_file = 'CVA6_count.xlsx'
result.to_excel(output_file, index=False)

print(f"Saved to {output_file}")