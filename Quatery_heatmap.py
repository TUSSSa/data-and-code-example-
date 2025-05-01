import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the time series dataset
file_path = "CVA6_count.xlsx"
data = pd.read_excel(file_path, parse_dates=['Month'], index_col='Month')

# Select data from 2010 to 2022
data = data.loc['2005-01-01':'2023-12-31']

# Create a quarter column
data['Year'] = data.index.year
data['Quarter'] = data.index.quarter

# Aggregate data by year and quarter
pivot_data = data.pivot_table(values='Count', index='Year', columns='Quarter', aggfunc='sum', fill_value=0)

# Customize the maximum value for the legend
vmax_value = 1702  # Adjust this value as needed

# Set font sizes
plt.rcParams.update({
    'font.size': 10,  # Base font size
    'axes.labelsize': 16,  # Axis label font size
    'xtick.labelsize': 18,  # X-axis tick label font size
    'ytick.labelsize': 18  # Y-axis tick label font size
})

# Plot the heatmap
plt.figure(figsize=(10, 8))  # Adjust the figure width
sns.heatmap(pivot_data, annot=False, fmt='.0f', cmap='YlGnBu')
plt.title('CVA6', fontsize=40)

# Remove axis labels
plt.xlabel('')
plt.ylabel('')

# Adjust subplot layout
plt.subplots_adjust(top=0.88, bottom=0.11, left=0.22, right=0.9, hspace=0.2, wspace=0.2)

# Save the image with specified resolution
plt.savefig('heatmap.png', dpi=300)  # Set dpi to 300, adjust as needed

plt.show()

