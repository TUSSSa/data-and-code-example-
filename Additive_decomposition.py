import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.seasonal import seasonal_decompose

# Load time series dataset
file_path = "CVA6_count.xlsx"
data = pd.read_excel(file_path, parse_dates=['Month'], index_col='Month')

# Select data from 2010 to 2019
data = data.loc['2010-01-01':'2019-12-31']

# Additive decomposition
result_add = seasonal_decompose(data['Count'], model='additive', period=12)

# Custom plotting
fig, axes = plt.subplots(4, 1, figsize=(16, 20), sharex=True)

# Original data
axes[0].plot(data['Count'], label='Original Count', linewidth=3, color='red')  # Set line width and color
axes[0].set_ylabel('Count', fontsize=30)  # Set y-axis font size
axes[0].tick_params(axis='both', which='major', labelsize=24)  # Set tick label font size

# Trend
axes[1].plot(result_add.trend, label='Trend', linewidth=3, color='green')  # Set line width and color
axes[1].set_ylabel('Trend', fontsize=30)  # Set y-axis font size
axes[1].tick_params(axis='both', which='major', labelsize=24)  # Set tick label font size

# Seasonality
axes[2].plot(result_add.seasonal, label='Seasonal', linewidth=3, color='blue')  # Set line width and color
axes[2].set_ylabel('Seasonal', fontsize=30)  # Set y-axis font size
axes[2].tick_params(axis='both', which='major', labelsize=24)  # Set tick label font size

# Residual
axes[3].plot(result_add.resid, label='Residual', linewidth=3, color='black')  # Set line width and color
axes[3].set_ylabel('Residual', fontsize=30)  # Set y-axis font size
axes[3].tick_params(axis='both', which='major', labelsize=24)  # Set tick label font size

plt.suptitle('CVA6', fontsize=50)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()