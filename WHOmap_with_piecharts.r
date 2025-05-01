# Set your working directory
setwd("C:\\Users\\ASUS\\Desktop\\ 代码\\Map_with_piecharts")

# Load necessary libraries
library(ggplot2)
library(scales)
library(readxl)    # Used to read Excel files
library(dplyr) 
library(lubridate) 
library(ggspatial)  # Ensure ggspatial package is installed
library(scatterpie)  # Ensure scatterpie package is installed and loaded
library(ggnewscale)  # Used to use multiple color scales in the same plot
library(sf)          # Used for handling spatial data
library(stringr)     # Used for string processing
library(ggimage)  # Add ggimage package for loading background images

# Read Excel file
data_path <- "yearly_proportions_with_continent(WHO).xlsx"  # Replace with your file path
data <- read_excel(data_path, sheet = 1)
data <- data %>% filter(Date != 0)

# Aggregate data by year and continent
data_summary <- data %>%
  group_by(Continent, Date) %>%
  summarise(across(everything(), ~ first(.)))  # Retain information from all columns

# Ensure Count column is numeric
data_summary$Count <- as.numeric(as.character(data_summary$Count))

# Read continent center coordinates
continent_centers_path <- "continent_centers (WHO).xlsx"
continent_centers <- read_excel(continent_centers_path)

# Merge data
data_summary <- data_summary %>%
  left_join(continent_centers, by = c("Continent" = "Continent"))

# Check merge result
print(head(data_summary))

# Get world map data
world_map <- map_data('world')

# Define file_columns
file_columns <- c("CVA2", "CVA4", "CVA6", "CVA10", "CVA16", "CVB3", "CVB5", "EV-D68", "EV-A71")  # Replace with actual column names

# Define time periods
time_periods <- list(
  "Before 2000" = data_summary %>% filter(Date < 2000),
  "Before 2005" = data_summary %>% filter(Date < 2006),
  "Before 2010" = data_summary %>% filter(Date < 2011),
  "Before 2015" = data_summary %>% filter(Date < 2016),
  "Before 2020" = data_summary %>% filter(Date <= 2020),
  "2001-2010" = data_summary %>% filter(Date <= 2010),
  "2011-2020" = data_summary %>% filter(Date <= 2020),
  "Before 2023" = data_summary %>% filter(Date <= 2023)
)

# Loop to generate and save heatmaps and pie charts
for (period_name in names(time_periods)) {
  # Extract current data
  current_data <- time_periods[[period_name]]
  
  # Aggregate data by continent
  aggregated_data <- current_data %>%
    group_by(Continent, Longitude.y, Latitude.y) %>%
    summarise(across(all_of(file_columns), sum, na.rm = TRUE), 
              Count = sum(Count, na.rm = TRUE)) %>%
    ungroup()
  
  # Check column names of aggregated_data
  colnames(aggregated_data)
  
  # Calculate percentages for each sample
  proportion_data <- aggregated_data %>%
    mutate(across(all_of(file_columns), ~ . / Count)) %>%
    mutate(Radius = pmax(log2(Count) / 3, 1))  # Modify radius calculation to log2(Count) / 3
  
  # Check column names of proportion_data
  colnames(proportion_data)
  
  # Plot gray map and pie charts
  combined_plot <- ggplot() +
    labs(title = paste0("Serotype   ", period_name)) +
    
    # Set theme: Remove panel grid and background, retain only white background
    theme(
      panel.background = element_rect(colour = NA, fill = 'transparent'),  # Set panel background color to transparent
      panel.grid.major = element_blank(),               # Remove major grid lines
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 36 ),
      legend.text = element_text(size = 14),            # Increase legend text size
      legend.title = element_text(size = 18),           # Increase legend title size
      legend.key.size = unit(1.4, "lines"),             # Increase legend key size
      panel.spacing = unit(0.1, "lines"),               # Reduce panel spacing
      plot.margin = margin(1, 1, 1, 1, unit = "lines"),  # Reduce margin space
      plot.background = element_rect(fill = 'transparent', color = NA)  # Set background to transparent
    ) +
    
    coord_fixed(ratio = 1) +
    
    # Plot pie charts
    geom_scatterpie(aes(x = Longitude.y, y = Latitude.y, r = Radius),  # Use radius
                    data = proportion_data, 
                    cols = file_columns, 
                    color = 'black', 
                    size = 0.3,
                    legend_name = "Serotype",
                    alpha = 1) +
    
    # Set pie chart theme
    theme(
      plot.background = element_rect(fill = 'transparent', color = NA),  # Set background to transparent
      panel.background = element_rect(fill = 'transparent', color = NA),  # Set panel background to transparent
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_blank(),  # Remove axis lines
      axis.text = element_blank(),  # Remove axis text
      axis.ticks = element_blank(),  # Remove axis ticks
      axis.title = element_blank()  # Remove axis titles
    )
  
  # Generate file name
  combined_file_name <- paste0("combined_plot1_", period_name, ".png")
  
  # Save combined plot
  tryCatch({
    ggsave(combined_file_name, plot = combined_plot, width = 3400, height = 2400, units = "px", device = "png")
    print(paste("Saved combined plot for", period_name))
  }, error = function(e) {
    print(paste("Error saving combined plot for", period_name, ":", e$message))
  })
}
#Use PS to put piecharts on WHOmap

