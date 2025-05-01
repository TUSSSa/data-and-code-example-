setwd("C:\\Users\\ASUS\\Desktop\\ 代码\\Map_with_piecharts")#set your workplace
library(ggplot2)
library(scales)
library(readxl)    # Used for reading Excel files
library(dplyr) 
library(lubridate) 
library(ggspatial)  # Ensure the ggspatial package is installed
library(scatterpie)  # Ensure the scatterpie package is installed and loaded
library(ggnewscale)  # For using multiple color scales in one plot
library(sf)          # For handling spatial data
library(stringr)     # For string processing

# Read Excel file
data_path <- "yearly_proportions.xlsx"  # Replace with your file path
data <- read_excel(data_path, sheet = 1)
data <- data %>% filter(Date != 0)

# Aggregate data by year and country
data_summary <- data %>%
  group_by(Nation, Date) %>%
  summarise(across(everything(), ~ first(.)))  # Retain information from all columns

# Standardize country names
data_summary <- data_summary %>% 
  mutate(Nation = case_when(
    Nation == "United Kingdom" ~ "United Kingdom",
    Nation == "Viet Nam" ~ "Vietnam",
    Nation == "Central African Republic" ~ "Central African Rep.",
    Nation == "Democratic Republic of the Congo" ~ "Dem. Rep. Congo",
    Nation == "Solomon Islands" ~ "Solomon Is.",
    Nation == "USA" ~ "United States of America",
    Nation == "UK" ~ "United Kingdom",
    TRUE ~ Nation  # Keep other values unchanged
  ))

# Ensure the Count column is numeric
data_summary$Count <- as.numeric(as.character(data_summary$Count))

# Read country center coordinates
country_centers_path <- "country_centers.xlsx"
country_centers <- read_excel(country_centers_path)

# Merge data
data_summary <- data_summary %>%
  left_join(country_centers, by = c("Nation" = "Nation"))

# Check merge results
print(head(data_summary))

# Get world map data
world_map <- map_data('world')

# Define file_columns
file_columns <- c("A2", "A4", "A6", "A10", "A16", "B3", "B5", "D68", "EV71")  # Replace with actual column names

# Create auxiliary data frame for legend
legend_data <- data.frame(
  Count = c(1, 10, 100, 1000),
  Radius = c(log10(1) + 1, log10(10) + 1, log10(100) + 1, log10(1000) + 1),
  Longitude = c(-100, -100, -100, -100),  # Choose a fixed position
  Latitude = c(50, 45, 40, 35)        # Choose a fixed position
)

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

# Loop to generate and save heat maps and pie charts
for (period_name in names(time_periods)) {
  # Extract current data
  current_data <- time_periods[[period_name]]
  
  # Aggregate data by country
  aggregated_data <- current_data %>%
    group_by(Nation, Longitude.y, Latitude.y) %>%
    summarise(across(all_of(file_columns), sum, na.rm = TRUE), 
              Count = sum(Count, na.rm = TRUE)) %>%
    ungroup()
  
  # Check column names of aggregated_data
  colnames(aggregated_data)
  
  # Calculate percentages for each sample
  proportion_data <- aggregated_data %>%
    mutate(across(all_of(file_columns), ~ . / Count)) %>%
    mutate(Radius = pmax(log10(Count) + 1, 1))  # Ensure radius is always positive and at least 1
  
  # Check column names of proportion_data
  colnames(proportion_data)
  
  # Merge data with map data
  merged_data <- left_join(world_map, aggregated_data, by = c("region" = "Nation"))
  
  # Plot gray map and pie charts
  combined_plot <- ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "gray90", color = 'black', size = 0.1) +
    labs(title = paste0("Enterovirus   ", period_name)) +
    
    # Set theme: Remove panel grids and background, retain only white background
    theme(
      panel.background = element_rect(colour = NA, fill = 'white'),  # Set panel background color
      panel.grid.major = element_blank(),               # Remove major grid lines
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 36 ),
      legend.text = element_text(size = 14),            # Increase legend text size
      legend.title = element_text(size = 18),           # Increase legend title size
      legend.key.size = unit(1.4, "lines"),             # Increase legend key size
      panel.spacing = unit(0.1, "lines"),               # Reduce panel spacing
      plot.margin = margin(1, 1, 1, 1, unit = "lines")  # Reduce margin space
    ) +
    
    coord_fixed(ratio = 1) +
    
    # Plot pie charts
    geom_scatterpie(aes(x = Longitude.y, y = Latitude.y, r = Radius),  # Use radius
                    data = proportion_data, 
                    cols = file_columns, 
                    color = 'black', 
                    size = 0.2,
                    legend_name = "Enterovirus",
                    alpha = 1) +
    
    # Plot legend
    geom_scatterpie_legend(legend_data$Radius, x = 180, y = 40, n = 4, labeller = function(x) 10^(x - 1)) +
    
    # Add country names
    #geom_text(data = proportion_data, aes(x = Longitude.y, y = Latitude.y, label = Nation), 
    #         color = "black", size = 3, vjust = -1) +
    
    # Set pie chart theme
    theme(
      plot.background = element_rect(fill = 'white', color = NA),  # Set background to white
      panel.background = element_rect(fill = 'white', color = NA),  # Set panel background to white
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
    ggsave(combined_file_name, plot = combined_plot, width = 7680, height = 4320, units = "px", device = "png")
    print(paste("Saved combined plot for", period_name))
  }, error = function(e) {
    print(paste("Error saving combined plot for", period_name, ":", e$message))
  })
}
