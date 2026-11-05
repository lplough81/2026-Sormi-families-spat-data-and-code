# Processing code from Ariana Huffmyer's github page 

# Load libraries. 
```{r}
library(MASS) 
library(tidyverse)
library(ggplot2)
library(readxl)
library(cowplot)
library(lme4)
library(lmerTest)
library(car)
library(effects)
library(emmeans)
library(tools)
library(pracma)
```

# Resazurin metabolic rates

## Load data 

```{r}
# Set the folder path
folder_path <- "path/to/files"  

# List all txt files in the folder
file_list <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE, recursive=TRUE)

# Check if any files are found
if (length(file_list) == 0) {
  stop("No txt files found in the specified folder. Check the folder path or file extension.")
}

# Initialize an empty list to store processed data
data_list <- list()

# Loop through each file and load the data
for (file in file_list) {
  # Ensure file is a character string
  if (!is.character(file)) next
  
  # Extract file name without extension
  file_name <- file_path_sans_ext(basename(file))
  
  # Read the text file, skipping the first row (metadata)
  data <- read.delim(file, header = FALSE, skip = 2)
  
  # Assign column names: first column is row labels (A–H), others are 1–12
  colnames(data) <- c("Row", as.character(1:12))
  
  #remove last column
  data<-data%>%dplyr::select(-ncol(data))
  
  # Convert to long format
  data_long <- data %>%
    pivot_longer(cols = -Row, names_to = "Column", values_to = "Value") %>%
    mutate(
      Column = sprintf("%02d", as.integer(Column)),  # Ensure two-digit column numbers
      Well_ID = paste0(Row, Column),                # Format well ID as "A01", "B02", etc.
      FileName = file_name,                          # Add file name
      plate = str_extract(file_name, "plate\\d+"),  # Extract "plateX"
      date = str_extract(file_name, "^\\d{8}"),     # Extract 8-digit date
      timepoint = str_extract(file_name, "T\\d+") %>% 
        str_remove("T") %>% 
        as.numeric()                                # Convert timepoint to numeric
    ) %>%
    dplyr::select(FileName, Well_ID, Value, date, plate, timepoint)  # Select relevant columns
  
  # Store the processed data in the list
  data_list[[file_name]] <- data_long
}

# Print an example of processed data
head(data_list[[file_name]])


# Combine all data frames into a single data frame (optional)
combined_data <- bind_rows(data_list, .id = "Source")

# Print the first few rows of the combined data (optional)
head(combined_data)

# Rename columns
combined_data<-combined_data%>%
  rename("well"=Well_ID, resazurin_counts=`Value`)%>%
  mutate(timepoint=as.character(timepoint))

head(combined_data)
```

Load in metadata. 

```{r}
metadata<-read_xlsx(path="/Users/louis.plough/Family Resazurin Trials/Fam42 Rez/MetaData42.xlsx")%>%
  mutate(date=as.character(date))
```

Join with data frame. 

```{r}
str(combined_data)
str(metadata)

full_data<-left_join(combined_data, metadata, by=c("date", "plate", "well"))%>%
  filter(!is.na(type))

head(full_data)
```

Load in size data. 

```{r}
size<-read_xlsx("/Users/louis.plough/Family Resazurin Trials/Fam42 Rez/Size_42.xlsx")%>%
  mutate(date=as.character(date))
```

Join with data. 
####Size_42.xlsx needs to have col name Plate changed to plate...
```{r}
full_data<-left_join(full_data, size, by=c("date", "plate", "well"))

###full_data<-full_data%>%
  mutate(date=as.character(date))
```

```{r}
head(full_data)
```

full_data%>%
  ggplot(aes(x=timepoint, y=resazurin_counts, colour=Family, group=interaction(date, plate, well)))+
  geom_point()+
  geom_line()+
  theme_classic()
```

Calculate fluorescence at each time point normalized to the starting value at time 0. 
```{r}
full_data<-full_data%>%
  group_by(date, plate, well, Family)%>%
  arrange(date, plate, well)%>%
  mutate(fluorescence.norm=resazurin_counts/first(resazurin_counts))
```
full_data<-full_data%>%
  mutate(value.mm=fluorescence.norm/area.mm)

###new correct for family...
family_curve <- full_data %>%
  group_by(timepoint, Family) %>%
  summarise(mean_fluor = mean(value.mm, na.rm = TRUE),
            se_fluor = sd(value.mm, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")




family_curve <- full_data %>%
  group_by(timepoint, Family) %>%
  summarise(mean_fluor = mean(fluorescence.norm, na.rm = TRUE),
            se_fluor = sd(fluorescence.norm, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# Plot mean curves with SE ribbon
ggplot(family_curve, aes(x = timepoint, y = mean_fluor, color = Family, group = Family)) +
  geom_ribbon(aes(ymin = mean_fluor - se_fluor, ymax = mean_fluor + se_fluor,
                   fill = Family), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  labs(x = "Timepoint", y = "Normalized fluorescence (mean ± SE)",
       title = "Mean fluorescence curves by Family") +
  theme_minimal()

# AUC DATA....
auc_data<-full_data%>%
  mutate(time=as.numeric(timepoint))%>%
  group_by(date, plate, well, temperature) %>%
  summarise(AUC = trapz(time, value.mm))

head(auc_data)

aucs <- full_data %>%
  mutate(time=as.numeric(timepoint))%>%
  group_by(date, plate, well, Family) %>%
  summarise(AUC = trapz(time, value.mm))

plot10a<-aucs%>%
  group_by(Family)%>%
  summarise(mean=mean(AUC, na.rm=TRUE), se=(sd(AUC, na.rm=TRUE))/sqrt(length(AUC)))%>%
  
plot10a<-aucs%>%
  group_by(Family)%>%
  summarise(mean=mean(AUC, na.rm=TRUE), se=(sd(AUC, na.rm=TRUE))/sqrt(length(AUC)))%>%
  ggplot(aes(x = reorder(Family, mean), y = mean)) +
  geom_point()+
  scale_fill_manual(values=c("darkorange", "purple2"), name="")+
  scale_colour_manual(values=c("darkorange", "purple2"), name="")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se))+
  labs(title = "C", x = "Family", y = "Metabolism (AUC)") +
  theme_classic();plot10a

library(dplyr)
library(ggplot2)

# Summary stats
auc_summary <- aucs %>%
  group_by(Family) %>%
  summarise(
    n    = sum(!is.na(AUC)),
    mean = mean(AUC, na.rm = TRUE),
    se   = sd(AUC, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# Order families by mean, applied to both datasets - add jittered points.

fam_order <- auc_summary$Family[order(auc_summary$mean)]
auc_summary <- auc_summary %>% mutate(Family = factor(Family, levels = fam_order))
aucs_plot   <- aucs %>% filter(!is.na(AUC)) %>%
  mutate(Family = factor(Family, levels = fam_order))

plot10a <- ggplot() +
  geom_jitter(data = aucs_plot, aes(x = Family, y = AUC),
              width = 0.15, height = 0, alpha = 0.4, colour = "grey40") +
  geom_errorbar(data = auc_summary,
                aes(x = Family, ymin = mean - se, ymax = mean + se),
                width = 0.2) +
  geom_point(data = auc_summary, aes(x = Family, y = mean), size = 3) +
  labs(x = "Day set", y = "Metabolism (AUC)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot10a
