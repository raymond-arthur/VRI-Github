Calculate Statistics (_calcstats_)
================
A fast and streamlined R function for comprehensive statistical analysis of numeric and character data columns, with focus on normality assessment for large datasets.

- [Features](#Features)
- [Installation](#Installation)
- [Usage](#Usage)
- [Documentation](#Documentation)
- [Code](#Code)
  - [Code aim and goals](#code-aims-and-goals:)
- [Usage and information](#Usage)


## Features

- Automatic handling of both numeric and character data
- Comprehensive statistical summaries
- Normality assessment (skewness, kurtosis, Shapiro-Wilk, Kolmogorov-Smirnov tests)
- Outlier detection using IQR and standard deviation methods
- Visualization tools (QQ plots, boxplots, bar charts)
- Progress tracking and time measurement
- Option to save results and plots to files


## Installation

```r
# Install required packages if not already installed
install.packages(c("tidyverse", "data.table", "scales", "nortest", "e1071", "fastqq"))
```

# Then source the function directly from GitHub:
```r
source("https://raw.githubusercontent.com/raymond-arthur/VRI-Github-stuff/main/calculate_statistics.R")
```

<br>
<br>

## Usage
# For numeric data
calculate_statistics(your_data$numeric_column, save = TRUE, plots = TRUE)

# For character data
calculate_statistics(your_data$character_column, save = FALSE, plots = TRUE)

<br>
<br>
Parameters

    data_column: The column to analyze (numeric or character)
    
    save: Logical, whether to save results to files (default: FALSE)
    
    plots: Logical, whether to generate plots (default: TRUE)

## Documentation

Documentation for _calcstats_ can be found here:
<https://github.com/raymond-arthur/VRI-Github-stuff/> or with the R `?help` function

```r
?help calcstats
```

<br> <br>

## Code:


The entire codebase for this package can be found here:
<https://github.com/raymond-arthur/VRI-Github-stuff/blob/main/calculate_statistics.R>
<br> <br>


### Code aims and goals:

The goal for this package is to create a self-contained and easily usable function that allows for the determination of normality of a large dataset.

We achieve this through a combination of the following:

1/ Heuristic determination of normality

We achieve this by visually examining the Q-Q plot of the data, determinating the quantiles and IQR, determinating of number of outliers, and visually inspection the data.

<br> <br>
An example of this output can be found here:

![](calculate_statistics.png)<!-- -->


