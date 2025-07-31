#' @export
calculate_statistics <- function(data_column, save = FALSE, plots = TRUE) {
  # Load required packages
  library(tidyverse)
  library(data.table)
  library(scales)
  library(nortest)
  library(e1071)
  library(fastqq)
  
  # Start time of analysis
  start_time <- Sys.time()
  
  # Check data type
  column_type <- class(data_column)
  
  # Extract the column name
  column_name <- deparse(substitute(data_column))
  
  # Initialize the progress bar
  pb <- txtProgressBar(min = 0, max = 3, style = 3)
  
  # Helper functions
  calculate_numeric_statistics <- function(data_column, column_name, save, plots, pb) {
    
    # Total number of data points
    total_data_points <- length(data_column)
    
    # Remove NA values
    na_count <- sum(is.na(data_column))
    data_column <- na.omit(data_column)
    non_na_data_points <- total_data_points - na_count
    
    # Percentage of NA values
    percent_na <- (na_count / total_data_points) * 100
    
    # Generate a sample of the data for normality tests + Q-Q plot
    #set.seed(42)
    data_sample <- sample(data_column, 
                          min(300, length(data_column))
                          )
    
    # Calculations
    Q1 <- quantile(data_column, 0.25)
    Q3 <- quantile(data_column, 0.75)
    IQR <- Q3 - Q1
    
    min_IQR <- Q1 - 1.5 * IQR
    max_IQR <- Q3 + 1.5 * IQR
    
    median_value <- median(data_column)
    sd_value <- sd(data_column)
    
    min_sd <- median_value - 3 * sd_value
    max_sd <- median_value + 3 * sd_value
    
    outliers_IQR <- sum(data_column < min_IQR | data_column > max_IQR)
    outliers_sd <- sum(data_column < min_sd | data_column > max_sd)
    
    percent_outliers_IQR <- (outliers_IQR / non_na_data_points) * 100
    percent_outliers_sd <- (outliers_sd / non_na_data_points) * 100
    
    skewness_value <- skewness(data_column)
    kurtosis_value <- kurtosis(data_column)
    excess_kurtosis_value <- kurtosis_value - 3
    
    # Normality tests if small data size
    if (non_na_data_points <= 5000) {
      ks_test_result <- lillie.test(data_sample)
      pearson_test_result <- pearson.test(data_sample)
      #only run shapiro wilks test if non-na-data-points =< 5000, else do nothing
      shapiro_test_result <- shapiro.test(data_column)
      }
    
    
    # Print results
    cat("\n\n")
    cat("## Summary Statistics for", column_name, ":\n\n")
    cat("Number of NA values:", na_count, "\n")
    cat("Percentage of NA values:", percent_na, "%\n\n")
    print(summary(data_column))
    cat("\n")
    cat("Standard Deviation (na.rm = TRUE):", sd(data_column, na.rm = TRUE), "\n\n")
    
    # Print test results
    cat("\n\n")
    cat("## Normality info for", column_name, ": \n")
    cat("Skewness (normal distribution between -1 and 1):", skewness_value, "\n")
    cat("Excess kurtosis (normal distribution between -1 and 1):", excess_kurtosis_value, "\n")
    #Print shapiro-wilks test result if non-na-data-points =< 5000, else print "Shapiro-Wilks test not performed, sample size too large"
    if (non_na_data_points <= 5000) {
      cat("Shapiro-Wilks test: W =", shapiro_test_result$statistic, ", p-value =", shapiro_test_result$p.value, "\n")
      print(ks_test_result)
      print(pearson_test_result)
    } else {
      cat("Most normality tests (Shapiro-Wilks, Kolmogorov–Smirnov) not performed due to sample size too large\n")
    }

    cat("\n\n\n")

    
    cat("## Outlier information: ## \n\n")
    cat("min_IQR:", min_IQR, "\n")
    cat("max_IQR:", max_IQR, "\n")
    cat("Number of outliers (IQR):", outliers_IQR, "\n")
    cat("Percentage of outliers (IQR):", percent_outliers_IQR, "%\n\n")
    cat("min_sd:", min_sd, "\n")
    cat("max_sd:", max_sd, "\n")
    cat("Number of outliers (SD):", outliers_sd, "\n")
    cat("Percentage of outliers (SD):", percent_outliers_sd, "%\n")
    
    if (save) {
      # Create a subfolder in the working directory if it doesn't exist
      if (!file.exists("data-analysis")) {
        dir.create("data-analysis")
      }
      # Save results to a text file in the subfolder
      save_to_file <- function(...) {
        capture.output(..., file = file.path("data-analysis", paste0(column_name, "_numeric_statistics", ".txt")))
      }
      save_to_file(
        cat("\n\n Summary Statistics for", column_name, ":\n\n"),
        cat("Number of NA values:", na_count, "\n"),
        cat("Percentage of NA values:", percent_na, "%\n\n"),
        capture.output(print(summary(data_column)))
      )
    }
    
    # Generate and save plots if plots is TRUE
    if (plots) {
      ## Boxplot
      boxplot(data_column, main = paste("Boxplot of", column_name), ylab = "Values", col = "lightblue")
      
      ### Add lines for min_IQR and max_IQR
      abline(h = min_IQR, col = "red", lty = 2)
      abline(h = max_IQR, col = "red", lty = 2)
      text(x = 1.1, y = min_IQR, labels = "min_IQR", col = "red", pos = 4)
      text(x = 1.1, y = max_IQR, labels = "max_IQR", col = "red", pos = 4)
      
      ### Add lines for min_sd and max_sd
      abline(h = min_sd, col = "blue", lty = 3)
      abline(h = max_sd, col = "blue", lty = 3)
      text(x = 1.2, y = min_sd, labels = "min_sd", col = "blue", pos = 4)
      text(x = 1.2, y = max_sd, labels = "max_sd", col = "blue", pos = 4)
      
      # Generate a Q-Q plot with a sample of the data
      fastqq::qqnorm(data_column, main = paste("Q-Q Plot of", column_name))
      qqline(data_column, col = "red")
      
      if (save) {
        # Create a subfolder in the working directory if it doesn't exist
        if (!file.exists("data-analysis")) {
          dir.create("data-analysis")
        }
        # Save plots as PNG files in the data-analysis folder
        png(filename = file.path("data-analysis", paste0(column_name, "_Boxplot",  ".png")))
        boxplot(data_column, main = paste("Boxplot of", column_name), ylab = "Values", col = "lightblue")
        abline(h = min_IQR, col = "red", lty = 2)
        abline(h = max_IQR, col = "red", lty = 2)
        text(x = 1.2, y = min_IQR, labels = "min_IQR", col = "red", pos = 4)
        text(x = 1.2, y = max_IQR, labels = "max_IQR", col = "red", pos = 4)
        abline(h = min_sd, col = "blue", lty = 3)
        abline(h = max_sd, col = "blue", lty = 3)
        text(x = 1.2, y = min_sd, labels = "min_sd", col = "blue", pos = 4)
        text(x = 1.2, y = max_sd, labels = "max_sd", col = "blue", pos = 4)
        dev.off()
        
        png(filename = file.path("data-analysis", paste0(column_name,"_QQplot",  ".png")))
        fastqq::qqnorm(data_column, main = paste("Q-Q Plot of", column_name))
        qqline(data_column, col = "red")
        dev.off()
        
        cat("\nSaved Boxplot and Q-Q Plot as PNG files in the 'data-analysis' folder.\n")
      }
    }
    
    setTxtProgressBar(pb, 3)
  }
  
  calculate_character_statistics <- function(data_column, column_name, save = FALSE, plots = TRUE, pb = NULL) {
    # Total number of data points
    total_data_points <- length(data_column)
    
    # Remove NA values
    na_count <- sum(is.na(data_column))
    data_column <- na.omit(data_column)
    non_na_data_points <- total_data_points - na_count
    
    # Percentage of NA values
    percent_na <- (na_count / total_data_points) * 100
    
    # Frequency table
    frequency_table <- table(data_column)
    
    # Sort frequency table by frequency and reverse the order for descending order
    frequency_table <- frequency_table[order(-frequency_table)]
    
    # Limit the frequency table to the top 500 entries if the number of unique entries is greater than 200
    if (length(unique(data_column)) > 500) {
      frequency_table <- head(frequency_table, 500)
    }
    
    # Calculate percentage frequency
    percent_frequency <- prop.table(frequency_table) * 100
    
    # Create a new variable for legend combining category name and count
    legend_labels <- paste(names(frequency_table), " (", format(as.numeric(frequency_table), big.mark = ",", scientific = FALSE), ")", sep = "")
    
    # Remove extra spaces in legend labels
    legend_labels <- gsub("\\(\\s+", "(", legend_labels)
    
    # Manually reorder levels of factor for correct legend order
    frequency_levels <- factor(names(frequency_table), levels = names(frequency_table)[order(-frequency_table)])
    legend_labels <- factor(legend_labels, levels = legend_labels[order(-frequency_table)])
    
    # Print results
    cat("\n\n Summary Statistics for", column_name, ":\n\n")
    cat("Number of NA values:", na_count, "\n")
    cat("Percentage of NA values:", percent_na, "%\n\n")
    cat("\n\nFrequency Table:\n\n")
    print(head(frequency_table, 200))
    
    if (save) {
      # Create a subfolder in the working directory if it doesn't exist
      if (!file.exists("data-analysis")) {
        dir.create("data-analysis")
      }
      # Save results to a text file in the subfolder
      save_to_file <- function(...) {
        capture.output(..., file = file.path("data-analysis", paste0(column_name, "_character_statistics",  ".txt")))
      }
      save_to_file(
        cat("\n\n Summary Statistics for", column_name, ":\n\n"),
        cat("Number of NA values:", na_count, "\n"),
        cat("Percentage of NA values:", percent_na, "%\n\n"),
        cat("\n\nFrequency Table:\n\n"),
        capture.output(print(head(frequency_table, 200)))
      )
    }
    
    # Generate and save plots if plots is TRUE
    if (plots) {
      ## Bar chart
      bar_chart <- ggplot(data = data.frame(Category = frequency_levels, Frequency = as.numeric(frequency_table)), aes(x = Category, y = Frequency, fill = legend_labels)) +
        geom_bar(stat = "identity") +
        labs(title = paste("Bar Chart of", column_name), x = column_name, y = "Frequency") +
        scale_y_continuous(labels = comma) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              panel.grid = element_blank()) +
        if (length(unique(data_column)) <= 15) {
          guides(fill = guide_legend(title = "Frequency Table"))
        } else {
          guides(fill = "none")
        }
      
      print(bar_chart)
      
      if (save) {
        # Create a subfolder in the working directory if it doesn't exist
        if (!file.exists("data-analysis")) {
          dir.create("data-analysis")
        }
        # Save plot as a PNG file in the data-analysis folder
        ggsave(filename = file.path("data-analysis", paste0(column_name, "_BarChart", ".png")), plot = bar_chart)
        
        cat("\nSaved Bar Chart as a PNG file in the 'data-analysis' folder.\n")
      }
      
    }
    
    if (!is.null(pb)) {
      setTxtProgressBar(pb, 3)
    }
  }
  
  # Check if the column type is numeric
  if (column_type == "numeric") {
    # Run calculations for numeric type
    setTxtProgressBar(pb, 1)
    cat("\n\nProcessing numeric type column.\n")
    calculate_numeric_statistics(data_column, column_name, save, plots, pb)
  } else if (column_type == "character") {
    # Run calculations for character type
    setTxtProgressBar(pb, 1)
    cat("\n\nProcessing character type column.\n")
    calculate_character_statistics(data_column, column_name, save, plots, pb)
  } else {
    # Unsupported data type
    cat("\n\nUnsupported data type. Please provide either numeric or character data.")
    close(pb)
  }
  
  # Close the progress bar
  close(pb)
  
  # Print analysis done with time taken
  cat("\n\n Analysis done. Time taken:", Sys.time() - start_time, "seconds.\n")
} 


