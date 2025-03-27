require(data.table)
library(ggplot2)

#outdir <- "dorado_MD-NOVA_5polyA_OutSorted/v2/"


#bam.all <- fread("dorado_MD-NOVA_OutSorted/v2/bam.all.tsv")


cols_to_keep <- c("seqnames", "aln_ID", "mapq", "strand", "start", "end", "tag.l3", "tag.l5", "tag.r3", "tag.r5")  

bam.all.tags <- bam.all[, ..cols_to_keep]


cols_to_clean <- c("tag.l3", "tag.l5", "tag.r3", "tag.r5")

#remove the numbers and commas
bam.all.tags[, (cols_to_clean) := lapply(.SD, function(col) {
  # Remove numbers and trailing commas
  col <- gsub("\\b\\d+\\.?\\d*\\b,?", "", col)
  # Remove leading/trailing commas and spaces
  col <- gsub("^,|,$", "", col)
  trimws(col)
}), .SDcols = cols_to_clean]


### check the potential values in tag columns
cols_to_check <- c("tag.l3", "tag.l5", "tag.r3", "tag.r5")

# Print unique values for each column
for (col in cols_to_check) {
  cat("\nUnique values in", col, ":\n")
  print(unique(bam.all.tags[[col]]))
}


# Combine the four columns into a single key for counting combinations

    #bam.all.tags[, combined := paste(tag.l3, tag.l5, tag.r3, tag.r5, sep = ", ")]

# Count the occurrences of each combination

    #combination_counts <- table(bam.all.tags$combined)
combination_counts <- bam.all.tags[,
                         .N,
                         .(tag.l3, tag.l5, tag.r3, tag.r5)
                         ][order(-N)]

# Combine the four columns into a single key for counting combinations
combination_counts[, combined := paste(tag.l3, tag.l5, tag.r3, tag.r5, sep = ", ")]


# Convert to data.table for better handling and readability
#combination_counts_dt <- as.data.table(as.data.frame(combination_counts))
setnames(combination_counts, "N", "Count")

output_file <- file.path(outdir, "combination_counts.tsv")
fwrite(combination_counts, output_file, sep = "\t")


combination_counts <- fread(output_file, sep = "\t")

### Create barplot
# number of top bars to display

num_bars <- 20 

# Sort the combinations by count (optional)
combination_counts <- combination_counts[order(-Count)]

# Limit the number of bars to the top 'num_bars' combinations
combination_counts   <- head(combination_counts, num_bars)

combination_counts_m <- melt.data.table(combination_counts, id.vars = c('Count', 'combined'), variable.name = 'end', value.name = 'tag')


# Create a bar plot of the combination counts and store it in a variable
plot <- ggplot(combination_counts, aes(x = reorder(combined, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(title = paste("Top", num_bars, "Most Abundant Combinations Across Columns"),
       x = "Combination of tag.l3, tag.l5, tag.r3, tag.r5",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 


plot

# Specify the output file path for the plot
output_file_plot <- file.path(outdir, "combination_counts_plot.png")

# Save the plot to a PNG file
ggsave(output_file_plot, plot = plot, width = 10, height = 6)
