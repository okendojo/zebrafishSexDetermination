# Example script to generate a repeat annotation summary
grep -v "^#" genome.fasta.out | cut -f 11 | sort | uniq -c > repeat_summary.txt


# Example R script to create a bar plot
library(ggplot2)

# Load repeat summary data
repeat_summary <- read.table("repeat_summary.txt", header = FALSE, col.names = c("Count", "Type"))

# Create a bar plot
ggplot(repeat_summary, aes(x = Type, y = Count, fill = Type)) +
  geom_bar(stat = "identity") +
  labs(title = "Repeat Element Distribution", x = "Repeat Type", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot as an image file (e.g., PNG)
ggsave("repeat_distribution.png", width = 8, height = 6)

