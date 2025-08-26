# Copyright 2025 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


#!/usr/bin/env Rscript

suppressMessages(library(qqman))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11) {
  cat("Usage: Rscript manhattan.R <input_file> <score_column> <top_fraction1> <top_fraction2> <use_absolute> <width> <height> <title> <output_table1> <output_table2> <output_plot>\n")
  quit(status = 1)
}

input_file <- args[1]
score_column <- args[2]
top_fraction1 <- as.numeric(args[3])
top_fraction2 <- as.numeric(args[4])
use_absolute <- tolower(args[5]) %in% c("true", "1", "yes")
plot_width <- as.numeric(args[6])
plot_height <- as.numeric(args[7])
title <- args[8]
output_table1 <- args[9]
output_table2 <- args[10]
output_plot <- args[11]

data <- read.table(input_file, header = TRUE)

score_used_col <- if (use_absolute) {
  abs_col <- paste0("abs_", score_column)
  data[[abs_col]] <- abs(data[[score_column]])
  abs_col
} else {
  score_column
}

data_sorted <- data[order(-data[[score_used_col]]), ]

top_n1 <- as.integer(nrow(data) * top_fraction1)
if (top_n1 < 1) stop("top_fraction1 too small; no candidates selected.")
top_candidates1 <- data_sorted[1:top_n1, ]
write.table(top_candidates1, file = output_table1, quote = FALSE, sep = "\t", row.names = FALSE)

top_n2 <- as.integer(nrow(data) * top_fraction2)
if (top_n2 < 1) stop("top_fraction2 too small; no candidates selected.")
top_candidates2 <- data_sorted[1:top_n2, ]
write.table(top_candidates2, file = output_table2, quote = FALSE, sep = "\t", row.names = FALSE)

threshold1 <- min(top_candidates1[[score_used_col]])
threshold2 <- min(top_candidates2[[score_used_col]])

png(output_plot, width = plot_width, height = plot_height, units = "px")
manhattan(data, p = score_used_col, logp = FALSE, genomewideline = threshold1,
          suggestiveline = threshold2, col = c("#56B4E9", "#F0E442"),
          ylab = if (use_absolute) paste0("|", score_column, "|") else score_column,
          main = title)
dev.off()
