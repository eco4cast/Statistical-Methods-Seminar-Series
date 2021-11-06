xx <- read.csv("Bolker_PollEV_Questions-1.csv")
names(xx) <- c("Q", "A", "junk")
xx <- xx[nchar(xx$A)>0,]
outfile <- "QA.md"


hdr <- "---
title: 'Questions and answers'
---
"

cat(hdr, "\n", file = outfile, append = FALSE)
for (i in 1:nrow(xx)) {
  cat(file = outfile, append = TRUE,
      i, ". ",
      xx$Q[i], "\n\n",
      "**Answer**: \n",
      xx$A[i], "\n\n")
}
