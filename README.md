# MDS
Multidimensional Scaling implementation and comparison with the smacof package in R.

### MDS_lib.R
Contains all the functions that perform the MDS optimalization. The final configuration is found through the majorization of the Kruskal's stress function.

### MDS_main.R
Shows the functionality of the provided library and compares the results witht the smacof package in R

### MDS_report.Rmd
Describes the comparison of our implementation vs. the `smacof` package in the R-markdown format

### MDS_report.pdf
Pdf version of MDS_report.Rmd

### basket.Rdata
Contains a toy dataset of 12 grocery store items and their co-purchases (similarity matrix)
