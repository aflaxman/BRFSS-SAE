This archive contains the small-area estimation (SAE) code developed
by Laura Dwyer-Lindgren at IHME (the same vintage as used for the
alcohol analysis), the first ten lines of the prepped BRFSS data, and
all of the other input files:

* merged_counties.rdata – a data frame mapping all FIPS that have
  existed to the merged fips set used for analysis (attached)

* census_regions_divisions.csv – a spreadsheet mapping states to
  census divisions (attached)

* county_adjancencies.rdata – a data frame listing every county and
  its neighbors (attached)

* covariates.rdata – a list in R where the first element is a
  data.frame of area-level covariates for males, and the second is
  area-level covariates for females (attached)

* brfss_microdata.rdata – a data frame with the prepped and compiled
  BRFSS data (first 10 lines attached as a csv)

Note that for all files except ‘merged_counties.rdata’, ‘fips’ means
merged fips (of which there are 3,127), not actual fips.

