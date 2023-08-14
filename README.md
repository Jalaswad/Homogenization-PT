# READ ME 

**To re-produce the figures and analyses in Al Aswad et al. (in prep), use the files in the following order:**

1."DownloadPbdb.r"

2."Dataprep.r"

3."Hexsim.r"

4."Survival_PT.r"

5."BC.r"

Each file depends on one or more products from the previous files to run properly, so they are not standalone.

**Some data files are also available, in case you have any difficulty creating them:**

The hexagonal cells dataframe, data.j, is created with "Dataprep.r "and uses a package ("DggridR") that is not compatible with Macbook M1 or M2 chips. In this case, the file is available here as "data.j.csv"

The "252.scotese_ploygon" shapefiles used in "Dataprep.r" to create the figures is available here as "ma.252.Rda". These shapefiles were created from the Scotese (2016) shapefiles in the GPlates software in 2021.
