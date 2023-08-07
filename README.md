# READ ME 

**To re-produce the figures and analyses in Al Aswad et al. (2023), use the files in the following order:**

1."DownloadPbdb.r"
2."Dataprep.r"
3."Hexsim.r"
4."Survival_PT.r"
5."BC.r"

Each file depends on one or more products from the previous files to run properly, so they are not standalone.

**Certain data files are also available:**

The Pbdb data that was downloaded using "DownloadPbdb.r" is available as "pbdb.jul2022.csv".

The hexagonal cells object, data.j, is created with "Dataprep.r "and uses a package ("DggridR") that is not compatible with Macbook M1 or M2 chips. In this case, the file is available here as "data.j.csv"

The shapefile used in "Dataprep.r" to create the figures is available here as "ma.252.Rda".
