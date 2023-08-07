# READ ME 

This folder contains the codes and files used for the analyses in the "Changes in marine ecosystems after the end-Permian mass extinction are linked to oceanic habitats," 2023, paper by Al Aswad et al.

To re-produce the figures and analyses in Al Aswad et al. (2023), use the files in the following order:

"DownloadPbdb.r"
"Dataprep.r"
"Hexsim.r"
"Survival_PT.r"
"BC.r"
Certain data files are also available:

The Pbdb data that was downloaded using "DownloadPbdb.r" is available as "pbdb.jul2022.csv".

The hexagonal cells object, data.j, is created with "Dataprep.r "and uses a package ("DggridR") that is not compatible with Macbook M1 or M2 chips. In this case, the file is available here as "data.j.csv"

The shapefile used in "Dataprep.r" to create the figures is available here as "ma.252.Rda".
