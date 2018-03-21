# simul_shad_spawn
This project aims at estimating the number of shad spawners from splashes recording on spawning grounds.
The approach is ABC.
The mechanistic model simulates individual shad spawning decisions:
- number of spawnings
- temporal distribution (within season and within night)
- spatial distribution
This model is informed by data collected in the Nivelle in 2017 and 2018, using radiotracking and accelerometers.
The model is run many times using a prior distribution for the number of spawners in the population. On each simulation, the actual sampling scheme is applied, summary statistics are computed and compared to summary statistics on actual data in order to infer the posterior distribution of the number of spawners.
The ultimate aim is to wrap all this in a Shiny app which could be used by the practitioners (Migado, Migradour...).
