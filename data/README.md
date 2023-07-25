This directory contains precipitation and simulation data used in the paper.

## Precipitation data

File `CONUS_precip_west_18.csv` consists of average precipitation data aggregated from the daily data over October 1, 2017 to September 30, 2018. The daily data are originally from the Global Historical Climate Network-Daily database, preprocessed by Dr. Mark Risser, and publicly available at his [personal website](https://sites.google.com/site/markdrisser/data-sets?authuser=0). We only focus on observations at 1939 stations located to the west of 90 degrees west.

The columns used in the paper are:
* `lcc_x`: x-coordinate of Lambert conformal conic projection.
* `lcc_y`: y-coordinate of Lambert conformal conic projection.
* `precip`: average precipitation data (in mm/day).

## Simulation

Simulation data is provided in `sim_input.RData`.

## Model outputs

Pre-computed model outputs are stored in other `RData` files.