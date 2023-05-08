
# `scripts` folder

This folder contains the R scripts containing all the analyses and generating
all the figures in the manuscript.

Folder contents: 

```
.
├── 00-shared-all.R
├── 01-field
│   ├── 00-field-shared.R
│   ├── 01-field-time-series.R
│   └── 02-field-maps.R
├── 02-prelims
│   ├── 01-alate-production.R
│   ├── 02-prelim-assays.R
│   └── 03-stage-structured-table.R
├── 03-experiments
│   ├── 01-exper-sims.R
│   └── 02-experiments.R
├── 04-stability
│   ├── 00-stability-shared.R
│   ├── 01-aphid-dispersal.R
│   ├── 02-wasp-dispersal.R
│   ├── 03-trajectories-experiments.R
│   ├── 04-stationary-points.R
│   ├── 05-trajectories-field-high-d.R
│   └── 06-trajectories-field.R
└── README.md
```


To run this code, first read the top-level `README.md` for `gameofclones-data`,
in which it will show you what R packages are required.
Next, check out file `scripts/00-shared-all.R` to set the objects `write_plots`
and `.n_threads`.
You should be set after that.

File `00-shared-all.R` contains code that is used by all scripts.
Other scripts are separated by folder into groups.
The individual scripts are named as `<order>-<description>`, 
where `<description>` can contain one or more hyphens.
The `<order>` field makes the scripts, when alphabetized, appear in the 
approximate order they show up in the manuscript.
There are four groups of scripts:


## 1. Field data

The first group (folder `01-field`) shows field data from 5–10 alfalfa fields at 
the Arlington Agricultural Research Station, Wisconsin, USA, 
through time and space.

* `00-field-shared.R` contains code that is shared among all `1-field`
  files, including reading and cleaning up the data frame containing 
  parasitism data.

* `01-field-time-series.R` creates figures 1A, 1B, and S14.
  These plots show proportion of aphids parasitized and infected with
  *Hamiltonella defensa* through time at our research site from 2011--2019.

* `02-field-maps.R` creates figure 1C where we plot parasitism
  through space at our research site.


## 2. Preliminary information

The second group (folder `02-prelims`) relates to preliminary information
we used for simulations and/or experiments.

* `01-alate-production.R` looks for affects of clonal line and density on the
  proportion of aphids that are alates, using observations on colonies
  we maintain in the lab. These results are described in the "Simulations" 
  section in the materials and methods.

* `02-prelim-assays.R` creates figures S4 and S5 that show differences
  differences among our experimental clones in how well they compete without
  parasitoids and how well they resist parasitism. These scripts also contain
  the analyses testing for these differences.

* `03-stage-structured-table.R` creates table S3
  ('Values used in simulations for stage-structured aphid model parameter')
  in LaTeX format.


## 3. Experiments and their simulations

The third group (folder `03-experiments`) relates to the experiments and the
simulations we did beforehand to develop predictions for them.


* `01-exper-sims.R` creates figures 2A, 2B, S2, and S3. These are the simulations
  I did to help understand the system before the experiments.
  They include how dispersal, starting resistance, and later perturbations
  affect the outcomes of the *in silico* experiment.

* `02-experiments.R` creates figure 2C, 2D, S6, and S7. These plots show
  the empirical results from the experiments.


## 4. Stability and equilibria

The fourth group (folder `04-stability`) explores the stability and equilibrium
outcomes of the system, with some trajectories shown to illustrate how and 
why stability/equilibrium changes occur.
These were written by Anthony R. Ives.

* `00-stability-shared.R` contains code that is shared among at least two
  `4-stability` files.

* `01-aphid-dispersal.R` makes figures 3A and S10. Both relate to how aphid dispersal
  changes the stability and equilibrium states of the system using simulations
  that match the experiments.

* `02-wasp-dispersal.R` makes figures 3B and S11. Both relate to how wasp dispersal
  changes the stability and equilibrium states of the system using simulations
  that match the field.

* `03-trajectories-experiments.R` creates figure S8:
  'Example trajectories for the model parameterized for the experiment (fig. 2)
  in which parasitoids occur only in patch 1.'

* `04-stationary-points.R` creates figure S9:
  'Illustration of three stationary points in the model parameterized for the
  lab experiment.'

* `05-trajectories-field-high-d.R` creates figure S13:
  'Consequences of increasing parasitoid dispersal on the persistence of the
  resistant clone.'

* `06-trajectories-field.R` creates figure S12:
  'For the model parameterized for the field data (fig. 3C,D), time
  trajectories of susceptible aphid clones (yellow), resistant clones (green),
  and parasitoids (purple).'


