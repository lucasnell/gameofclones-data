
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
│   ├── 02-field-maps.R
│   ├── 03-field-monitoring-analysis.R
│   └── 04-field-disturbance-analysis.R
├── 02-prelims
│   ├── 01-alate-production.R
│   ├── 02-prelim-assays.R
│   └── 03-stage-structured-table.R
├── 03-experiments
│   ├── 01-exper-sims.R
│   ├── 02-experiments.R
│   └── 03-exper-sims-stoch.R
├── 04-stability
│   ├── 00-stability-shared.R
│   ├── 01-aphid-dispersal.R
│   ├── 02-trajectories-experiments.R
│   ├── 03-stationary-points.R
│   └── 04-paras-disp-hetero.R
│   └── 05-paras-disp-hetero-no-aphid.R
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

* `01-field-time-series.R` creates figure 1A.
  This plot shows proportion of aphids parasitized and infected with
  *Hamiltonella defensa* through time at our research site from 2011--2019.

* `02-field-maps.R` creates figure 1B where we plot parasitism
  through space at our research site.

* `03-field-monitoring-analysis.R` creates figures S1--S3 and table S1 that
  relate to analysis of parasitism in field montoring data.

* `04-field-disturbance-analysis.R` creates figure S4 and tables S2--S4 that
  relate to analysis of parasitism in field experiment data.


## 2. Preliminary information

The second group (folder `02-prelims`) relates to preliminary information
we used for simulations and/or experiments.

* `01-alate-production.R` looks for affects of clonal line and density on the
  proportion of aphids that are alates, using observations on colonies
  we maintain in the lab. These results are described in the "Simulations" 
  section in the materials and methods.

* `02-prelim-assays.R` creates figures S9 and S10 that show differences
  differences among our experimental clones in how well they compete without
  parasitoids and how well they resist parasitism. These scripts also contain
  the analyses testing for these differences.

* `03-stage-structured-table.R` creates table S7
  ('Values used in simulations for stage-structured aphid model parameter')
  in LaTeX format.


## 3. Experiments and their simulations

The third group (folder `03-experiments`) relates to the experiments and the
simulations we did beforehand to develop predictions for them.


* `01-exper-sims.R` creates figures 2A, 2B, and S6--S8.
  These are the *a priori* simulations of the experiments.
  They include how dispersal, starting resistance, and later perturbations
  affect the outcomes of the *in silico* experiment.

* `02-experiments.R` creates figure 2C, 2D, and S11--S13. These plots show
  the empirical results from the experiments.

* `03-exper-sims-stoch.R` creates figure 2E and 2F.
  These plots show *post hoc*, stochastic simulations of the experiments.


## 4. Stability and equilibria

The fourth group (folder `04-stability`) explores the stability and equilibrium
outcomes of the system, with some trajectories shown to illustrate how and 
why stability/equilibrium changes occur.
These were written by Anthony R. Ives.

* `00-stability-shared.R` contains code that is shared among at least two
  `04-stability` files.

* `01-aphid-dispersal.R` makes figures 3A, S16, and S17.
  These relate to how aphid dispersal changes the stability and equilibrium 
  states of the system using simulations that match the experiments.

* `02-trajectories-experiments.R` creates figure S14:
  'Example trajectories for the model parameterized for the experiment (fig. 2)
  in which parasitoids occur only in patch 1.'

* `03-stationary-points.R` creates figure S15:
  'Illustration of three stationary points in the model parameterized for the
  lab experiment.'

* `04-paras-disp-hetero.R` creates figures 3B, S18, S19, and S21.
  These relate to how parasitoid dispersal heterogeneity (both 
  aphid-density-dependent and aphid-density-independent sources) changes the 
  stability and equilibrium states of the system using simulations 
  that match the field.

* `05-paras-disp-hetero-no-aphid.R` creates fig. S20.
  This figure is the same as fig. S19 except that aphid-density-dependent 
  parasitoid dispersal was set to zero.



