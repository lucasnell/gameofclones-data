# `data-raw` directory

This directory contains the raw data from the experiments and some preliminary
assays. 

Folder contents:

```
.
├── Arlington.geojson
├── README.md
├── alate-counts.csv
├── assays-population_growth.csv
├── assays-wasp_resistance_choice.csv
├── assays-wasp_resistance_no_choice.csv
├── disturbance.csv
├── experiments-main.csv
├── experiments-pesky_wasps.csv
├── hamiltonella-2012-2019.csv
└── parasitism-2011-2019.csv
```


Individual file descriptions are below.


## `alate-counts.csv`

Surveys of winged vs non-winged aphids in aphid colonies maintained in
the lab.

Column descriptions:

* `line`: name of aphid clonal line
* `sample_date`: date we sampled
* `date`: date the colony was started
* `apterous`: number of non-winged aphids
* `alate`: number of winged aphids



## `Arlington.geojson`

Outlines of fields at Arlington Agricultural Research Station,
Wisconsin, USA. L. Nell outlined these fields by hand on Google Earth,
then converted them to GeoJSON.


## `assays-population_growth.csv`

Competition assays between the two aphid clones used in this study.
For these, we added 4 fourth-instar aphids of each line
into a pot containing three fava bean plants and counted the number of
each line twice per week until the plants died. 

Column descriptions:

* `year`: year sampled
* `month`: month sampled
* `day`: day of month sampled
* `obs`: initials of observer that did counting
* `rep`: assay repetition
* `line`: name of clonal line
* `num`: number of aphids
* `max_temp`: maximum temperature over the preceding week



## `assays-wasp_resistance_choice.csv`

Results from the choice experiment that assayed the two clones' relative
resistance to parasitism.
In these assays, we added 10 juveniles
of each aphid line into a deli container with two fava bean leaves.
We then added 3 female wasps to the containers for about 2 hours.
These wasps were mated but had not been exposed to aphids because
successful attacks on aphids can change wasps’ future aphid-color
preferences (Langley et al. 2006).
After wasp exposures, we put the aphids into petri dishes separated by
clonal line, each dish containing two fava bean leaves with bases
inserted into agar gel.
We kept these at 20ºC for 10 days, after which we counted the number
of mummies and living aphids.

Column descriptions:

* `start-year`: year when the exposures first occurred
* `start-month`: month when the exposures first occurred
* `start-day`: day of month when the exposures first occurred
* `end-year`: year when response variables were measured
* `end-month`: month when response variables were measured
* `end-day`: day of month when response variables were measured
* `wasp_group`: factor indicating which group of wasps was used
* `line`: name of clonal aphid line
* `juv-exposed`: number of juveniles exposed to wasps
* `juv-assayed`: number of juveniles surviving exposure (sometimes juveniles 
  were crushed by leaves in the containers or otherwise disappeared during 
  the exposure)
* `adult-end`: adults counted on the sampling date
* `juvenile-end`: juveniles counted on the sampling date
* `mummy-end`: mummies counted on the sampling date


## `assays-wasp_resistance_no_choice.csv`

Results from the no choice experiment that looked for differences among clones
in resistance to parasitism.
The only difference from the choice assays is how aphids were exposed to 
parasitoids:
We first exposed 10 juvenile aphids of a single line to parasitoids for about 3
hours, then exposed 10 juveniles of the other line to the same parasitoids for
the same duration. We did this for both aphid lines being the first to
be exposed.

The columns are the same as in `assays-wasp_resistance_choice.csv`, 
except for `round`, which indicates whether the parasitoids were exposed
to that particular line first (value of `1`) or second (`2`).


## `disturbance.csv`

Field experimental dataset from
[Kishinevsky and Ives (2022)](https://doi.org/10.1002/ecs2.4050)
on disturbance and how it affects pea aphids, predators, parasitoids, and
hyperparasitoids.

Column descriptions:


* `Day_of_sampling`: Integer day when sampling occurred, where day zero is 
  the first sampling day.
* `date`: Date of sampling.
* `field_tr`: Number of field with character indicating treatment.
* `aphids_1sweep`: Number of pea aphids per sweep.
* `total_para`: Number of dissected pea aphids that were parasitized.
* `dissected`: Total number of dissected pea aphids.
* `logAphids`: Log-transformed pea aphids per sweep.
* `wasp`: Total Aphidius ervi wasps across all sweeps.
* `wasps_1sweep`: Number of Aphidius ervi wasps per sweep.
* `sweeps.wasps`: Number of sweeps used to collect Aphidius ervi wasps.
* `trt`: The treatment implemented where the sample was collected
  ("Strips", "Insecticide", or "Control").





## `experiments-main.csv`

The main dataset for the experiments.

Column descriptions:

* `treatment`: DISPERSAL or ISOLATED; former indicates dispersal between cages, 
  latter indicates no dispersal between cages
* `rep`: repetition number for the set of cages
* `start_date`: date at which the cages were started
* `red_line`: which red line was used? (should be WIA-5D)
* `green_line`: which green line was used? (should be UT3)
* `cage`: WASP or NO WASP cage
* `date`: the date the observations were made
* `observer`: the observer's initials
* `plant<i>_red`: number of red aphids on plant `i`; values of -999 
  indicate these counts were intentionally skipped
* `plant<i>_green`: number of green aphids on plant `i`; values of -999 
  indicate these counts were intentionally skipped
* `wasps`: number of wasps in the cage; values of -999 indicate these counts 
  were intentionally skipped
* `mummies`: number of mummies in the cage; values of -999 indicate these
  counts were intentionally skipped
* `alates_total_red`: total red alates in the cage on up to 25% of the 
  tops of plants and the tops/sides
* `alates_total_green`: total green alates in the cage on up to 25% of the 
  tops of plants and the tops/sides
* `alates_in_red`: total number of red alates input from the other 
  cage of the pair
* `alates_in_green`: total number of green alates input from the other 
  cage of the pair
* `replaced_plants`: comma-separated (no spaces) list of the plant numbers
  that were replaced
* `wasps_rm`: number of adult wasps removed from the cage
* `mumm_rm`: number of mummies removed from the cage (should only happen for
  no wasp cages containing pesky wasps)
* `notes`: other useful information


## `experiments-pesky_wasps.csv`


Data on pesky wasps (wasps present in no-wasp cages) from experiments.

Column descriptions:

* `rep`: repetition of cage pair
* `date`: date of observation
* `time`: time of observation; only included when we checked multiple times per
  day.
* `mummies`: number of mummies found and removed
* `adult females`: number of adult wasps found and removed that were 
  identified as female
* `adult males`: number of adult wasps found and removed that were 
  identified as male
* `adults unk.`: number of adult wasps found and removed that were not 
  identified by sex
* `notes`: other information, such as removals that occurred without recording
  numbers



## `parasitism-2011-2019.csv`

These are raw data on parasitism from field surveys.
Data from 2011-2016 are from 
[Ives et al. (2020)](https://doi.org/10.1038/s41559-020-1155-0).

Column descriptions:

* `Year`: year of sampling date
* `Month`: month of sampling date
* `Day`: day of month of sampling date
* `Field`: field where sampling was done
* `Height_In`: height of alfalfa in inches
* `Cycle`: harvest cycle when sampling occurred
* `GrUnpara`: green unparasitized aphids
* `GrPara`: green parasitized aphids
* `RedUnpara`: red unparasitized aphids
* `RedPara`: red parasitized aphids
* `SweepCount`: number of sweeps used to collect pea aphids
* `SweepsPeaAphids`: number of pea aphids counted from sweeps



## `hamiltonella-2012-2019.csv`

These are raw data on *Hamiltonella defensa* infection from field surveys.
Data from 2012-2016 are from 
[Ives et al. (2020)](https://doi.org/10.1038/s41559-020-1155-0).

Column descriptions:

* `year`: year sampled
* `season`: season sampled
* `date`: date sampled
* `field`: field sampled
* `clone`: clone name (often new) sampled
* `ham`: 1 or 0 for whether this clone was infected with *H. defensa*


