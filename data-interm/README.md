# `data-interm` directory

This directory contains intermediate files containing simulation code 
that takes a while to run.

Folder contents:

```
.
├── README.md
├── paras-disp-hetero-traj.csv
├── paras-disp-hetero.csv
├── stable-sims-aphid_d.csv
└── stable-sims-wasp_d.csv
```


* `paras-disp-hetero-traj.csv`: Simplified dataset generated in simulations to
  test how parasitoid dispersal heterogeneity affects stability. See
  `scripts/04-stability/paras-disp-hetero.R` for more.
* `paras-disp-hetero.csv`: Full dataset generated in simulations to test how
  parasitoid dispersal heterogeneity affects stability. See
  `scripts/04-stability/paras-disp-hetero.R` for more.
* `stable-sims-aphid_d.csv`: Simulations of how aphid dispersal affect the
  final states of the system.
  See `scripts/04-stability/stable-sims-aphid_d.R` for more.
* `stable-sims-wasp_d.csv`: Simulations of how wasp dispersal affect the
  final states of the system.
  See `scripts/04-stability/stable-sims-wasp_d.R` for more.


