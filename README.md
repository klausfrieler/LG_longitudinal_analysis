# Longgold Longitudinal Analysis
Longitudinal analysis workbench for the LongGold-Project (Development of musicality in adolescence).

#### Measurement Error Simulations
To run the script for simulating measurement error open a console and enter

```
Rscript --vanilla run_simulations.R -n N -s K:L -o OUTDIR
```

where ```N``` is the number of simulated data points per parameter combination and ```K``` and ```L``` are integers with ```1<=K<=L<108```, indicating a subset of simulation parameters. ```OUTDIR``` is the subdirectory to store the results, default is ```data/simulations```. Note that the output directory must already exist otherwise the script will stop with an error.
