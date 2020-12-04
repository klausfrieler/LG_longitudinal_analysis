# Longgold Longitudinal Analysis
Longitudinal analysis workbench for the LongGold-Project (Development of musicality in adolescence).

#### Measurement Error Simulations
To run the script for simulating measurement error open a console and enter

```
Rscript --vanilla run_simulations.R [-n N] [-s [K[:L],...]] [-o OUTDIR] [-l LABEL] [-f [VAR=VALUE, ...]]
```

The option ```-n``` gives with ```N``` the number of simulated data points per parameter combination. 

The ```-s``` options allows to select simulation ids as a comma-separated list of ranges and single values, with a maximum value of 120. Needs quoting if more than one item is provided. Example: -s "1, 3:4, 7" will run simulations with parameter sets #1, #3, #4, and #7. If missing, the full range 1:120 will be included (if no filters are defined.)

The ```-o``` option gives in ```OUTDIR``` is the subdirectory to store the results, default is ```data/simulations```. Note that the output directory must already exist otherwise the script will stop with an error.

With the ```-l``` option you can provided a label for this simulation to be used in file names and output. Default is ```simu_<MIN_SIMUL_ID>_<MAX_SIMUL_ID>```.

You can also filter the simulation parameter sets using the ```-f``` option and providing a comma-separated list of variable/value pairs, e.g., ```size=250, error_level = medium``` specifies all simulations with samplings size 250 and a medium error level.  Usuables variables are ```size``` (allowed values: 250, 500, 1000, 2000, 4000), ```error_level``` (value "none", "small", "medium", "large"), ```with_na``` (value FALSE, TRUE) and ```method``` (value "mvt", "mvt-decor", "mvt-upcor"). Needs quoting, if containes more than one item. If filters and ids are both specified, the intersection of the two sets will be executed.
