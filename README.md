# A Nonstationary Soft Partitioned Gaussian Process Model via Random Spanning Trees

## Data

## Code

## Instructions for use

All workflow information to reproduce key results of the paper is contained in `reproduce.Rmd`. The steps for running the workflow are:

**Step 1**: Install all required R packages in `dependencies.txt`.

**Step 2** (Optional): Compile the source code for the BRISC library by running the following command in your terminal:
```
R CMD SHLIB BRISC/BRISC.cpp BRISC/BRISC_Predict.cpp BRISC/util.cpp BRISC/init.c BRISC/lbfgs.c
```
This step is optional and only required if you want to re-run SPGP model fitting and prediction.

**Step 3** (Optional): Modify the options in Line 19-23 of `reproduce.Rmd` accordingly. By default, key figures and tables are reproduced from pre-computed model outputs in `data/`. To re-generate model outputs, make the following change(s):
* To re-generate simulation input data, set `gen_sim_input = TRUE`.
* To re-run SPGP model in simulation, set `fit_SPGP_sim = TRUE`.
* To re-run TGP model in simulation, set `fit_TGP_sim = TRUE`.
* To re-run SPGP model in precipitation data application, set `fit_SPGP_app = TRUE`.
* To re-run TGP model in precipitation data application, set `fit_TGP_app = TRUE`.

Please note that re-running these models can take several hours.

**Step 4**: Run the workflow in `reproduce.Rmd` by compiling the R Markdown file. For example, this is can be done by running the following command in your terminal:
```
Rscript -e "rmarkdown::render('reproduce.Rmd')"
```
Reproducible results can be found in the generated `reproduce.html` file.
