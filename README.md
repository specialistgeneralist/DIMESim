# DIMESim - an agent-based model of mobilisation and radicalisation of social movements

This repository contains the code for **DIMESim**, a theoretically-
founded and empirically-informed agent-based model of mobilisation and radicalisation of social movements introduced and analysed in the paper:

&#x1F4D3; _From Mobilisation to Radicalisation: Isolating the Conditions for the Persistence and Radicalisation of Social Movements Using an Agent-Based Model_

&#x1F58B; By Emma F. Thomas, Mengbin Ye, Simon D. Angus, Tony J. Mathew, Winnifred Louis, Liam
Walsh, Silas Ellery, Morgana Lizzio-Wilson, and Craig McGarty.

<img src="figs/DIMESim.png" alt="Overview of the DIMESim model." style="display: block; margin-left: auto; margin-right: auto; max-width: 400px;">

## Requirements
DIMESim has been tested on MATLAB R2021a and requires no additional toolboxes. However, to run experiments in parallel using `SimRunner_par2()` (see below), you will require the  Parallel Computing Toolbox (`distrib_computing_toolbox`). However, the function checks for this at run time and reverts to serial processing if the toolbox is not available.

You can check if you have the Parallel Computing Toolbox by running the following command in MATLAB:
```matlab
>> license('checkout', 'distrib_computing_toolbox')
```
A response value of `1` indicates that you have the toolbox.

## Initalisation
In Matlab, navigate to the DIMESim director (using `cd`) and run the following command to initialise the model:
```matlab
>> setup
```
This will add the necessary paths.

To run an experiment, navigate to the `results` directory,
```matlab
>> cd results
```

## Running experiments

### A single model run
The main model function, `runModel.m` is found in the `model/` folder. To test the installation and run a single runfile (no experiments, no replicates) read the `constant_runfile.txt` and run it with `runModel.m` as follows:
```matlab
>> IN = read_constant_runfile('constant_runfile.txt');
>> res = runModel(IN);
```
this should not take more than 1s. The output `res` is a structure containing the model output. If you made no changes to the `constant_runfile.txt`, then you will have run the model with n=100 agents (rows) and 1000 timesteps (cols):
```matlab
>> res
res = 
  struct with fields:

                U: [1x1001 single]   .. Global broadcast signal (1=failure, -1=success)
       B_after_IR: [100x1001 single] .. Signal after individual reframing (1=failure, -1=success)
       B_after_CR: [100x1001 single] .. Signal after collective reframing (1=failure, -1=success)
                D: [100x1001 single] .. Disidentification values (in [0,100])
                I: [100x1001 single] .. Innovation values (in [0,100])
                M: [100x1001 single] .. Moralisation values (in [0,100])
                E: [100x1001 single] .. Energisation values (in [0,100])
                A: [100x1001 single] .. Action intention values (1=active, 0=inactive)
                C: [100x1001 single] .. Tactical orientation values (-1=radical, 1=conventional)
               xh: [100x1001 single] .. Last active action (-1=radical, 1=conventional)
                x: [100x1001 single] .. Current action (-1=radical, 1=conventional, 0=inactive)
           t_stop: 1000
stopping_cond_met: 0
        adjMatrix: [100x100 double]
```
To plot, for example, the time-series of average agent actions, you could run:
```matlab
>> plot(mean(res.x))
```


### Replication experiments 

To run factorial experiments across parameters as shown in the paper, we use the `SimRunner_par2()` function as a wrapper around `runModel.m`.

To generate the figures in the paper, navigate to the `results` directory and run one of the following commands as follows (qualitative runtimes are noted):
* Fig. 2a (quick): `SimRunner_par2('NetworkEvolution_C.txt')`
* Fig. 2b (quick): `SimRunner_par2('NetworkEvolution_LR.txt')`
* Fig. 2c-d (short): `SimRunner_par2('Timeseries_CDominant.txt')`
* Fig. 2e-f (short): `SimRunner_par2('Timeseries_LRDominant.txt')`
* Fig. 3-4, S1 (medium): `SimRunner_par2('ParameterSweep_p-F.txt')`
* Fig. 5, S2 (long): `SimRunner_par2('ParameterSweep_nu-R_p-F.txt')`
* Fig. S3-S4 (medium): `SimRunner_par2('ParameterSweep_nu-R.txt')`
* Fig. S5a,c (short): `SimRunner_par2('InitialActions_CScenario.txt');`
* Fig. S5b,d (short): `SimRunner_par2('InitialActions_LRScenario.txt')`
* Fig. S6 (quick): `NetworkTests()`

After running one of these experiments, `.mat` files will be saved in the `results` directory. To generate the plots in the paper, run the following commands (note that main text figures have been stylized):
* Fig. 2a: `res_timeseries_plot('NetworkEvolution_C.mat', [], true, true)`
* Fig. 2b: `res_timeseries_plot('NetworkEvolution_LR.mat', [], true, true)`
* Fig. 2c-d: `res_timeseries_plot('Timeseries_CDominant.mat', [], true, true)`
* Fig. 2e-f: `res_timeseries_plot('Timeseries_LRDominant.mat', [], true, true)`
* Fig. 3: `res_contourf([], [], [], 'paper_figure', 3)`
* Fig. 4: `res_contourf([], [], [], 'paper_figure', 4)`
* Fig. 5: `res_contourf([], [], [], 'paper_figure', 5)`
* Fig. S1: `res_contourf([], [], [], 'paper_figure', -1)`
* Fig. S2: `res_contourf([], [], [], 'paper_figure', -2)`
* Fig. S3: `res_contourf([], [], [], 'paper_figure', -3)`
* Fig. S4: `res_contourf([], [], [], 'paper_figure', -4)`
* Fig. S5a,c: `res_bargraph('InitialActions_CScenario.mat','initial_action')`
* Fig. S5b,d: `res_bargraph('InitialActions_LRScenario.mat','initial_action')`
* Fig. S6: `NetworkTests()`

## Acknowledgements
We make use of the following packages to run and analyse our model:
* **Distinguishable Colors**: Simple colour generator for maximum perceptual distinction.
  - Author: Timothy E. Holy (Washington) ([@timholy](https://github.com/timholy))
  - Ref: https://au.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
* **Hatchfill**: Fills an area with hatching or speckling.
    - Author: Takeshi Ikuma (LSU) ([@hokiedsp](https://github.com/hokiedsp))
    - Ref: https://github.com/hokiedsp/matlab-hatchfill2/
* **Matlab Axis Label Alignment tools**: Align axis labels with axes in 3d plots.
    - Author: Ligong Han (Rutgers) ([@phymhan](https://github.com/phymhan))
    - Ref: https://github.com/phymhan/matlab-axis-label-alignment
* **SimRunner**: Scientific full-factorial experiment wrapper for numerical simulation.
  - Author: Simon D. Angus (Monash) ([@specialistgeneralist](https://github.com/specialistgeneralist))
  - Ref: https://github.com/specialistgeneralist/SimRunner


## Cite
-- BibTex entry forthcoming --
