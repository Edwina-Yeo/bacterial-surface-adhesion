# Surface adhesion of dilute  of swimmers in shear flow for paper: 'A shear-induced limit on bacterial surface adhesion in fluid flow'
Authors: Edwina Yeo, contact edwina.yeo.14@ucl.ac.uk
Date: 05/25

#Dependencies:
Versions used: Python 3.12.4, Matlab R2023b, Julia version 1.10.4+0.aarch64.apple.darwin14. 
Packages and toolboxes: Julia: Random, LinearAlgebra, DelimitedFiles, Trapz, StatsBase, Matlab Optimization Toolbox, Python: 

To create data the appropriate bash processing script should be run, these scripts use main.jl passing the appropriate model parameters. 

#Processing scripts:
1. process_adhesion.sh: runs batched simulations to calculate data of adhesion of bacteria in varying wall shear rates using main.jl (Fig 1b and Fig 2c).
2. process_density.sh: runs the simulation required to create bacterial density plot Fig 2b.

#Computation scripts:
1. main.jl: Julia script which solves the each bacteria trajectory using Euler_maryumara solver returns position and time (xb,tb) if they collide with base.(Fig 1b and Fig 2c).
2. main_density.jl Julia script which solves the each bacteria trajectory using Euler_maryumara solver returning the position of each swimmer at time intervals (Fig 2b)
3. gen_data_adhesion.py: imports (xb,tb) datasets calculated using process-adhesion.sh from calculates the average adhesion and SD at varying shear rates


#Plotting scripts:
1. fig_adhesion.m takes flux and adhesion calculated by gen_data_adhesion.py and creates Fig1b and Fig2c
2. fig_density.py takes location data calculated by process.density.sh and creates Fig 2b. 
3. fig_regimes.m creates regime diagram Fig 2d
4. figure1-agent-based.svg combines Fig 1b with schematic Fig 1a. 
5. figure2-continuum.svg combines schematic Fig2a with plots Fig2b-d. 