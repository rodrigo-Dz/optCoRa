# --- Main script for CoRa analysis ---
# --- Rodrigo Aguilar
# --- September 2025


# --- Required libraries ---
using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using Plots
using NLsolve
using Sobol
using DelimitedFiles
using Statistics

# The code is structured in the following way:
# 1. Settings: Here you can edit the main arguments for the analysis, 
        # the perturbation details, and the details for each type of analysis 
        # (exploration, optimization, dynamics, curve, sensitivity).
# 2. Main functions: Here we include the main functions for the CoRa analysis, 
        # which are defined in the Library/FN_CoRa.jl file.
# 3. Model: Here is the model file, which is defined in the Library
        # folder. The should contain the ODEs of the systems and the 
        # function to compute the steady state for the analogous system.
# 4. Core parameters - Experiment: Here we include the core parameters for the model 


# --- Settings ---
               # ARGS_model_experiment_Setting.jl
Settings_file = "ARGS_FDPv1_p01_Fig2_Settings.jl"
include(string("./InputFiles/", Settings_file))

# --- Main functions ---
fn = include(string("Library/FN_CoRa.jl"))

# --- Model ---
mm = include(string("Library/Md_",iARG.mm,".jl"))

# --- Core parameters ---
include(string("InputFiles/ARGS_",iARG.mm,"_",iARG.ex,"_Par.jl"))

# ------

# We include 5 different types of analysis:
#  - Explore: Explore the parameter space to find regions of good performance.
#  - Optimize: Optimize the parameters to find the best performance.
#  - Dynamics: Simulate the dynamics of the system under perturbation.
#  - Curve: Generate the CoRa curve for a range of conditions.
#  - Sensitivity: Compute the sensitivity of the system to parameter changes. 

# To run a specific analysis, simply edit the "an" argument in the main
# arguments section of the settings file. The results will be saved in the Output folder.

if (iARG.an == "Explore")
        fn.explore(
                iARG,       # Main arguments     
                mm,         # Model
                p,          # Core parameters
                pert,       # Perturbation details
                expl,       # Exploration details
                u0)         # Initial conditions
elseif (iARG.an == "Optimize")
        fn.optimize(
                iARG,       # Main arguments     
                mm,         # Model
                p,          # Core parameters
                pert,       # Perturbation details
                opt,        # Optimization details
                u0)         # Initial conditions
elseif (iARG.an == "Dynamics")
        fn.dynamics(
                mm,         # Model
                u0,         # Initial conditions
                p,          # Core parameters
                pert,       # Perturbation details
                dyn,        # Dynamics details
                iARG)       # Main arguments
elseif (iARG.an == "Curve")
        fn.curve(
                p,          # Core parameters
                u0,         # Initial conditions
                mm,         # Model
                pert,       # Perturbation details
                iARG)       # Main arguments
elseif (iARG.an == "Sensitivity")
        fn.sensitivity(
                p,          # Core parameters
                u0,         # Initial conditions
                mm,         # Model
                pert,       # Perturbation details
                iARG,       # Main arguments
                sens)       # Sensitivity details   
else
        println("Error: Analysis type not recognized")
end




