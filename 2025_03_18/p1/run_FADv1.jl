## Running in julia terminal
## julia -t 4
	#cd("C:\\MyLibrary\\Location\\")
	#using Pkg; Pkg.activate(".");
	using Pkg;
	using CSV;
	using DelimitedFiles;

	iARG = (mm = "FADv1",		# Label for motif file
			ex = "Fig1",		# Label for parameters file
			pp = :mY,			# Label for perturbation type
			ax = :mY,			# Label for condition/environment
			an = "OptCoRa");		# Chose analysis type (Options: ExSSs, ExDyn, CoRams, OptCoRa)
	include("optCoRa_main.jl")
 