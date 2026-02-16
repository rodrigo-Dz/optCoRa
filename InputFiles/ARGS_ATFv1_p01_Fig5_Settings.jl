
# --- Main arguments ---
iARG = (mm = "ATFv1_p01",       # Label for motif file
        ex = "Fig5",            # Label for parameters file
        pp = :mY,		# Label for perturbation type
        ax = :mY,               # Label for condition/environment
        an = "Sensitivity"		# Label for analysis type (Explore, Optimize, Dynamics, Curve, Sensitivity)
)

# --- Perturbation details ---
pert = (p   = iARG.pp,	        # Parameter to be perturbed
        d   = 1.05,		# Perturbation size (Delta rho)
        c   = iARG.ax,	        # Condition parameter
        r   = [-3.0,3.0],	# Range of conditions
        eps = 0.1,              # Threshold
        coras = 50,             # Number of conditions for each curve
)     

# --- Edit the analysis you want: ---
     
# Exploration details
expl  = (pOp  = [:mU,  :mW, :eP],	        # Parameters to optimize
        pMin = [-3.0, -3.0, -3.0],      # Min (log scale)
        pMax = [ 3.0,  3.0, 3.0],       # Max (log scale)   
        n_points = 236,                # number of points to evaluate (power of 2 recommended)
        prtD =1)		        # flag for printing full curve 

# Optimization details
opt  = (pOp  = [:mU,   :mW,  :eP ], #,  :gU,  :gW,  :e0,  :eM],	        # Parameters to optimize
        pMin = [-3.0, -3.0, -3.0],  # , -3.0, -3.0, -3.0, -3.0],      # Min (log scale)
        pMax = [ 3.0,  3.0,  3.0],  # ,  3.0,  3.0,  3.0,  3.0],         # Max (log scale)   
        iter = 5000,                   # number of points to evaluate
        cov = [0.1, 0.1, 0.1],      # ,0.1, 0.1, 0.1, 0.1],          # variance for each parameter
        M = 10,                         # mutation step size
        prtD =1,                        # flag for printing full curve
        rand = 1,                       # flag for random initial parameters
        mean = 50,                      # estabilizar la media de la distribucion (nMol), NaN para no estabilizarla  
        sd = 50)	

# Dynamics details
dyn = (tspan = (0.0, 200.0),   # time span for the simulation
        plot = 1,              # Index of variable to plot 
        pert_size = 1.5,       # Perturbation size (Delta rho)
        pert_time = 30.0       # Time of perturbation          
)


sens = (pSens = [:mU, :eP, :mW, :gU,  :gW,  :e0,  :eM],      # Parameters to analyze sensitivity
        r   = [-3.0,3.0],            # Range of conditions
        coras = 50,                  # Number of conditions for each curve
)