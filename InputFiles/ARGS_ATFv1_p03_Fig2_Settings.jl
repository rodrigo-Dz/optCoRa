
# --- Main arguments ---
iARG = (mm = "ATFv1_p03",       # Label for motif file
        ex = "Fig2",            # Label for parameters file
        pp = :mY,		# Label for perturbation type
        ax = :mY,               # Label for condition/environment
        an = "Optimize"		# Label for analysis type (Explore, Optimize, Dynamics, Curve)
)

# --- Perturbation details ---
pert = (p   = iARG.pp,	        # Parameter to be perturbed
        d   = 1.05,		# Perturbation size (Delta rho)
        c   = iARG.ax,	        # Condition parameter
        r   = [-3.0,3.0],	# Range of conditions
        eps = 0.1,              # Treshold
        coras = 30,             # Number of conditions for each curve
)     

# --- Edit the analysis you want: ---
     
# Exploration details
expl  = (pOp  = [:mU,:mW,:eP],	        # Parameters to optimize
        pMin = [-3.0, -3.0, -3.0],      # Min (log scale)
        pMax = [3.0,  3.0,  3.0],       # Max (log scale)   
        n_points = 2048,                # number of points to evaluate (power of 2 recommended)
        prtD =1)		        # flag for printing full curve 

# Optimization details
opt  = (pOp  = [:mU, :mW, :eP],	        # Parameters to optimize
        pMin = [-3.0, -3.0, -3.0],      # Min (log scale)
        pMax = [3.0, 3.0, 3.0],         # Max (log scale)   
        iter = 10000,                   # number of points to evaluate (power of 2 recommended)0
        cov = [0.1, 0.1, 0.1],          # variance for each parameter
        M = 10,                         # mutation step size
        prtD =1,                        # flag for printing full curve
        rand = 1,                       # flag for random initial parameters
        )	

# Dynamics details
dyn = (tspan = (0.0, 1000.0),   # time span for the simulation
        plot = 1,              # Index of variable to plot 
        pert_size = 2,       # Perturbation size (Delta rho)
        pert_time = 500.0       # Time of perturbation          
)
