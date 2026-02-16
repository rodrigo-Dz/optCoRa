using DataStructures

p = OrderedDict(
    :g   => 0.01,      # Dilution rate (e.g. [0.01,0.24] 1/min)
    :mY  => 0.125,     # Y maximum synthesis rate dependent of Y1 (nM/min)
	:kD  => 1.0,       # Activation threshold for Y synthesis dependent of Y1 (nM)
    :gY  => 0.1,       # Y degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mU  => 0.001,     # U synthesis rate dependent of Y (nM/min)
    :gU  => 0.0001,    # U degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mW  => 0.001,       # W constitutive synthesis rate (e.g. [0.005,0.02] nM/min)
    :gW  => 0.0001,    # W degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :e0  => 0.0001,    # U:W dissociation (unbinding) rate (e.g. [0.05,140] 1/min)
    :eP  => 0.001,    # U:W association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
    :eM  => 0.5,       # U:W (C) mutual anhilation/co-degradation (e.g. [0.03,2] 1/nM 1/min)
    :m0  => 1.25,      # Y0 synthesis rate dependent of W (nM/min)
    :m1  => 12.5,      # Y1 maximum synthesis rate dependent of Y0 (nM/min)
	:k1  => 1.0,       # Activation threshold for Y1 synthesis dependent of Y0 (nM)
    :mUs => NaN, 
)

# Initial conditions
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]     # Initial conditions for the ODE system