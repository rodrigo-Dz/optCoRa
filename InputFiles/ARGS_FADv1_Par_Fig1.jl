# Kinetic parameters
p = Dict([
    :g   => 0.0001,      # Dilution rate (e.g. [0.01,0.24] 1/min)
    :mY  => 0.125,       # Y synthesis rate dependent of W (nM/min)
    :gY  => 1,       # Y degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mU  => 0.002659529,    # U synthesis rate dependent of Y (nM/min)
    :gU  => 0.0001,      # U degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mW  => 748.209878,      # W constitutive synthesis rate (e.g. [0.005,0.02] nM/min)
    :gW  => 0.0001,    # W degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :e0  => 0.0001,    # U:W dissociation (unbinding) rate (e.g. [0.05,140] 1/min)
    :eP  => 430.317554,    # U:W association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
    :eM  => 0.5,       # U:W (C) mutual anhilation/co-degradation (e.g. [0.03,2] 1/nM 1/min)
    :mUs => NaN,       # LOCAL: U constitutive synthesis rate in the locally analogous system (mU*Yss)
]);

#Inital conditions
x0 = zeros(length(mm.odeFB.syms));
