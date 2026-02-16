using DataStructures

# Diccionario ordenado
p = OrderedDict(
    :g   => 0.001,      # Dilution rate (e.g. [0.01,0.24] 1/min)
    :mY  => 1.0,       # Y synthesis rate dependent of W (nM/min)
    :gY  => 1.0,       # Y degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mU  => 0.125,     # U synthesis rate dependent of Y (nM/min)
    :gU  => 0.0001,      # U degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mW  => 1.0,       # Maximum W synthesis rate dependent of W (nM/min)
    :kD  => 0.02,	   # W induction KD (nM)
    :gW  => 0.0001,    # W degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :e0  => 0.0001,    # U:W dissociation (unbinding) rate (e.g. [0.05,140] 1/min)
    :eP  => 0.0375,    # U:W association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
    :eM  => 0.5,       # U:W (C) mutual anhilation/co-degradation (e.g. [0.03,2] 1/nM 1/min)
    :mUs => NaN,       # LOCAL: U constitutive synthesis rate in the locally analogous system (mU*Yss)
)

u0 = [0.1, 0.1, 0.1, 0.1]     # Initial conditions for the ODE system