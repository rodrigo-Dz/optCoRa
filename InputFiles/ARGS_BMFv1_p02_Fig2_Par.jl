using DataStructures

# Diccionario ordenado
p = OrderedDict(
    :g   => 0.01,      # Dilution rate (e.g. [0.01,0.24] 1/min)
    :mY  => 0.001,       # Y synthesis rate dependent of W (nM/min)
    :gY  => 0.1,       # Y degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mA  => 40.57,    # A constitutive synthesis rate (e.g. [0.005,0.02] nM/min)
    :mB  => 755.5,     # B synthesis rate dependent of Y (1/min)
    :mU  => 431.2,       # U synthesis rate dependent of Y (nM/min)
    :e0  => 0.002924,    # U:W dissociation (unbinding) rate (e.g. [0.05,140] 1/min)
    :eP  => 0.03629,      # U:W association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
    :bA  => 0.03611,       # U to Us transition rate (1/(nM min))
    :bI  => 6.988,       # U:W association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
    :m0  => 1.25,      # Y0 synthesis rate dependent of W (nM/min)
    :m1  => 12.5,      # Y1 maximum synthesis rate dependent of Y0 (nM/min)
	:k1  => 1.0,       # Activation threshold for Y1 synthesis dependent of Y0 (nM)
    :kD  => 0.002569,       # Activation threshold for Y synthesis dependent of Y1 (nM)
    :mBs => 0.0,       # LOCAL: B constitutive synthesis rate in the locally analogous system (mB*Yss)
);

u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]     # Initial conditions for the ODE system