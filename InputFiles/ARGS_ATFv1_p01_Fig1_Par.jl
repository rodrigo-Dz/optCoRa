using DataStructures

# Diccionario ordenado
p = OrderedDict(
    :g   => 0.0001,    
    :mY  => 0.125,     
    :gY  => 1,       
    :mU  => 0.125,      
    :gU  => 0.0001,    
    :mW  => 0.1,       
    :gW  => 0.0001,    
    :e0  => 0.0001,    
    :eP  => 0.0375,     
    :eM  => 0.5,       
    :mUs => NaN
)

u0 = [0.0, 0.0, 0.0, 0.0]     # Initial conditions for the ODE system