using DataStructures

# Diccionario ordenado
p = OrderedDict(
    :g   => 0.0001,    
    :mY  => 0.125,     
    :gY  => 0.1,       
    :mU  => 0.01,      
    :gU  => 0.0001,    
    :mW  => 100,       
    :gW  => 0.0001,    
    :e0  => 0.0001,    
    :eP  => 0.001,     
    :eM  => 0.5,       
    :mUs => NaN
)

u0 = [0.0, 0.0, 0.0, 0.0]     # Initial conditions for the ODE system