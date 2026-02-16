using DataStructures

# Diccionario ordenado
p = OrderedDict(
    :g   => 0.0001,    
    :mY  => 0.1,     
    :gY  => 0.1,       
    :mU  => 0.0030,      
    :gU  => 0.0001,    
    :mW  => 839.1278,       
    :gW  => 0.0001,    
    :e0  => 0.0001,    
    :eP  => 0.0781,     
    :eM  => 0.5,       
    :mUs => NaN
)

u0 = [0.0, 0.0, 0.0, 0.0]     # Initial conditions for the ODE system