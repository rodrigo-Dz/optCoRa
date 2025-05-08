###CoRa_v2.5 - Optimización Bayesiana para Minimizar Evaluaciones Costosas

using DelimitedFiles;
using Distributions;
using CSV;
using Sobol;
using GaussianProcesses;
using Optim;
using ModelingToolkit;


# Definir parámetros de entrada iARG:
iARG = (mm = "FADv1",       # Label for motif file
        ex = "Fig3",        # Label for parameters file
        pp = :mY,            # Label for perturbation type
        ax = :mY,            # Label for condition/environment
        an = "OptCoRa");

# Cargar funciones y parámetros:
@independent_variables t # Definir la variable independiente t para el modelo
mm = include(string("Library/Md_",iARG.mm,".jl"));
fn = include(string("Library/FN_CoRa_v2.jl"));
include(string("InputFiles/ARGS_",iARG.mm,"_Pert_",iARG.ex,".jl")); # Detalles de perturbación
include(string("InputFiles/ARGS_",iARG.mm,"_Par_",iARG.ex,".jl")); # Parámetros principales
include(string("InputFiles/ARGS_",iARG.mm,"_OptCoRa_",iARG.ex,".jl"));


# Definir parámetros del modelo:
n = 20;
gap_size = 10;
gap_tol = 0.5;
rtol = 1e-10; # Corregir rtol (era 1-e10, ahora 1e-10)

#p0 = [0.125, 0.1, 0.0375]; # Convertir pO a un vector

ranP = 0;
if(ranP==1)
    fn.mrwR(mrw,p);
end

p0 = [p[:mU], p[:mW], p[:eP]]

x0 = zeros(length(mm.odeFB.syms));
strict = true; # Definir 'strict'
eps = 0.2; # Definir 'eps'
ruN = 1
    open(string("OUT_OptCoRa_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,"_",ruN,".txt"), "w") do io
        writedlm(io, [vcat("Run","Iteration",[string(param) for param in mrw.pOp],string("|CoRa<=",pert.eps,"|"),"min(CoRa)",10 .^ collect(pert.r[1]:((abs(pert.r[1]) + abs(pert.r[2]))/n):pert.r[2]))], '\t')
        # Función objetivo para optimizar op1:
    function objective_function(p0)
        p[:mU] = p0[1]
        p[:mW] = p0[2]
        p[:eP] = p0[3]
        p[:mUs] = NaN
        CoRa1 = fn.CoRaCurve(p, pert, mm, n, gap_size, gap_tol, x0, rtol, strict);
        CoRa1 = replace(CoRa1, NaN => 0.0)
        pavove = sum(CoRa1 .>= eps)/length(CoRa1)
        op, mi, cost = fn.CoRam(CoRa1, eps);
        writedlm(io, [vcat(ruN,0,[p[i] for i in mrw.pOp],op,mi,CoRa1)],'\t')
        flush(io)
        penalty = 0
        for i in 1:length(p0)
            if p0[i] < 0.01 || p0[i] > 1000
                penalty += 1e6 # Penalización alta para forzar la restricción
            end
        end
        return pavove + (1/n * cost) + penalty;
        #return pavove + penalty
        #return cost;
        #return maximum(filter(!isnan,CoRa1)) 
    end

    # Optimizar la función objetivo para minimizar op1:
    result = optimize(objective_function, p0, NelderMead());
    println("holi??")
    end
    println("Resultado de la optimización: ", result);
    println("Mejor conjunto de parámetros encontrado: ", result.minimizer);
