using DelimitedFiles;
using Distributions;
using CSV;
using Sobol;
using GaussianProcesses;
using ModelingToolkit;

# Definir parámetros de entrada iARG:
iARG = (mm = "FADv1",       # Label for motif file
        ex = "Fig1",        # Label for parameters file
        pp = :mY,           # Label for perturbation type
        ax = :mY,           # Label for condition/environment
        an = "OptCoRa");

# Cargar funciones y parámetros:
@independent_variables t  # Definir la variable independiente t para el modelo
mm = include(string("Library/Md_", iARG.mm, ".jl"));
fn = include(string("Library/FN_CoRa_v2.jl"));
include(string("InputFiles/ARGS_", iARG.mm, "_Pert_", iARG.ex, ".jl")); # Detalles de perturbación
include(string("InputFiles/ARGS_", iARG.mm, "_Par_", iARG.ex, ".jl"));  # Parámetros principales
include(string("InputFiles/ARGS_", iARG.mm, "_OptCoRa_", iARG.ex, ".jl"));

# Definir parámetros del modelo:
n = 30;
gap_size = 20;
gap_tol = 0.5;
rtol = 1e-10;  # Corregir rtol (era 1-e10, ahora 1e-10)

ranP = 0;
if(ranP==1)
    fn.mrwR(mrw,p);
end

p0 = [p[:mU], p[:mW], p[:eP]]
#fn.mrwR(mrw,p);
p1 = [p[:mU], p[:mW], p[:eP]]
#fn.mrwR(mrw,p);
p2 = [p[:mU], p[:mW], p[:eP]]
#p0 = [151.99240991582812, 999.9500652357119, 374.98127446339197]

x0 = zeros(length(mm.odeFB.syms));
strict = true;  # Definir 'strict'
eps = 0.2;      # Definir 'eps'
ruN = 1

open(string("OUT_OptCoRa_", iARG.mm, "_", iARG.ex, "_", iARG.pp, "_", iARG.ax, "_", ruN, ".txt"), "w") do io
    writedlm(io, [vcat("Run", "Iteration", [string(param) for param in mrw.pOp], 
                      string("|CoRa<=", pert.eps, "|"), "min(CoRa)", 
                      10 .^ collect(pert.r[1]:((abs(pert.r[1]) + abs(pert.r[2]))/n):pert.r[2]))], '	')

    # Implementación del algoritmo Nelder-Mead
    function nelder_mead(objective_function, initial_simplex, 
                         alpha=1.0, gamma=2.0, rho=0.5, sigma=0.5, max_iter=1000, tol=1e-6)
        n = size(initial_simplex, 2) - 1  # Dimensión del problema
        simplex = initial_simplex     # Inicializar el simplex
        values = [objective_function(simplex[:, i]) for i in 1:n+1]  # Valores iniciales

        for iter in 1:max_iter
            # Ordenar el simplex según los valores de la función objetivo
            perm = sortperm(values)
            simplex = simplex[:, perm]
            values = values[perm]

            # Centroid del simplex sin incluir el peor punto
            centroid = mean(simplex[:, 1:end-1], dims=2)

            # Reflexión
            xr = centroid .+ alpha .* (centroid .- simplex[:, end])
            fr = objective_function(xr)

            if fr < values[1]  # Mejor que el mejor
                xe = centroid .+ gamma .* (xr .- centroid)
                fe = objective_function(xe)
                if fe < fr
                    simplex[:, end] = xe
                    values[end] = fe
                else
                    simplex[:, end] = xr
                    values[end] = fr
                end
            elseif fr < values[end-1]  # Mejor que el segundo peor
                simplex[:, end] = xr
                values[end] = fr
            else
                if fr < values[end]
                    xc = centroid .+ rho .* (xr .- centroid)
                else
                    xc = centroid .+ rho .* (simplex[:, end] .- centroid)
                end
                fc = objective_function(xc)
                if fc < values[end]
                    simplex[:, end] = xc
                    values[end] = fc
                else
                    for i in 2:n+1
                        simplex[:, i] = simplex[:, 1] .+ sigma .* (simplex[:, i] .- simplex[:, 1])
                        values[i] = objective_function(simplex[:, i])
                    end
                end
            end

            # Convergencia
            if maximum(abs.(values .- mean(values))) < tol
                break
            end
        end

        return (simplex[:, 1], values[1])  # Mejor punto y valor
    end

    # Definir función objetivo para optimizar
    function objective_function(p0)
        p[:mU] = p0[1]
        p[:mW] = p0[2]
        p[:eP] = p0[3]
        p[:mUs] = NaN
        x0 = zeros(length(mm.odeFB.syms));
        CoRa1 = fn.CoRaCurve(p, pert, mm, n, gap_size, gap_tol, x0, rtol, strict)
        
        #if isempty(CoRa1) || all(isnan, CoRa1)
        #    return Inf  # Penalización alta para soluciones inválidas
        #end
        
        if isempty(filter(!isnan, CoRa1))
            return Inf  # Retornar un valor alto si no hay valores válidos en CoRa1
        end
        #valid_CoRa1 = filter(!isnan, CoRa1)
        #pavove = sum(valid_CoRa1 .>= eps) / length(valid_CoRa1)
        nas = count(isnan, CoRa1)
        #+nas
        #CoRa1 = [isnan(x) ? 0.0 : x for x in CoRa1]

        pavove = (sum(filter(!isnan, CoRa1) .>= eps)) / length(CoRa1)
        op, mi, cost = fn.CoRam(CoRa1, eps)
        writedlm(io, [vcat(ruN, 0, [p[i] for i in mrw.pOp], op, mi, CoRa1)], '	')
        flush(io)
        penalty = 0
        for i in 1:length(p0)
            if p0[i] < 0.001 || p0[i] > 1000
                penalty = Inf # Penalización alta para forzar la restricción
            end
        end
        return pavove + (cost) + nas*1000000 + penalty;  # Optimizar el mínimo CoRa
        #return pavove + (1/n * cost) +  penalty;  # Optimizar el mínimo CoRa
    end

    # Configurar simplex inicial
    #initial_simplex = hcat(p0, p1, p2)
    initial_simplex = hcat(p0, p0 .+ 0.1 .* p0, p0 .+ 0.2 .* p0)
    # Ejecutar Nelder-Mead
    best_params, best_value = nelder_mead(objective_function, initial_simplex)

    println("Resultado de la optimización: ", best_value)
    println("Mejor conjunto de parámetros encontrado: ", best_params)
end
