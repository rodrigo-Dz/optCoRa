using DelimitedFiles
using Distributions
using CSV
using Sobol
using GaussianProcesses
using ModelingToolkit

# Definir parámetros de entrada iARG:
iARG = (mm = "FADv1",       # Label for motif file
        ex = "Fig3",        # Label for parameters file
        pp = :mY,           # Label for perturbation type
        ax = :mY,           # Label for condition/environment
        an = "OptCoRa");

# Cargar funciones y parámetros:
@independent_variables t  # Definir la variable independiente t para el modelo
mm = include(string("Library/Md_", iARG.mm, ".jl"))
fn = include(string("Library/FN_CoRa_v2.jl"))
include(string("InputFiles/ARGS_", iARG.mm, "_Pert_", iARG.ex, ".jl")) # Detalles de perturbación
include(string("InputFiles/ARGS_", iARG.mm, "_Par_", iARG.ex, ".jl"))  # Parámetros principales
include(string("InputFiles/ARGS_", iARG.mm, "_OptCoRa_", iARG.ex, ".jl"))

# Definir parámetros del modelo:
n = 20
gap_size = 10
gap_tol = 0.5
rtol = 1e-10

ranP = 1
if ranP == 1
    fn.mrwR(mrw, p)
end

x0 = zeros(length(mm.odeFB.syms))
strict = true  # Definir 'strict'
eps = 0.2      # Definir 'eps'
ruN = 1

open(string("OUT_OptCoRa_", iARG.mm, "_", iARG.ex, "_", iARG.pp, "_", iARG.ax, "_", ruN, ".txt"), "w") do io
    writedlm(io, [vcat("Run", "Iteration", [string(param) for param in mrw.pOp], 
                      string("|CoRa<=", pert.eps, "|"), "min(CoRa)", 
                      10 .^ collect(pert.r[1]:((abs(pert.r[1]) + abs(pert.r[2]))/n):pert.r[2]))], '\t')

    # Definir función objetivo
    function objective_function(p0)
        p[:mU], p[:mW], p[:eP] = p0
        CoRa1 = fn.CoRaCurve(p, pert, mm, n, gap_size, gap_tol, x0, rtol, strict)

        if isempty(filter(!isnan, CoRa1))
            return Inf
        end

        pavove = (sum(filter(!isnan, CoRa1) .>= eps) + count(isnan, CoRa1)) / length(CoRa1)
        op, mi, cost = fn.CoRam(CoRa1, eps)
        writedlm(io, [vcat(ruN, 0, [p[i] for i in mrw.pOp], op, mi, CoRa1)], '\t')
        flush(io)

        penalty = sum(p0 .< 0.01) * Inf + sum(p0 .> 1000) * Inf
        return pavove + (1/n * cost) + count(isnan, CoRa1) * 10 + penalty
    end

    # Implementación de SAMIN
    function samin(objective_function, bounds, max_iters=1000, temp_init=1.0, alpha=0.9, max_attempts=100)
        dim = length(bounds)
        current_position = [rand(bounds[i][1]:0.1:bounds[i][2]) for i in 1:dim]  # Posición inicial
        current_score = objective_function(current_position)

        best_position = copy(current_position)
        best_score = current_score

        temp = temp_init  # Temperatura inicial

        # Para registrar la trayectoria
        trajectory = [copy(current_position)]
        scores = [current_score]

        for iter in 1:max_iters
            for attempt in 1:max_attempts
                # Generar nuevo candidato
                candidate = [current_position[i] + (rand() - 0.5) * (bounds[i][2] - bounds[i][1]) * temp
                             for i in 1:dim]
                # Asegurarse de que el candidato está en los límites
                candidate = [clamp(candidate[i], bounds[i][1], bounds[i][2]) for i in 1:dim]

                candidate_score = objective_function(candidate)

                # Aceptar el candidato con probabilidad basada en la temperatura
                Δ = candidate_score - current_score
                if Δ < 0 || rand() < exp(-Δ / temp)
                    current_position = copy(candidate)
                    current_score = candidate_score

                    # Actualizar mejor solución
                    if current_score < best_score
                        best_position = copy(current_position)
                        best_score = current_score
                    end

                    # Registrar la posición actual
                    push!(trajectory, copy(current_position))
                    push!(scores, current_score)
                    break
                end
            end

            # Enfriar la temperatura
            temp *= alpha

            println("Iteración $iter: Mejor costo = $best_score")
            if temp < 1e-6
                break
            end
        end

        return (best_position, best_score, trajectory, scores)
    end

    # Configurar SAMIN
    bounds = [(0.1, 1000), (0.1, 1000), (0.1, 1000)]  # Límites para p[:mU], p[:mW], p[:eP]
    best_position, best_score, trajectory, scores = samin(objective_function, bounds)

    println("Resultado de la optimización: ", best_score)
    println("Mejor conjunto de parámetros encontrado: ", best_position)
end
