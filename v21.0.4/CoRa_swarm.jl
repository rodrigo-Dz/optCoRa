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
rtol = 1e-6

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

    # Implementación de PSO
    function pso(objective_function, bounds, swarm_size=5, max_iters=100, ω=0.5, c1=1.5, c2=1.5)
        dim = length(bounds)
        particles = [rand(bounds[i][1]:0.1:bounds[i][2], swarm_size) for i in 1:dim]
        velocities = [zeros(swarm_size) for _ in 1:dim]
        personal_best_positions = deepcopy(particles)
        personal_best_scores = fill(Inf, swarm_size)
        global_best_position = zeros(dim)
        global_best_score = Inf

        for iter in 1:max_iters
            for i in 1:swarm_size
                position = [particles[d][i] for d in 1:dim]
                score = objective_function(position)

                if score < personal_best_scores[i]
                    personal_best_scores[i] = score
                    for d in 1:dim
                        personal_best_positions[d][i] = particles[d][i]
                    end
                end

                if score < global_best_score
                    global_best_score = score
                    global_best_position = position
                end
            end

            # Actualizar posiciones y velocidades
            for i in 1:swarm_size
                for d in 1:dim
                    r1, r2 = rand(), rand()
                    velocities[d][i] = ω * velocities[d][i] +
                                       c1 * r1 * (personal_best_positions[d][i] - particles[d][i]) +
                                       c2 * r2 * (global_best_position[d] - particles[d][i])
                    particles[d][i] += velocities[d][i]

                    # Respetar los límites
                    particles[d][i] = clamp(particles[d][i], bounds[d][1], bounds[d][2])
                end
            end

            println("Iteración ", iter, ": Mejor costo global = ", global_best_score)
        end

        return (global_best_position, global_best_score)
    end

    # Configurar PSO
    bounds = [(0.1, 1000), (0.1, 1000), (0.1, 1000)]  # Límites para p[:mU], p[:mW], p[:eP]
    best_position, best_score = pso(objective_function, bounds)

    println("Resultado de la optimización: ", best_score)
    println("Mejor conjunto de parámetros encontrado: ", best_position)
end
