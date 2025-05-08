using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using Plots
using NLsolve

# Define el sistema de ODEs
function ode_system!(du, u, p, t)
    Y, U, W, C = u
    g, mY, gY, mU, gU, mW, gW, e0, eP, eM = p

    du[1] = (mY * W) - ((g + gY) * Y)  # dY/dt
    du[2] = (mU * Y) - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)  # dU/dt
    du[3] = mW - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)  # dW/dt
    du[4] = -((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)  # dC/dt
end

# Parámetros
params = [
        0.0001,  # g
        1.0,   # mY (se variará)
        1,   # gY
        1, # mU (se variará)
        0.0001,  # gU
        1,  # mW
        0.0001, # gW
        0.0001, # e0
        1, # eP
        0.5    # eM
]

# Condiciones iniciales
u0 = [0.0, 0.0, 0.0, 0.0]  # Y, U, W, C

# Rango de tiempo
tspan = (0.0, 1000.0)  # Tiempo largo para alcanzar el estado estacionario

# Encontrar puntos de equilibrio
function find_equilibrium(p)
    function equilibrium_condition(u)
        du = similar(u)
        ode_system!(du, u, p, 0.0)
        return du
    end
    result = nlsolve(equilibrium_condition, u0, autodiff=:forward)
    return result.zero
end

# Calcular el Jacobiano en un punto de equilibrio
function compute_jacobian(u_eq, p)
    # Función que retorna las derivadas como un arreglo
    function system_derivatives(u)
        du = similar(u)
        ode_system!(du, u, p, 0.0)
        return du
    end
    ForwardDiff.jacobian(system_derivatives, u_eq)
end

# Bifurcation analysis: Variar mY y mU en escala logarítmica
mY_range = 10 .^ range(-3, 3, length=50)  # Rango logarítmico para mY: 0.01 a 100
mU_range = 10 .^ range(-3, 3, length=50)  # Rango logarítmico para mU: 0.01 a 100
bifurcation_data = []

for mY in mY_range, mU in mU_range
    params[2] = mY  # Actualizar mY
    params[6] = mU  # Actualizar mU
    u_eq = find_equilibrium(params)
    J = compute_jacobian(u_eq, params)
    eigenvalues, eigenvectors = eigen(J)  # Calcular eigenvalues y eigenvectores

    # Encontrar el eigenvalue más grande y el más pequeño
    max_real_positive = (maximum(real(λ) for λ in eigenvalues))
    min_real = (minimum(real(λ) for λ in eigenvalues))

    max_imaginary = maximum(abs(imag(λ)) for λ in eigenvalues; init=0.0)
    

    # Encontrar los índices de los eigenvalues más grande y más pequeño
    idx_max = argmax(real.(eigenvalues))
    idx_min = argmin(real.(eigenvalues))

    # Obtener los eigenvectores correspondientes
    eigenvector_max = eigenvectors[:, idx_max]
    eigenvector_min = eigenvectors[:, idx_min]

    # Extraer la contribución de Y en los eigenvectores
    contribution_Y_max = abs(eigenvector_max[1])  # Primer componente (Y)
    contribution_Y_min = abs(eigenvector_min[1])  # Primer componente (Y)
    contribution_U_max = abs(eigenvector_max[2])  # Primer componente (Y)
    contribution_U_min = abs(eigenvector_min[2])  # Primer componente (Y)contribution_Y_max = abs(eigenvector_max[1])  # Primer componente (Y)
    contribution_W_max = abs(eigenvector_max[3])  # Primer componente (Y)contribution_Y_max = abs(eigenvector_max[1])  # Primer componente (Y)
    contribution_W_min = abs(eigenvector_min[3])  # Primer componente (Y)
    contribution_C_max = abs(eigenvector_max[4])  # Primer componente (Y)contribution_Y_max = abs(eigenvector_max[1])  # Primer componente (Y)
    contribution_C_min = abs(eigenvector_min[4]) 
    push!(bifurcation_data, (mY, mU, max_real_positive, min_real, contribution_Y_max, contribution_Y_min,  contribution_U_max, contribution_U_min, contribution_W_max, contribution_W_min, contribution_C_max, contribution_C_min, max_imaginary, sign(real(max_real_positive)), sign(real(min_real))))
end

# Preparar datos para los heatmaps
mY_values = [data[1] for data in bifurcation_data]
mU_values = [data[2] for data in bifurcation_data]
max_real_positive_values = [data[3] for data in bifurcation_data]
min_real_values = [data[4] for data in bifurcation_data]
contribution_Y_max_values = [data[5] for data in bifurcation_data]
contribution_Y_min_values = [data[6] for data in bifurcation_data]
contribution_U_max_values = [data[7] for data in bifurcation_data]
contribution_U_min_values = [data[8] for data in bifurcation_data]
contribution_W_max_values = [data[9] for data in bifurcation_data]
contribution_W_min_values = [data[10] for data in bifurcation_data]
contribution_C_max_values = [data[11] for data in bifurcation_data]
contribution_C_min_values = [data[12] for data in bifurcation_data]
max_imaginary_values = [data[13] for data in bifurcation_data]
sign_lambda1 = [data[14] for data in bifurcation_data]
sign_lambda2 = [data[15] for data in bifurcation_data]

# Reshape data for heatmaps
mY_grid = reshape(mY_values, (length(mY_range), length(mU_range)))
mU_grid = reshape(mU_values, (length(mY_range), length(mU_range)))
max_real_positive_grid = reshape(max_real_positive_values, (length(mY_range), length(mU_range)))
min_real_grid = reshape(min_real_values, (length(mY_range), length(mU_range)))
max_imaginary_grid = reshape(max_imaginary_values, (length(mY_range), length(mU_range)))

contribution_Y_max_grid = reshape(contribution_Y_max_values, (length(mY_range), length(mU_range)))
contribution_Y_min_grid = reshape(contribution_Y_min_values, (length(mY_range), length(mU_range)))
contribution_U_max_grid = reshape(contribution_U_max_values, (length(mY_range), length(mU_range)))
contribution_U_min_grid = reshape(contribution_U_min_values, (length(mY_range), length(mU_range)))
contribution_W_max_grid = reshape(contribution_W_max_values, (length(mY_range), length(mU_range)))
contribution_W_min_grid = reshape(contribution_W_min_values, (length(mY_range), length(mU_range)))
contribution_C_max_grid = reshape(contribution_C_max_values, (length(mY_range), length(mU_range)))
contribution_C_min_grid = reshape(contribution_C_min_values, (length(mY_range), length(mU_range)))
sign_lambda1_grid = reshape(sign_lambda1, (length(mY_range), length(mU_range)))
sign_lambda2_grid = reshape(sign_lambda2, (length(mY_range), length(mU_range)))
# Gráfica 1: Contribución de Y en el eigenvector del eigenvalue más grande
heatmap(mY_range, mU_range, contribution_Y_max_grid, xlabel="mY", ylabel="mW",
        title="Contribución de Y (Eigenvector del λ más grande)", xscale=:log10, yscale=:log10,
        color=:heat, colorbar_title="Contribución de Y")
savefig("Contribution_Y_max.png")

# Gráfica 2: Contribución de Y en el eigenvector del eigenvalue más pequeño
heatmap(mY_range, mU_range, contribution_Y_min_grid, xlabel="mY", ylabel="mW",
        title="Contribución de Y (Eigenvector del λ más pequeño)", xscale=:log10, yscale=:log10,
        color=:viridis, colorbar_title="Contribución de Y")
savefig("Contribution_Y_min.png")

# Gráfica 1: Contribución de Y en el eigenvector del eigenvalue más grande
heatmap(mY_range, mU_range, contribution_U_max_grid, xlabel="mY", ylabel="mW",
        title="Contribución de U (Eigenvector del λ más grande)", xscale=:log10, yscale=:log10,
        color=:heat, colorbar_title="Contribución de U")
savefig("Contribution_U_max.png")

# Gráfica 2: Contribución de Y en el eigenvector del eigenvalue más pequeño
heatmap(mY_range, mU_range, contribution_U_min_grid, xlabel="mY", ylabel="mW",
        title="Contribución de U (Eigenvector del λ más pequeño)", xscale=:log10, yscale=:log10,
        color=:viridis, colorbar_title="Contribución de U")
savefig("Contribution_U_min.png")

# Gráfica 1: Contribución de Y en el eigenvector del eigenvalue más grande
heatmap(mY_range, mU_range, contribution_W_max_grid, xlabel="mY", ylabel="mW",
        title="Contribución de W (Eigenvector del λ más grande)", xscale=:log10, yscale=:log10,
        color=:heat, colorbar_title="Contribución de W")
savefig("Contribution_W_max.png")

# Gráfica 2: Contribución de Y en el eigenvector del eigenvalue más pequeño
heatmap(mY_range, mU_range, contribution_W_min_grid, xlabel="mY", ylabel="mW",
        title="Contribución de W (Eigenvector del λ más pequeño)", xscale=:log10, yscale=:log10,
        color=:viridis, colorbar_title="Contribución de W")
savefig("Contribution_W_min.png")

heatmap(mY_range, mU_range, contribution_C_max_grid, xlabel="mY", ylabel="mW",
        title="Contribución de C (Eigenvector del λ más grande)", xscale=:log10, yscale=:log10,
        color=:heat, colorbar_title="Contribución de C")
savefig("Contribution_C_max.png")

# Gráfica 2: Contribución de Y en el eigenvector del eigenvalue más pequeño
heatmap(mY_range, mU_range, contribution_C_min_grid, xlabel="mY", ylabel="mW",
        title="Contribución de C (Eigenvector del λ más pequeño)", xscale=:log10, yscale=:log10,
        color=:viridis, colorbar_title="Contribución de C")
savefig("Contribution_C_min.png")





heatmap(mY_range, mU_range, max_real_positive_grid, xlabel="mY", ylabel="mW",
        title="Parte Real Maxima de Eigenvalues", xscale=:log10, yscale=:log10,
        color=:heat, colorbar_title="Máx Re(λ) > 0")
        savefig("Bif_max_real.png")

heatmap(mY_range, mU_range, min_real_grid, xlabel="mY", ylabel="mW",
    title="Parte Real Minima de Eigenvalues", xscale=:log10, yscale=:log10,
    color=:viridis, colorbar_title="Min Re(λ) > 0")
    savefig("Bif_min_real.png")
# Gráfica 2: Parte imaginaria de los eigenvalues
heatmap(mY_range, mU_range, max_imaginary_grid, xlabel="mY", ylabel="mW",
        title="Parte Imaginaria de Eigenvalues", xscale=:log10, yscale=:log10,
        color=:viridis, colorbar_title="Máx |Im(λ)|")
        savefig("Bif_imag.png")

        # Gráfica 3: Signo del eigenvalue con mayor valor absoluto
heatmap(mY_range, mU_range, sign_lambda1_grid, xlabel="mY", ylabel="mW",
title="Signo del Eigenvalue con Mayor Valor Absoluto", xscale=:log10, yscale=:log10,
color=:coolwarm, colorbar_title="Signo de λ₁")
savefig("Lambda1_sign.png")

# Gráfica 4: Signo del segundo eigenvalue con mayor valor absoluto
heatmap(mY_range, mU_range, sign_lambda2_grid, xlabel="mY", ylabel="mW",
title="Signo del Segundo Eigenvalue con Mayor Valor Absoluto", xscale=:log10, yscale=:log10,
color=:coolwarm, colorbar_title="Signo de λ₂")
savefig("Lambda2_sign.png")