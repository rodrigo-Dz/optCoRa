using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using NLsolve
using Plots
using DifferentialEquations, ForwardDiff


mm = include(string("Library/Md_ATFv1.jl"))
include(string("InputFiles/ARGS_ATFv1_Par_Fig1.jl"))

# Extraer el modelo ODE del módulo mm
# Extraer el modelo ODE del módulo mm
odeFB = mm.odeFB

# Vector de estado inicial (ajustar según el problema específico)
u0 = [1.0, 1.0, 1.0, 1.0]  # Valores arbitrarios
params = [p[:g], p[:mY], p[:gY], p[:mU], p[:gU], p[:mW], p[:gW], p[:e0], p[:eP], p[:eM], p[:mUs]]

# Función para encontrar el punto de equilibrio
# Función para encontrar el punto de equilibrio
function equilibrium_point(odeFB, params)
    # Resolver odeFB = 0 para encontrar el equilibrio (estado estacionario)
    # Reemplaza esta parte por la llamada adecuada dependiendo de cómo esté definida la función en odeFB
    # Si f es la función que define el sistema ODE, usa esa.
    function system_eq(x)
        du = similar(x)
        odeFB.f(du, x, params, 0.0)  # Asegúrate de usar la función correcta aquí
        return du
    end
    result = nlsolve(system_eq, u0, autodiff=:forward)
    # Resolver el sistema
    return result.zero
end
# Encontrar el punto de equilibrio

# Calcular la matriz Jacobiana en el punto de equilibrio
def_jacobian(u) = ForwardDiff.jacobian(x -> begin
    du = similar(x)
    odeFB.f(du, x, params, 0.0)
    du
end, u)

# Bifurcation analysis: Variar mY y mW en escala logarítmica
mY_range = 10 .^ range(-3, 3, length=50)  # Rango logarítmico para mY: 0.001 a 1000
mW_range = 10 .^ range(-3, 3, length=50)  # Rango logarítmico para mW: 0.001 a 1000
bifurcation_data = []

for mY in mY_range, mW in mW_range
    params[2] = mY  # Actualizar mY
    params[6] = mW  # Actualizar mW
    # Calcular la matriz Jacobiana en el estado inicial
    equilibrium = equilibrium_point(odeFB, params)
    jac_matrix_at_equilibrium = def_jacobian(equilibrium)
    eigenvalues, eigenvectors = eigen(jac_matrix_at_equilibrium)  # Calcular eigenvalues y eigenvectores

    #l_max = (maximum(real(λ) for λ in eigenvalues))
    #l_min = (minimum(real(λ) for λ in eigenvalues))

    sorted_eigenvalues = sort(eigenvalues, by=abs, rev=true)
    l_slow = sorted_eigenvalues[4]  # Eigenvalue con menor valor absoluto
    l_fast = sorted_eigenvalues[1]  # Eigenvalue con mayor valor absoluto

    #l_slow = eigenvalues[1]
    #l_fast = eigenvalues[4]

    #l_lento = log10(minimum(abs(real(λ)) for λ in eigenvalues))
    #l_fast = log10(maximum(abs(real(λ)) for λ in eigenvalues))

    l_im = maximum(abs(imag(λ)) for λ in eigenvalues; init=0.0)


    idx_max = argmax(abs.(real.(eigenvalues)))
    idx_min = argmin(abs.(real.(eigenvalues)))

    # Obtener los eigenvectores correspondientes
    eigenvector_slow = eigenvectors[:, idx_min]
    eigenvector_fast = eigenvectors[:, idx_max]

    #eigenvector_slow = eigenvectors[:, l_slow]
    #eigenvector_fast = eigenvectors[:, l_fast]


    contribution_Y_slow = (eigenvector_slow[1])  
    contribution_Y_fast = (eigenvector_fast[1])  
    contribution_U_slow = (eigenvector_slow[2])  
    contribution_U_fast = (eigenvector_fast[2])  
    contribution_W_slow = (eigenvector_slow[3])  
    contribution_W_fast = (eigenvector_fast[3])  
    contribution_C_slow = (eigenvector_slow[4])  
    contribution_C_fast = (eigenvector_fast[4]) 

    push!(bifurcation_data, (mY, mW, 
    log10(abs(real(l_slow))), sign(real(l_slow)),
    log10(abs(real(l_fast))), sign(real(l_fast)),
    l_im,
    abs(contribution_Y_slow), abs(contribution_Y_fast),
    abs(contribution_U_slow), abs(contribution_U_fast),
    abs(contribution_W_slow), abs(contribution_W_fast),
    abs(contribution_C_slow), abs(contribution_C_fast)
    ))
end


# Preparar datos para los heatmaps
mY_values = [data[1] for data in bifurcation_data]
mW_values = [data[2] for data in bifurcation_data]
l_slow = [data[3] for data in bifurcation_data]
sign_l_slow = [data[4] for data in bifurcation_data]
l_fast = [data[5] for data in bifurcation_data]
sign_l_fast = [data[6] for data in bifurcation_data]
l_im = [data[7] for data in bifurcation_data]

contribution_Y_slow_values = [data[8] for data in bifurcation_data]
contribution_Y_fast_values = [data[9] for data in bifurcation_data]
contribution_U_slow_values = [data[10] for data in bifurcation_data]
contribution_U_fast_values = [data[11] for data in bifurcation_data]
contribution_W_slow_values = [data[12] for data in bifurcation_data]
contribution_W_fast_values = [data[13] for data in bifurcation_data]
contribution_C_slow_values = [data[14] for data in bifurcation_data]
contribution_C_fast_values = [data[15] for data in bifurcation_data]



# Reshape data for heatmaps
mY_grid = reshape(mY_values, (length(mY_range), length(mW_range)))
mW_grid = reshape(mW_values, (length(mY_range), length(mW_range)))
l_slow_grid = reshape(l_slow, (length(mY_range), length(mW_range)))
sign_l_slow_grid = reshape(sign_l_slow, (length(mY_range), length(mW_range)))
l_fast_grid = reshape(l_fast, (length(mY_range), length(mW_range)))
sign_l_fast_grid = reshape(sign_l_fast, (length(mY_range), length(mW_range)))
l_im_grid = reshape(l_im, (length(mY_range), length(mW_range)))

contribution_Y_slow_grid = reshape(contribution_Y_slow_values, (length(mY_range), length(mW_range)))
contribution_Y_fast_grid = reshape(contribution_Y_fast_values, (length(mY_range), length(mW_range)))
contribution_U_slow_grid = reshape(contribution_U_slow_values, (length(mY_range), length(mW_range)))
contribution_U_fast_grid = reshape(contribution_U_fast_values, (length(mY_range), length(mW_range)))
contribution_W_slow_grid = reshape(contribution_W_slow_values, (length(mY_range), length(mW_range)))
contribution_W_fast_grid = reshape(contribution_W_fast_values, (length(mY_range), length(mW_range)))
contribution_C_slow_grid = reshape(contribution_C_slow_values, (length(mY_range), length(mW_range)))
contribution_C_fast_grid = reshape(contribution_C_fast_values, (length(mY_range), length(mW_range)))



# Gráfica 2: Segundo eigenvalue con mayor valor absoluto
heatmap(mY_range, mW_range, (l_slow_grid), xlabel="mY", ylabel="mW",
        title="Eigenvalie mas lento", xscale=:log10, yscale=:log10,
        color=:heat, colorbar_title="log|λs|")
savefig("./Lambda_lento_abs.png")

# Gráfica 4: Signo del segundo eigenvalue con mayor valor absoluto
heatmap(mY_range, mW_range, sign_l_slow_grid, xlabel="mY", ylabel="mW",
        title="Signo del Eigenvalue mas lento", xscale=:log10, yscale=:log10,
        color=:coolwarm, colorbar_title="Signo de λs")
savefig("./Lambda_lento_sign.png")



heatmap(mY_range, mW_range, (l_fast_grid), xlabel="mY", ylabel="mW",
        title="Eigenvalue mas rapido", xscale=:log10, yscale=:log10,
        color=:heat, colorbar_title="log|λf|")
savefig("./Lambda_fast.png")

heatmap(mY_range, mW_range, sign_l_fast_grid, xlabel="mY", ylabel="mW",
        title="Signo del Eigenvalue mas rapido", xscale=:log10, yscale=:log10,
        color=:coolwarm, colorbar_title="Signo de λf")
savefig("./Lambda_fast_sign.png")

heatmap(mY_range, mW_range, l_im_grid, xlabel="mY", ylabel="mW",
        title="Imaginario", xscale=:log10, yscale=:log10,
        color=:viridis, colorbar_title="|im(λ₂)|")
savefig("./Lambda_im.png")


# Contibucion de Y
heatmap(mY_range, mW_range, contribution_Y_slow_grid, xlabel="mY", ylabel="mW",
        title="Contribución de Y slow", xscale=:log10, yscale=:log10,
        color=:heat, colorbar_title="Contribución de Y")
savefig("./Contribution_Y_slow.png")

heatmap(mY_range, mW_range, contribution_Y_fast_grid, xlabel="mY", ylabel="mW",
        title="Contribución de Y fast", xscale=:log10, yscale=:log10,
        color=:viridis, colorbar_title="Contribución de Y")
savefig("./Contribution_Y_fast.png")

# Contribucion de U
heatmap(mY_range, mW_range, contribution_U_slow_grid, xlabel="mY", ylabel="mW",
        title="Contribución de U slow", xscale=:log10, yscale=:log10,
        color=:heat, colorbar_title="Contribución de U")
savefig("./Contribution_U_slow.png")

heatmap(mY_range, mW_range, contribution_U_fast_grid, xlabel="mY", ylabel="mW",
        title="Contribución de U fast", xscale=:log10, yscale=:log10,
        color=:viridis, colorbar_title="Contribución de U")
savefig("./Contribution_U_fast.png")

# Contribucion de W
heatmap(mY_range, mW_range, contribution_W_slow_grid, xlabel="mY", ylabel="mW",
        title="Contribución de W slow", xscale=:log10, yscale=:log10,
        color=:heat, colorbar_title="Contribución de W")
savefig("./Contribution_W_slow.png")

heatmap(mY_range, mW_range, contribution_W_fast_grid, xlabel="mY", ylabel="mW",
        title="Contribución de W fast", xscale=:log10, yscale=:log10,
        color=:viridis, colorbar_title="Contribución de W")
savefig("./Contribution_W_fast.png")

# Contribucion de C
heatmap(mY_range, mW_range, contribution_C_slow_grid, xlabel="mY", ylabel="mW",
        title="Contribución de C slow", xscale=:log10, yscale=:log10,
        color=:heat, colorbar_title="Contribución de C")
savefig("./Contribution_C_slow.png")

heatmap(mY_range, mW_range, contribution_C_fast_grid, xlabel="mY", ylabel="mW",
        title="Contribución de C fast", xscale=:log10, yscale=:log10,
        color=:viridis, colorbar_title="Contribución de C")
savefig("./Contribution_C_fast.png")
