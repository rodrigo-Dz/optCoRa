using ForwardDiff, NLsolve, LinearAlgebra

# Definir el sistema de ecuaciones (compatible con ForwardDiff)
function sistema!(du, u, p)
    g, mY, gY, mU, gU, mW, gW, e0, eP, eM = p
    Y, U, W, C = u

    du[1] = (mY * W) - ((g + gY) * Y)  # dY/dt
    du[2] = (mU * Y) - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)  # dU/dt
    du[3] = mW - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)  # dW/dt
    du[4] = -((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)  # dC/dt

    return du
end

# Función para encontrar el estado estacionario
function estado_estacionario(p)
    # Condiciones iniciales para Y, U, W, C
    u0 = [1.0, 1.0, 1.0, 1.0]

    # Resolver el sistema no lineal
    result = nlsolve((du, u) -> sistema!(du, u, p), u0)
    return result.zero  # Devuelve el estado estacionario
end

# Calcular la matriz jacobiana en el estado estacionario
function jacobiana_estado_estacionario(p)
    u_estacionario = estado_estacionario(p)
    ForwardDiff.jacobian(u -> sistema!(similar(u), u, p), u_estacionario)
end

# Parámetros del sistema
p = (
    0.0001,  # g
    0.01,    # mY  ---
    0.1,     # gY
    0.1,     # mU 
    0.0001,  # gU
    100,    # mW ---
    0.0001,  # gW
    0.0001,  # e0
    0.1,     # eP
    0.5      # eM
)

# Encontrar el estado estacionario
u_estacionario = estado_estacionario(p)
println("Estado estacionario:")
println(u_estacionario)

# Calcular la matriz jacobiana en el estado estacionario
J = jacobiana_estado_estacionario(p)
println("Matriz Jacobiana en el estado estacionario:")
println(J)

# Calcular eigenvalores y eigenvectores
eigen_result = eigen(J)
println("Eigenvalores:")
println(eigen_result.values)
println("Eigenvectores:")
println(eigen_result.vectors)