using Distributed
using DataFrames
using CSV

# Añadir procesos para el procesamiento paralelo
addprocs(4)

# Crear el DataFrame inicial con valores en la columna x
df = DataFrame(x = rand(1:10, 100))  # Genera 100 valores aleatorios entre 1 y 10

# Definir la función que se aplicará a cada valor de x
@everywhere function calcular_cuadrado(fila)
    return fila.(x^2+x)
end

# Ejecutar el procesamiento en paralelo usando pmap
filas = eachrow(df)  # Divide el DataFrame en filas para procesarlas individualmente
resultados = pmap(calcular_cuadrado, filas)

# Agregar los resultados al DataFrame como una nueva columna
df.y = resultados

# Guardar el DataFrame con los resultados en un archivo CSV
CSV.write("resultados_cuadrados.csv", df)

# Mostrar el DataFrame resultante
println(df)
