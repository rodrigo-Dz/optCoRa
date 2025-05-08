using DelimitedFiles  # Para leer archivos delimitados
using Plots           # Para graficar

# Leer el archivo .txt
filename = "OUT_OptCoRa_FADv1_Fig3_mY_mY_2.txt"  # Cambia esto al nombre de tu archivo
data = readdlm(filename, '\t', skipstart=1)

# Extraer columnas 3 a 5
col3 =  data[:, 3]  # Columna 3
col4 = data[:, 4]  # Columna 4
col5 = data[:, 5]  # Columna 5
# Crear el plot
plot(1:length(col3), col3, label="Columna 3", yscale=:log10)  # Graficar columna 3
plot!(1:length(col4), col4, label="Columna 4", yscale=:log10) # Agregar columna 4 al mismo plot
plot!(1:length(col5), col5, label="Columna 5", yscale=:log10) # Agregar columna 5 al mismo plot

# Personalizar el gráfico
title!("Gráfico de las columnas 3 a 5")
xlabel!("Índice")
ylabel!("Valores")
legend()
