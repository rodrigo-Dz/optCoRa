using DelimitedFiles
using Plots
gr()  # Selecciona el backend GR

# Leer el archivo "data.txt"
data = readdlm("OUT_ExSSs_FADv1_Fig2_mY_mY.txt", '\t', Float64, header=true)

# Extraer las columnas necesarias
mY = data[1][:, 1]          # Primera columna (mY)
CoRa_mY = data[1][:, end]   # Última columna (CoRa(mY))

# Verificar que los datos se leyeron correctamente
println("mY: ", mY)
println("CoRa(mY): ", CoRa_mY)

# Crear la gráfica
plot(mY, CoRa_mY, seriestype = :line, xscale = :log10,
     xlabel = "mY", ylabel = "CoRa(mY)", ylims=[0,1],
     label = "CoRa(mY)", legend = :topright, grid = true)

# Guardar la gráfica en un archivo (opcional)
savefig("grafica.png")

# Mostrar la gráfica
display(plot)