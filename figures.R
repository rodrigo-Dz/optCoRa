library(dbscan)
library(plotly)
library(dplyr)
library(readr)
library(scales)
library(geometry)
library(viridisLite)  
library(htmlwidgets)
library(tidyr)
library(dplyr)
library(viridis)
library(patchwork)
library(ggplot2)

setwd("~/Desktop/optimizeCoRa_29sep/optimizeCoRa")

### ---- Explore ----
 # Read file
curves <- read_tsv("/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/OUT_ExplCoRa_FADv1_p01_Fig1_mY_mY.txt", col_names = TRUE)
model = "FADv1_p01"
exp = "Fig1"
 # Different colors in case of other behaviors in the simulations
cont_colors <- viridis(100)[as.numeric(cut(curves$robustness, breaks = 100))]
final_colors <- ifelse(curves$oscilations != 0, "purple",
                       ifelse(curves$negative_sol != 0, "red",
                              ifelse(curves$other_errors != 0, "orange",
                                     ifelse(curves$negative_CoRa!= 0, "blue",
                                            ifelse(curves$Greater_CoRa != 0, "pink",
                                                   ifelse(curves$not_same != 0, "green",
                                     cont_colors))))))
curves$color <- final_colors
curves$is_special <- with(curves,
                          oscilations != 0 |
                            negative_sol != 0 |
                            other_errors != 0 |
                            negative_CoRa != 0 |
                            Greater_CoRa != 0 |
                            not_same != 0
)


### --- Heatmap Robustness 2D

p1 <- "mW"
p2 <- "eP"

p_heatmap <- plot_ly(
  data = curves,
  x = curves[[p1]],
  y = curves[[p2]],
  z = ~robustness,
  type = "contour",
  colorscale = "Viridis",
  hoverinfo = "text",
  text = ~paste0(
    p1, "=", signif(.data[[p1]], 4),
    "<br>", p2, "=", signif(.data[[p2]], 4),
    "<br>robustness=", signif(robustness, 4)
  ),
  zmin = 0,
  zmax = 1,
  showscale = FALSE,
  contours = list(showlabels = TRUE,
                  labelfont = list(size = 15, color = "white"))
) %>%
  add_markers(
    data = curves %>% filter(is_special),
    x = ~ .data[[p1]],
    y = ~ .data[[p2]],
    marker = list(color = ~color, size = 7, opacity = 0.95),
    inherit = FALSE,
    showlegend = FALSE
  ) %>%
  layout(
    xaxis = list(title = list(text = p1, font = list(size = 18)),
                 type = "log", tickfont = list(size = 14)),
    yaxis = list(title = list(text = p2, font = list(size = 18)),
                 type = "log", tickfont = list(size = 14))
  )

p_heatmap



### Heatmap 2D - Steady State
p_heatmap <- plot_ly(
  data = curves,
  x = curves[[p1]],
  y = curves[[p2]],
  z = ~log10(steady_state),
  type = "contour",
  hoverinfo = "text",
  text = ~paste0(
    p1, "=", signif(.data[[p1]], 4),
    "<br>", p2, "=", signif(.data[[p2]], 4),
    "<br>robustness=", signif(robustness, 4)
  ),
  colorscale = "Viridis",
  showscale = FALSE,
  contours = list(showlabels = TRUE,
                  labelfont = list(size = 15, color = "white"))
) %>%
  add_markers(
    data = curves %>% filter(is_special),
    x = ~ .data[[p1]],
    y = ~ .data[[p2]],
    marker = list(color = ~color, size = 7, opacity = 0.95),
    inherit = FALSE,
    showlegend = FALSE
  ) %>%
  layout(
    xaxis = list(title = list(text = p1, font = list(size = 18)),
                 type = "log", tickfont = list(size = 14)),
    yaxis = list(title = list(text = p2, font = list(size = 18)),
                 type = "log", tickfont = list(size = 14))
  )

p_heatmap



###### Comparaciones entre distintos parametros

n_par = 3   # numero de parametros en la exploracion
curves[1:n_par] <- log10(curves[1:n_par])
param_combinations <- combn(names(curves)[1:n_par], 2, simplify = FALSE)


for (i in seq_along(param_combinations)) {
  p <- ggplot(curves) +
    geom_point(aes(
      x = .data[[param_combinations[[i]][1]]],
      y = .data[[param_combinations[[i]][2]]],
      color = robustness
    ), size = 0.4) +  
    scale_color_viridis_c(name = "robustness") +
    xlim(-3, 3) +
    ylim(-3, 3) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste(param_combinations[[i]][1], "vs", param_combinations[[i]][2]),
      x = param_combinations[[i]][1],
      y = param_combinations[[i]][2]
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
  file_name <- paste0(
    "./Output/OUT_ExplCoRa_" , model, "_", exp, "_",
    param_combinations[[i]][1], "_vs_",
    param_combinations[[i]][2], ".png"
  )
  ggsave(
    filename = file_name,
    plot = p,
    width = 6,
    height = 5,
    dpi = 300,
    bg = "white"
  )
}


### HEATMAP 3D ROBUSTNESS 
p1 <- "mW"
p2 <- "mU"
p3 <- "eP"
p <- plot_ly(curves, 
             x = curves[[p1]], 
             y = curves[[p2]], 
             z = curves[[p3]], 
             marker = list(
               size = 10, 
               color = ~robustness,
               colorscale = "Viridis",
               colorbar = list(title = "Robustness"),
               showscale = TRUE
             ),
             hoverinfo = "text",
             text = ~paste0(
               p1, "=", signif(.data[[p1]], 4),
               "<br>", p2, "=", signif(.data[[p2]], 4),
               "<br>", p3, "=", signif(.data[[p3]], 4),
               "<br>robustness=", signif(robustness, 4)
             ),
             type = "scatter3d", 
             mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = p1, type = "log", range = c(log10(0.001), log10(1000))),
    yaxis = list(title = p2, type = "log", range = c(log10(0.001), log10(1000))),
    zaxis = list(title = p3, type = "log", range = c(log10(0.001), log10(1000)))
  )) %>%
  add_markers(
    data = curves %>% filter(is_special),
    x = ~ .data[[p1]],
    y = ~ .data[[p2]],
    z = ~ .data[[p3]],
    marker = list(color = ~color, size = 7, opacity = 0.95),
    inherit = FALSE,
    showlegend = FALSE
  )
p



### HEATMAP 3D SS 
p1 <- "mW"
p2 <- "mU"
p3 <- "eP"
p <- plot_ly(curves, 
             x = curves[[p1]], 
             y = curves[[p2]], 
             z = curves[[p3]], 
             marker = list(
               size = 10, 
               color = ~log10(steady_state),
               colorscale = "Viridis",
               colorbar = list(title = "Robustness"),
               showscale = TRUE
             ),
             hoverinfo = "text",
             text = ~paste0(
               p1, "=", signif(.data[[p1]], 4),
               "<br>", p2, "=", signif(.data[[p2]], 4),
               "<br>", p3, "=", signif(.data[[p3]], 4),
               "<br>robustness=", signif(robustness, 4)
             ),
             type = "scatter3d", 
             mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = p1, type = "log", range = c(log10(0.001), log10(1000))),
    yaxis = list(title = p2, type = "log", range = c(log10(0.001), log10(1000))),
    zaxis = list(title = p3, type = "log", range = c(log10(0.001), log10(1000)))
  )) %>%
  add_markers(
    data = curves %>% filter(is_special),
    x = ~ .data[[p1]],
    y = ~ .data[[p2]],
    z = ~ .data[[p3]],
    marker = list(color = ~color, size = 7, opacity = 0.95),
    inherit = FALSE,
    showlegend = FALSE
  )
p



# ---- Optimize ----

optimization_history <- read_tsv("./Output/OUT_OptCoRa_FADv1_p01_Fig1_mY_mY.txt", col_names = TRUE)
n_par = 3
p1 = "mU"
p2 = "mW"
p3 = "eP"

hist(optimization_history$robustness)
hist(optimization_history$min_CoRa)
hist(optimization_history$mU, xlim=c(-3,3))


# Random Walk 3D
p <- plot_ly(optimization_history, 
             x = optimization_history[[p1]], 
             y = optimization_history[[p2]],
             z = optimization_history[[p3]],
             marker = list(
               size = 10, 
               color = ~robustness,
               colorscale = "Viridis",
               colorbar = list(title = "Robustness"),
               showscale = TRUE
             ),
             type = "scatter3d", 
             mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = p1, type = "log", range = c(log10(0.001), log10(1000))),
    yaxis = list(title = p2, type = "log", range = c(log10(0.001), log10(1000))),
    zaxis = list(title = p3, type = "log", range = c(log10(0.001), log10(1000)))
  ))
p

# Distribucion de SS
best <- optimization_history[optimization_history$robustness > 0.8,]
last <- optimization_history[8000:10000, ]
ggplot(last, aes(x = log10(steady_state))) +
  geom_histogram(
    bins = 40,
    fill = "#215EA1",   # azul profundo, puedes cambiarlo
    color = "white",
    alpha = 0.9
  ) +
  labs(
    x = expression(log[10](steady_state)),
    y = "Count",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )





optimization_history[2:(1+n_par)] <- log10(optimization_history[2:(1+n_par)])
param_combinations <- combn(names(optimization_history)[2:(1+n_par)], 2, simplify = FALSE)


for (i in seq_along(param_combinations)) {
  p <- ggplot(optimization_history) +
    geom_point(aes(
      x = .data[[param_combinations[[i]][1]]],
      y = .data[[param_combinations[[i]][2]]],
      color = robustness
    ), size = 0.4) +  
    scale_color_viridis_c(name = "robustness") +
    xlim(-3, 3) +
    ylim(-3, 3) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste(param_combinations[[i]][1], "vs", param_combinations[[i]][2]),
      x = param_combinations[[i]][1],
      y = param_combinations[[i]][2]
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
  file_name <- paste0(
    "./Output/OUT_OptCoRa_" , model, "_", exp, "_robustness_",
    param_combinations[[i]][1], "_vs_",
    param_combinations[[i]][2], ".png"
  )
  ggsave(
    filename = file_name,
    plot = p,
    width = 6,
    height = 5,
    dpi = 300,
    bg = "white"
  )
}



##### Seady state


for (i in seq_along(param_combinations)) {
  p <- ggplot(optimization_history) +
    geom_point(aes(
      x = .data[[param_combinations[[i]][1]]],
      y = .data[[param_combinations[[i]][2]]],
      color = log10(steady_state)
    ), size = 0.4) +  
    scale_color_viridis_c(name = "log10(SS)") +
    xlim(-3, 3) +
    ylim(-3, 3) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste(param_combinations[[i]][1], "vs", param_combinations[[i]][2]),
      x = param_combinations[[i]][1],
      y = param_combinations[[i]][2]
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
  file_name <- paste0(
    "./Output/OUT_OptCoRa_" , model, "_", exp, "_SS_",
    param_combinations[[i]][1], "_vs_",
    param_combinations[[i]][2], ".png"
  )
  ggsave(
    filename = file_name,
    plot = p,
    width = 6,
    height = 5,
    dpi = 300,
    bg = "white"
  )
}



######  Density

optimization_history <- optimization_history %>% filter(robustness>0.8)
p1 = "mU"
p_hist <- ggplot(optimization_history, aes(x = optimization_history[[p1]])) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),  # convierte a proporción (suma 1)
    bins = 100,
    fill = "#215EA1",
    alpha = 1,
    color = "white"
  ) +
  scale_x_log10(
    limits = c(1e-3, 1e3),
    breaks = 10^seq(-3, 3, by = 1),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  labs(x = p1, y = "Proportion") +
  theme_minimal() +
  theme(
    panel.grid = element_line(color = "#E5ECF6"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

p_hist
file_name_hist <- paste0("./Output/OUT_OptCoRa_" , model, "_", exp,"_",p1,"_hist_log.png")
ggsave(
  filename = file_name_hist,
  plot = p_hist,
  width = 8,
  height = 5,
  dpi = 300
)





########### 
library(tidyverse)

# Figura 1: ATFv1 exp04
optimization_history <- read_tsv("./Output/OUT_OptCoRa_ATFv1_p01_Fig4_mY_mY.txt", col_names = TRUE)
model = "ATFv1_p01"
exp = "Fig4"


p_density <- optimization_history %>%
  pivot_longer(c(mU, mW, eP, gU, gW, e0, eM),
               names_to = "Parameter", values_to = "Value") %>%
  mutate(
    Value = log10(Value),
    Parameter = factor(Parameter, levels = c("mU", "mW", "eP", "gU", "gW", "e0", "eM"))
  ) %>%
  ggplot(aes(x = Value, fill = Parameter)) +
  geom_density(
    aes(y = after_stat(density)),
    alpha = 0.6,
    color = NA,
    adjust = 1
  ) +
  xlim(-3, 3) +
  scale_fill_manual(
    values = c(
      "mU" = "#215EA1",
      "mW" = "#22BAB9",
      "eP" = "#94E2EA",
      "gU" = "grey20",
      "gW" = "grey40",
      "e0" = "grey60",
      "eM" = "grey80"
    ),
    name = ""
  ) +
  labs(x = expression(log[10]("Parameter value")), y = "Density") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid = element_line(color = "#E5ECF6"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

p_density
file_name <- paste0("./Output/OUT_OptCoRa_" , model, "_", exp,"_density.png")
ggsave(
  filename = file_name,
  plot = p_density,
  width = 8,
  height = 5,
  dpi = 300
)


p_density <- optimization_history %>%
  pivot_longer(c(mU, mW, eP, gU, gW, e0, eM),
               names_to = "Parameter", values_to = "Value") %>%
  mutate(
    Value = log10(Value),
    Parameter = factor(Parameter, levels = c("mU", "mW", "eP", "gU", "gW", "e0", "eM"))
  ) %>%
  ggplot(aes(x = Value, fill = Parameter)) +
  geom_histogram(
    position = "identity",
    alpha = 0.5,
    bins = 100, size = 0
  ) +
  xlim(-3, 3) +
  scale_fill_manual(
    values = c(
      "mU" = "#215EA1",
      "mW" = "#22BAB9",
      "eP" = "#94E2EA",
      "gU" = "grey20",
      "gW" = "grey40",
      "e0" = "grey60",
      "eM" = "grey80"
    ),
    name = ""
  ) +
  labs(x = expression(log[10]("Parameter value")), y = "Density") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid = element_line(color = "#E5ECF6"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

p_density
file_name <- paste0("./Output/OUT_OptCoRa_" , model, "_", exp,"_hist.png")
ggsave(
  filename = file_name,
  plot = p_density,
  width = 8,
  height = 5,
  dpi = 300
)

p_density_violin_horizontal <- optimization_history %>%
  pivot_longer(c(mU, mW, eP, gU, gW, e0, eM),
               names_to = "Parameter", values_to = "Value") %>%
  mutate(
    Value = log10(Value),
    Parameter = factor(Parameter, levels = c("mU", "mW", "eP", "gU", "gW", "e0", "eM"))
  ) %>%
  ggplot(aes(x = Value, y = Parameter, fill = Parameter)) +
  geom_violin(
    alpha = 0.7,
    trim = TRUE,
    scale = "width",
    width = 0.8,
    orientation = "y",
    color = "black",
    size = 0.2
  ) +
  geom_boxplot(
    width = 0.15,
    alpha = 0.3,
    outlier.shape = NA,
    fill = "white",
    orientation = "y",
    color = "black",  # Mantiene el contorno negro del boxplot
    size = 0.2  # Grosor del contorno del boxplot
  ) +
  xlim(-3, 3) +
  scale_fill_manual(
    values = c(
      "mU" = "#215EA1",
      "mW" = "#22BAB9",
      "eP" = "#94E2EA",
      "gU" = "grey20",
      "gW" = "grey40",
      "e0" = "grey60",
      "eM" = "grey80"
    ),
    guide = "none",
    name = ""
  ) +
  labs(
    x = "Value (log scale)",
    y = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_line(color = "#E5ECF6"),
    panel.grid.minor.x = element_line(color = "#E5ECF6", linetype = "dotted"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

p_violin_horizontal

file_name <- paste0("./Output/OUT_OptCoRa_" , model, "_", exp,"_violin.png")
ggsave(
  filename = file_name,
  plot = p_violin_horizontal,
  width = 8,
  height = 5,
  dpi = 300
)





#### Iteraciones 
optimization_history <- read_tsv("./Output/OUT_OptCoRa_ATFv1_p01_Fig4_mY_mY.txt", col_names = TRUE)
p_optimized = 7
optimization_history[1:p_optimized+1] <- log10(optimization_history[1:p_optimized+1])
#optim <- optimization_history[1:p_optimized+1]

library(ggplot2)
library(tidyr)
library(patchwork)

historial_df <- optimization_history[1:500, ]

#c(mU, mW, eP)
#c(mU, mW. eP, gU, gW, e0, eM)
#c(mA, mB, mU, e0, eP, bA, bI, kD)
p1 <- ggplot(historial_df %>% 
               pivot_longer(c(mU, mW, eP, gU, gW, e0, eM)), 
             aes(x = Iteration, y = value, color = name)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c(
      "mU" = "#215EA1",
      "mW" = "#22BAB9",
      "eP" = "#94E2EA",
      "gU" = "grey20",
      "gW" = "grey40",
      "e0" = "grey60",
      "eM" = "grey80"
    ),
    name = "Parameter"
  ) +
  labs(y = "Parameter") +
  theme_minimal()

p1
p2 <- ggplot(historial_df, aes(x = Iteration)) +
  geom_line(aes(y = robustness, color = "Robustness"), linewidth = 1) +
  geom_line(aes(y = min_CoRA, color = "min(CoRa)"), linewidth = 1) +
  labs (y = "f(x)") +
  scale_color_manual(values = c("Robustness" = "#22BAB9", "min(CoRa)" = "#215EA1"), name = "Metrics") +
  theme_minimal()+
  ylim(0, 1)

p2

iteraciones_seleccionadas <- c(1, 50, 100, 150, 300)

sub <- historial_df %>% filter(Iteration %in% iteraciones_seleccionadas)

library(ggplot2)
library(dplyr)
library(tidyr)

plot_curves2D <- function(data, titulo, iteraciones) {
  data_processed <- data %>%
    mutate(curve_id = row_number()) %>%
    pivot_longer(
      cols = 15:42,
      names_to = "curve",
      values_to = "CoRa"
    ) %>%
    mutate(curve = as.numeric(curve)) %>%
    group_by(Iteration) %>%
    mutate(unique_par = first(robustness)) %>%
    ungroup() %>%
    # Asegurar que las iteraciones estén en el orden correcto
    mutate(Iteration_factor = factor(Iteration, levels = iteraciones))
  
  n_iter <- length(iteraciones)
  colores_gris_azul <- colorRampPalette(c( "lightgray", "darkgray","#94E2EA", "#22BAB9","#215EA1"))(n_iter)
  names(colores_gris_azul) <- as.character(iteraciones)
  
  ggplot(data_processed, aes(x = log10(curve), y = CoRa, group = Iteration_factor, 
                             color = Iteration_factor)) +
    geom_line(alpha = 0.8, linewidth = 1) +
    scale_color_manual(
      values = colores_gris_azul,
      name = "Iteration"
    ) +
    labs(
      x = "mY",
      y = "CoRa",
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "#E5ECF6"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.key.height = unit(0.5, "cm"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    ) +
    ylim(0, 1)
}

fig1 <- plot_curves2D(sub, "F1", iteraciones_seleccionadas)


ff <- p1/p2/fig1

file_name <- paste0("./Output/OUT_OptCoRa_" , model, "_", exp,"_randomwalk.png")
ggsave(
  filename = file_name,
  plot = ff,
  width = 10,
  height = 5,
  dpi = 300
)



