# Helpers Functions for plotting Columbus results

## Layout Plotting Functions
plot_plate <-
  function(table,
           property,
           discrete_property = T,
           legend,
           title) {
    
    layout_plot <- ggplot(table,
                          aes(
                            x = column,
                            y = row,
                            fill = {{property}},
                          )) + 
      geom_tile(color = "grey30") +
      scale_y_reverse(breaks = 1:16, labels = LETTERS[1:16]) +
      scale_x_continuous(breaks = 0:24) +
      coord_equal() +
      facet_wrap(vars(screen))+
      xlab("Column") +
      ylab("Row") +
      ggtitle(title)
    
    if (discrete_property) {
      layout_plot + scale_fill_tableau(palette = "Tableau 20",
                                       name = legend)
    }
    else{
      layout_plot + scale_fill_viridis_c(legend)
    }
  }

## Density Plotting Functions
plot_density <-
  function(table,
           property,
           property_legend,
           tint,
           tint_legend,
           title) {
    
    density_plot <- ggplot(table,
                           aes(
                             x = {{property}},
                             y = ..density..,
                             color = {{tint}},
                             fill = {{tint}}))
    
    density_plot + geom_density(alpha = 0.3) +
      scale_x_continuous(property_legend, 
                         trans = "log2") +
      scale_color_few(name = tint_legend) +
      scale_fill_few(name = tint_legend) +
      facet_grid(vars(marker),
                 vars(screen),
                 scales = "free_x") +
      ggtitle(title)
    
  }

### Crossbar plotting functions

plot_crossbars<- function(table,
                           property,
                           property_legend,
                           title) {
  

    crossbar_plot <- ggplot(table, 
                            aes(x = cell_line,
                                y = {{property}}))
    
    crossbar_plot + 
      geom_point(alpha = 0.5,
                 position = position_jitter(width = 0.2)) +
      stat_summary(fun.data = 'mean_sdl',
                   fun.args = list(mult = 1),
                   geom = "crossbar",
                   width = 0.75,) +
      facet_wrap(vars(marker),
                 scales = "free_y") +
      expand_limits(y = 0) +
      ylab(property_legend) +
      xlab("Cell Line") +
      ggtitle(title)
}