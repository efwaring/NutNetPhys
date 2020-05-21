# Script using least-cost to explore Narea response at nutnet

## load packages

## load functions
source('optimal_vcmax_R/calc_optimal_vcmax.R')
sourceDirectory('optimal_vcmax_R/functions')
source('n_from_gas_exchange/n_from_gas_exchange.R')

## create hypothesis figure
### leaf N and chi by N supply with different differences in N demand
n_supply_trend = calc_optimal_vcmax(beta = seq(100, 300, 10))
n_supply_trend$photo_n = fvcmax25_nrubisco(n_supply_trend$vcmax) + fjmax25_nbioe(n_supply_trend$jmax)
n_supply_trend$photo_n_nobetachange = n_supply_trend$photo_n[10]

hypothesis_data = data.frame(cbind(c(n_supply_trend$beta, n_supply_trend$beta), 
                        c(n_supply_trend$photo_n, n_supply_trend$photo_n_nobetachange),
                        c(rep('no_change', 21), rep('change', 21))))
colnames(hypothesis_data) = c('beta', 'Narea', 'demand')
hypothesis_data$beta = as.numeric(as.character(hypothesis_data$beta))
hypothesis_data$Narea = as.numeric(as.character(hypothesis_data$Narea))

hypothesis_plot = ggplot(data = hypothesis_data, 
                         aes(x = (1/beta), y = Narea, col = demand)) +
  theme(legend.position = 'none', 
        axis.title.y=element_text(size=rel(3), colour = 'black'),
        axis.title.x=element_text(size=rel(3), colour = 'black'),
        axis.text.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.y=element_text(size=rel(2.5), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line()
