# how to convert d13C to ci/ca (Cernusak et al., 2013)

air = -8
plant = seq(-28, -15, 1)
delta = (air - plant) / (1 + plant * 0.001)
a = 4.4
b = 27
chi = (delta - a) / (b - a)
chi
