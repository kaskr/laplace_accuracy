# Poisson with aggregated latent intensities

Model:

u1 ~ N(mu1, sd1)
u2 ~ N(mu2, sd2)

N|u ~ Pois(exp(u1)+exp(u2))
