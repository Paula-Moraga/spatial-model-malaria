# Installing INLA

install.packages("INLA", repos = "https://inla.r-inla-download.org/R/stable", dep = TRUE)

# Loading R packages

library(INLA)
library(leaflet)
library(viridis)
library(ggplot2)
library(cowplot)

# Reading data

d <- read.csv("https://raw.githubusercontent.com/Paula-Moraga/spatial-model-malaria/master/d.csv")
dp <- read.csv("https://raw.githubusercontent.com/Paula-Moraga/spatial-model-malaria/master/dp.csv")

# Fitting model

coo <- cbind(d$longitude, d$latitude)
mesh <- inla.mesh.2d(loc = coo, max.edge = c(0.1, 5), cutoff = 0.01)
spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
indexs <- inla.spde.make.index("s", spde$n.spde)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

coop <- cbind(dp$longitude,dp$latitude)
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

stk.e <- inla.stack(tag = "est", data = list(y = d$positive,
                    numtrials = d$examined), A = list(1, A),
                    effects = list(data.frame(b0 = 1, alt = d$alt,
                    temp = d$temp, prec = d$prec, hum = d$hum,
                    pop = d$pop, aqua=d$dist_aqua), s = indexs))

stk.p <- inla.stack(tag = "pred", data = list(y = NA, numtrials = NA),
                    A = list(1, Ap), effects = list(data.frame(b0 = 1,
                    alt = dp$altitude, temp = dp$temp,
                    prec = dp$prec, hum = dp$hum, pop = dp$pop,
                    aqua=dp$dist_aqua), s = indexs))

stk.full <- inla.stack(stk.e, stk.p)

formula <- y ~ 0 + b0 + alt + temp + prec + hum + pop + aqua + f(s, model = spde)

res <- inla(formula, family = "binomial", Ntrials = numtrials,
            control.family = list(link = "logit"),
            data = inla.stack.data(stk.full),
            control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full)))


# Results

summary(res)

index <- inla.stack.index(stack = stk.full, tag = "pred")$data

prev_mean <- res$summary.fitted.values[index, "mean"]
prev_ll <- res$summary.fitted.values[index, "0.025quant"]
prev_ul <- res$summary.fitted.values[index, "0.975quant"]

pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent")

# Mean prevalence predicted values
leaflet() %>% addTiles() %>%
  addCircles(lng = coop[, 1], lat = coop[, 2], color = pal(prev_mean)) %>%
  addLegend("bottomright", pal = pal, values = prev_mean, title = "Mean Prev") %>%
  addScaleBar(position = c("bottomleft"))

# 0.025quant prevalence for predicted values
leaflet() %>% addTiles() %>%
  addCircles(lng = coop[, 1], lat = coop[, 2], color = pal(prev_ll)) %>%
  addLegend("bottomright", pal = pal, values = prev_ll, title = "0.025Quant Prev") %>%
  addScaleBar(position = c("bottomleft"))

# 0.975quant prevalence for predicted values
leaflet() %>% addTiles() %>%
  addCircles(lng = coop[, 1], lat = coop[, 2], color = pal(prev_ul)) %>%
  addLegend("bottomright", pal = pal, values = prev_ul, title = "0.975Quant Prev") %>%
  addScaleBar(position = c("bottomleft"))


# Exceedance probabilities

excprob <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1-inla.pmarginal(q = 0.2, marginal = marg)})

pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent")

leaflet() %>% addTiles() %>%
  addCircles(lng = coop[, 1], lat = coop[, 2], color = pal(excprob)) %>%
  addLegend("bottomright", pal = pal, values = excprob, title = "P(prev>0.2)" ) %>%
  addScaleBar(position = c("bottomleft"))

excprob <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1-inla.pmarginal(q = 0.4, marginal = marg)})

leaflet() %>% addTiles() %>%
  addCircles(lng = coop[, 1], lat = coop[, 2], color = pal(excprob)) %>%
  addLegend("bottomright", pal = pal, values = excprob, title = "P(prev>0.4)" ) %>%
  addScaleBar(position = c("bottomleft"))

excprob <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1-inla.pmarginal(q = 0.6, marginal = marg)})

leaflet() %>% addTiles() %>%
  addCircles(lng = coop[, 1], lat = coop[, 2], color = pal(excprob)) %>%
  addLegend("bottomright", pal = pal, values = excprob, title = "P(prev>0.6)" ) %>%
  addScaleBar(position = c("bottomleft"))

excprob <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1-inla.pmarginal(q = 0.8, marginal = marg)})

leaflet() %>% addTiles() %>%
  addCircles(lng = coop[, 1], lat = coop[, 2], color = pal(excprob)) %>%
  addLegend("bottomright", pal = pal, values = excprob, title = "P(prev>0.8)" ) %>%
  addScaleBar(position = c("bottomleft"))


# Plots for spatial field

rang <- apply(mesh$loc[, c(1, 2)], 2, range)
proj <- inla.mesh.projector(mesh, xlim = rang[, 1], ylim = rang[, 2], dims = c(300, 300))
mean_s <- inla.mesh.project(proj, res$summary.random$s$mean)
sd_s <- inla.mesh.project(proj, res$summary.random$s$sd)
df <- expand.grid(x = proj$x, y = proj$y)
df$mean_s <- as.vector(mean_s)
df$sd_s <- as.vector(sd_s)

gmean <- ggplot(df, aes(x = x, y = y, fill = mean_s)) + geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

gsd <- ggplot(df, aes(x = x, y = y, fill = sd_s)) + geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

plot_grid(gmean, gsd)




