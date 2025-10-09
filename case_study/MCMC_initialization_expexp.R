# set starting values
par.start_expexp <- list(mu.0 = 0.4, sigma.0 = 1,
    mu.1 = 0, sigma.1 = 1, mu.2 = 0, sigma.2 = 1,
    mu.3 = 0, sigma.3 = 1, beta.d = 0, mu.d1 = 0,
    sigma.d1 = 1, mu.d2 = 0, sigma.d2 = 1, mu.d3 = 0,
    sigma.d3 = 1, mu.y1nd = 0, sigma.y1nd = 1, mu.y1d = 0,
    sigma.y1d = 1, mu.y0nd = 0, sigma.y0nd = 1,
    mu.y0d = 0, sigma.y0d = 1, mu.y1 = 0, sigma.y1 = 1,
    mu.y2 = 0, sigma.y2 = 1, mu.y3 = 0, sigma.y3 = 1,
    mu.delta = 0, sigma.delta = 1)

# set the parameters priors' values
theta.prior_expexp <- list(mu.0 = 0, sigma2.0 = 1,
    mu.1 = 0, sigma2.1 = 5^2, mu.2 = 0, sigma2.2 = 5^2,
    mu.3 = 0, sigma2.3 = 5^2, mu.d = 0, sigma2.d = 10^2,
    mu.d1 = 0, sigma2.d1 = 5^2, mu.d2 = 0, sigma2.d2 = 5^2,
    mu.d3 = 0, sigma2.d3 = 5^2, mu.y1nd = 0, sigma2.y1nd = 10^2,
    mu.y1d = 0, sigma2.y1d = 10^2, mu.y0nd = 0,
    sigma2.y0nd = 10^2, mu.y0d = 0, sigma2.y0d = 10^2,
    mu.y1 = 0, sigma2.y1 = 5^2, mu.y2 = 0, sigma2.y2 = 5^2,
    mu.y3 = 0, sigma2.y3 = 5^2, mu.delta = 0, sigma2.delta = 10^2)

# set the proposals' values
proposal_expexp <- list(sd.0 = 0.2, sd.1 = 0.15,
    sd.2 = 0.15, sd.3 = 0.15, sd.d = 0.2, sd.d1 = 0.15,
    sd.d2 = 0.15, sd.d3 = 0.15, sd.y1nd = 0.2, sd.y1d = 0.2,
    sd.y0nd = 0.2, sd.y0d = 0.2, sd.y1 = 0.1, sd.y2 = 0.1,
    sd.y3 = 0.1, sd.delta = 0.1)
