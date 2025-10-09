par.start_weibweib <- list(mu.0 = 0.3, sigma.0 = 0.5,
    mu.1 = 0, sigma.1 = 1, mu.2 = 0, sigma.2 = 1,
    mu.3 = 0, sigma.3 = 1, alpha.d = 1.2, beta.d = 0,
    mu.d1 = 0, sigma.d1 = 1, mu.d2 = 0, sigma.d2 = 1,
    mu.d3 = 0, sigma.d3 = 1, shape.y1nd = 1.5, scale.y1nd = 1,
    mu.y1nd = -4, sigma.y1nd = 1, shape.y1d = 1.5,
    scale.y1d = 1, mu.y1d = -4, sigma.y1d = 1, shape.y0nd = 1.5,
    scale.y0nd = 1, mu.y0nd = -4, sigma.y0nd = 1,
    shape.y0d = 1.5, scale.y0d = 1, mu.y0d = -4,
    sigma.y0d = 1, mu.y1 = 0, sigma.y1 = 1, mu.y2 = 0,
    sigma.y2 = 1, mu.y3 = 0, sigma.y3 = 1, mu.delta = 0,
    sigma.delta = 1)

# set the parameters priors' values
theta.prior_weibweib <- list(mu.0 = 0, sigma2.0 = 1,
    mu.1 = 0, sigma2.1 = 5^2, mu.2 = 0, sigma2.2 = 5^2,
    mu.3 = 0, sigma2.3 = 5^2, a.d = 1, b.d = 1,
    mu.d = 0, sigma2.d = 10^2, mu.d1 = 0, sigma2.d1 = 5^2,
    mu.d2 = 0, sigma2.d2 = 5^2, mu.d3 = 0, sigma2.d3 = 5^2,
    a.y1nd = 1, b.y1nd = 1, mu.y1nd = 0, sigma2.y1nd = 10^2,
    a.y1d = 1, b.y1d = 1, mu.y1d = 0, sigma2.y1d = 10^2,
    a.y0nd = 1, b.y0nd = 1, mu.y0nd = 0, sigma2.y0nd = 10^2,
    a.y0d = 1, b.y0d = 1, mu.y0d = 0, sigma2.y0d = 10^2,
    mu.y1 = 0, sigma2.y1 = 5^2, mu.y2 = 0, sigma2.y2 = 5^2,
    mu.y3 = 0, sigma2.y3 = 5^2, mu.delta = 0, sigma2.delta = 10^2)

# set the proposals' values
proposal_weibweib <- list(sd.0 = 0.1, sd.1 = 0.1,
    sd.2 = 0.1, sd.3 = 0.1, scale.d = 0.004, sd.d = 0.1,
    sd.d1 = 0.1, sd.d2 = 0.1, sd.d3 = 0.1, scale.y1nd = 0.01,
    sd.y1nd = 0.1, scale.y1d = 0.01, sd.y1d = 0.1,
    scale.y0nd = 0.01, sd.y0nd = 0.1, scale.y0d = 0.01,
    sd.y0d = 0.1, sd.y1 = 0.1, sd.y2 = 0.1, sd.y3 = 0.1,
    sd.delta = 0.1)

