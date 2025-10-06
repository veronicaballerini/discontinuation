dweib <- function(x, a, b) {

    {
        a * x^{
            a - 1
        }
    } * exp(b - exp(b) * x^{
        a
    }) * (1 - (x < 0))

}

hweib <- function(x, a, b) {

    {
        a * x^{
            a - 1
        }
    } * exp(b)

}


Sweib <- function(x, a, b) {

    exp(-exp(b) * x^{
        a
    })

}

Eweib <- function(alpha, beta) {
    ev <- exp(-beta/alpha) * gamma(1 + 1/alpha)
    return(ev)
}


dtweib <- function(x, sh, sc, a, b = Inf) {

    tt <- rep(0, length(x))

    Ga <- pweibull(a, shape = sh, scale = sc)
    Gb <- pweibull(b, shape = sh, scale = sc)

    tt[x >= a & x <= b] <- dweibull(x[x >= a & x <= b], shape = sh, scale = sc)/(Gb - Ga)

    return(tt)
}

qtweibull <- function(p, shape, scale, a = 0, b = Inf) {
    stopifnot(all(p >= 0 & p <= 1) & all(scale > 0) & all(shape > 0))
    Fa <- pweibull(a, shape, scale)
    Fb <- pweibull(b, shape, scale)
    pNew <- p * (Fb - Fa) + Fa
    x <- qweibull(pNew, shape, scale)
    return(x)
}

# qtweib <- function(p, sh, sc, a, b){ G <- pweibull(a, shape=sh, scale=sc) + p*{pweibull(b, shape=sh, scale=sc)-pweibull(a, shape=sh, scale=sc)} Gin <-
# qweibull(G, shape=sh, scale=sc) result <- pmin(pmax(a, Gin), b) return(result) } Stweib <- function(x, sh, sc, a){ 1-ptweibull(x,sh,sc,a) }

Stweib <- function(x, a, b, l) {

    exp(exp(b) * (l^{
        a
    } - x^{
        a
    }))

}

ptweibull <- function(q, shape, scale, a = 0, b = Inf) {
    stopifnot(all(scale > 0) & all(shape > 0))
    Fa <- pweibull(a, shape, scale)
    Fb <- pweibull(b, shape, scale)
    p <- (pweibull(q, shape, scale) - Fa)/(Fb - Fa)
    inda <- which(q < a)
    indb <- which(q > b)
    if (length(inda) > 0)
        p[inda] <- 0
    if (length(indb) > 0)
        p[indb] <- 1
    return(p)
}

rtweibull <- function(n, shape, scale, a = 0, b = Inf) {
    stopifnot(n > 0 & all(scale > 0) & all(shape > 0))
    x <- runif(n)
    Fa <- pweibull(a, shape, scale)
    Fb <- pweibull(b, shape, scale)
    y <- (1 - x) * Fa + x * Fb
    return(qweibull(y, shape, scale))
}

# rtweib <- function(n, sh, sc, a, b=Inf){ u <- runif(n, min = 0, max = 1) x <- qtweib(u, sh, sc, a, b) return(x) }

Igamma <- function(z, u) pgamma(u, z) * gamma(z)

Weib.Trunc.Moments <- function(order = 1, shape, scale = 1, location = 0, left = 0, right = Inf, central = FALSE) {
    ## Parameters order: moment order shape: shape paameter scale: scale parameter location: location parameter, if location=0 two-parameter Weibull moment
    ## left: left bound right: right bound central: if TRUE compute the central moment, if FALSE compute the noncentral moment
    if (left < location) {
        left <- location
    }
    if (left >= right)
        stop("'left' must be strictly greater than 'right'")
    m <- NULL
    n = order
    ## noncentral moment
    if (!central) {
        vect1 <- choose(n, 0:n) * (location^(0:n)) * scale^(n:0) * (Igamma((n:0)/shape + 1, ((right - location)/scale)^shape) - Igamma((n:0)/shape + 1, ((left -
            location)/scale)^shape))
        const1 <- exp(-((left - location)/scale)^shape) - exp(-((right - location)/scale)^shape)
        m = sum(vect1)/const1
    }
    ## central moment
    if (central) {
        vect1 <- choose(n, 0:n) * ((Igamma(1/shape + 1, ((right - location)/scale)^shape) - Igamma(1/shape + 1, ((left - location)/scale)^shape))/(exp(-((right -
            location)/scale)^shape) - exp(-((left - location)/scale)^shape)))^(0:n)
        vect2 <- Igamma((n:0)/shape + 1, ((right - location)/scale)^shape) - Igamma((n:0)/shape + 1, ((left - location)/scale)^shape)
        const1 <- scale^n/(exp(-((left - location)/scale)^shape) - exp(-((right - location)/scale)^shape))
        m = sum(vect1 * vect2) * const1
    }
    return(m)
}

myload <- function(x) {
    x <- load(x)
    get(x)
}

create_covariate_grid_general <- function(continuous_vars = list(x1 = c(-2.9, 2.9, 0.1)), binary_count = 2) {
    # Generate sequences for each continuous variable
    cont_seq_list <- lapply(continuous_vars, function(r) seq(from = r[1], to = r[2], by = r[3]))

    # Name the continuous variables if not already named
    if (is.null(names(cont_seq_list))) {
        names(cont_seq_list) <- paste0("x", seq_along(cont_seq_list))
    }

    # Create grid of continuous covariates
    cont_grid <- expand.grid(cont_seq_list)

    # Create all combinations of binary covariates (0/1)
    if (binary_count > 0) {
        binary_combinations <- expand.grid(rep(list(0:1), binary_count))
        names(binary_combinations) <- paste0("x", (ncol(cont_grid) + 1):(ncol(cont_grid) + binary_count))

        # Repeat continuous grid for each binary combination and bind
        final_grid <- do.call(rbind, lapply(1:nrow(binary_combinations), function(i) {
            cbind(cont_grid, binary_combinations[rep(i, nrow(cont_grid)), ])
        }))
    } else {
        final_grid <- cont_grid
    }

    rownames(final_grid) <- NULL
    return(final_grid)
}

hdi99 <- function(obj) {
    hdi(obj, 0.99)
}
hdi95 <- function(obj) {
    hdi(obj, 0.95)
}
