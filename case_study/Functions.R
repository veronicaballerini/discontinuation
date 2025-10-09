Sexp <- function(x, rate) {

    exp(-rate * x)

}

Stexp <- function(x, rate, a = 0, b = Inf) {

    (exp(-rate * x) - exp(-rate * b))/(exp(-rate *
        a) - exp(-rate * b))

}

Etexp <- function(rate, a = 0, b = Inf) {

    (a * exp(-rate * a) - b * exp(-rate * b))/(exp(-rate *
        a) - exp(-rate * b)) + 1/rate

}

rtexp <- function(n, rate, a = 0, b = Inf) {
    stopifnot(n > 0 & all(rate > 0))
    x <- runif(n)
    Fa <- pexp(a, rate)
    Fb <- pexp(b, rate)
    y <- (1 - x) * Fa + x * Fb
    return(qexp(y, rate))
}

dtexp <- function(x, rate, a = 0, b = Inf) {
    stopifnot(all(rate > 0))
    Fa <- pexp(a, rate)
    Fb <- pexp(b, rate)
    y <- dexp(x, rate)
    inda <- which(x < a)
    indb <- which(x > b)
    if (length(inda) > 0)
        y[inda] <- 0
    if (length(indb) > 0)
        y[indb] <- 0
    return(y/(Fb - Fa))
}

dweib <- function(x, a, b) {

    (a * x^(a - 1)) * exp(b - exp(b) * x^a) * (1 -
        (x < 0))

}

# hweib <- function(x, a, b){
# (a*x^(a-1))*exp(b) }


Sweib <- function(x, a, b) {

    exp(-exp(b) * x^a)

}

# Eweib <- function(alpha, beta){ ev <-
# exp(-beta/alpha)*gamma(1+1/alpha) return(ev)
# }


dtweib <- function(x, sh, sc, a, b = Inf) {

    tt <- rep(0, length(x))

    Ga <- pweibull(a, shape = sh, scale = sc)
    Gb <- pweibull(b, shape = sh, scale = sc)

    tt[x >= a & x <= b] <- dweibull(x[x >= a & x <=
        b], shape = sh, scale = sc)/(Gb - Ga)

    return(tt)
}

Stweib <- function(x, a, b, l) {

    exp(exp(b) * (l^a - x^a))

}

rtgamma <- function(n, shape, scale, a = 0, b = Inf) {
    # CDF values at lower and upper bounds
    plower <- pgamma(a, shape = shape, scale = scale)
    pupper <- pgamma(b, shape = shape, scale = scale)

    # Draw uniform samples between CDF(lower)
    # and CDF(upper)
    u <- runif(n, min = plower, max = pupper)

    # Inverse CDF (quantile function) to get
    # the truncated samples
    qgamma(u, shape = shape, scale = scale)
}

dtgamma <- function(x, shape, scale, a = 0, b = Inf,
    log = FALSE) {
    # Density at x
    dens <- dgamma(x, shape = shape, scale = scale)
    # Normalization constant
    norm_const <- pgamma(b, shape = shape, scale = scale) -
        pgamma(a, shape = shape, scale = scale)
    # Truncated density (set to 0 outside
    # [lower, upper])
    dens[x < a | x > b] <- 0
    # Return (log-)density
    if (log) {
        out <- rep(-Inf, length(x))
        inside <- (x >= a & x <= b)
        out[inside] <- log(dgamma(x[inside], shape = shape,
            scale = scale)) - log(norm_const)
        return(out)
    } else {
        return(dens/norm_const)
    }
}

myload <- function(x) {
    x <- load(x)
    get(x)
}

hdi99 <- function(obj) {
    hdi(obj, 0.99)
}
hdi95 <- function(obj) {
    hdi(obj, 0.95)
}

create_covariate_grid_general <- function(continuous_vars = list(x1 = c(-2.9,
    2.9, 0.1)), binary_count = 2) {
    # Generate sequences for each continuous
    # variable
    cont_seq_list <- lapply(continuous_vars, function(r) seq(from = r[1],
        to = r[2], by = r[3]))

    # Name the continuous variables if not
    # already named
    if (is.null(names(cont_seq_list))) {
        names(cont_seq_list) <- paste0("x", seq_along(cont_seq_list))
    }

    # Create grid of continuous covariates
    cont_grid <- expand.grid(cont_seq_list)

    # Create all combinations of binary
    # covariates (0/1)
    if (binary_count > 0) {
        binary_combinations <- expand.grid(rep(list(0:1),
            binary_count))
        names(binary_combinations) <- paste0("x",
            (ncol(cont_grid) + 1):(ncol(cont_grid) +
                binary_count))

        # Repeat continuous grid for each
        # binary combination and bind
        final_grid <- do.call(rbind, lapply(1:nrow(binary_combinations),
            function(i) {
                cbind(cont_grid, binary_combinations[rep(i,
                  nrow(cont_grid)), ])
            }))
    } else {
        final_grid <- cont_grid
    }

    rownames(final_grid) <- NULL
    return(final_grid)
}

Igamma <- function(z, u) pgamma(u, z) * gamma(z)

Weib.Trunc.Moments <- function(order = 1, shape,
    scale = 1, location = 0, left = 0, right = Inf,
    central = FALSE) {
    if (left < location) {
        left <- location
    }
    if (left >= right) {
        stop("'left' must be strictly less than 'right'")
    }

    m <- NULL
    n <- order
    k_range <- 0:n

    # Non-central moment
    if (!central) {
        vect1 <- choose(n, k_range) * location^k_range *
            rep(scale, length(k_range))^(n - k_range) *
            (Igamma((n - k_range)/shape + 1, ((right -
                location)/scale)^shape) - Igamma((n -
                k_range)/shape + 1, ((left - location)/scale)^shape))

        const1 <- exp(-((left - location)/scale)^shape) -
            exp(-((right - location)/scale)^shape)

        m <- sum(vect1)/const1
    }

    # Central moment
    if (central) {
        mean_weibull <- (Igamma(1/shape + 1, ((right -
            location)/scale)^shape) - Igamma(1/shape +
            1, ((left - location)/scale)^shape))/(exp(-((right -
            location)/scale)^shape) - exp(-((left -
            location)/scale)^shape))

        vect1 <- choose(n, k_range) * mean_weibull^k_range

        vect2 <- Igamma((n - k_range)/shape + 1,
            ((right - location)/scale)^shape) -
            Igamma((n - k_range)/shape + 1, ((left -
                location)/scale)^shape)

        const1 <- rep(scale, length(k_range))^n/(exp(-((left -
            location)/scale)^shape) - exp(-((right -
            location)/scale)^shape))

        m <- sum(vect1 * vect2) * const1
    }

    return(m)
}

myload <- function(x) {
    x <- load(x)
    get(x)
}

hdi99 <- function(obj) {
    hdi(obj, 0.99)
}
hdi95 <- function(obj) {
    hdi(obj, 0.95)
}

create_covariate_grid_general <- function(continuous_vars = list(x1 = c(-2.9,
    2.9, 0.1)), binary_count = 2) {
    # Generate sequences for each continuous
    # variable
    cont_seq_list <- lapply(continuous_vars, function(r) seq(from = r[1],
        to = r[2], by = r[3]))

    # Name the continuous variables if not
    # already named
    if (is.null(names(cont_seq_list))) {
        names(cont_seq_list) <- paste0("x", seq_along(cont_seq_list))
    }

    # Create grid of continuous covariates
    cont_grid <- expand.grid(cont_seq_list)

    # Create all combinations of binary
    # covariates (0/1)
    if (binary_count > 0) {
        binary_combinations <- expand.grid(rep(list(0:1),
            binary_count))
        names(binary_combinations) <- paste0("x",
            (ncol(cont_grid) + 1):(ncol(cont_grid) +
                binary_count))

        # Repeat continuous grid for each
        # binary combination and bind
        final_grid <- do.call(rbind, lapply(1:nrow(binary_combinations),
            function(i) {
                cbind(cont_grid, binary_combinations[rep(i,
                  nrow(cont_grid)), ])
            }))
    } else {
        final_grid <- cont_grid
    }

    rownames(final_grid) <- NULL
    return(final_grid)
}
