"geoXY" <-
function(latitude, longitude,
         lat0 = min(latitude, na.rm=TRUE),
         lon0 = min(longitude, na.rm=TRUE),
         unit = 1.) {
    lat0 <- rep(lat0, length(latitude))
    lon0 <- rep(lon0, length(longitude))
    yDist <- geoDist(lat0, lon0, latitude, lon0)
    xDist <- geoDist(lat0, lon0, lat0, longitude)
    cbind(X= xDist * sign(longitude - lon0)/unit,
              Y = yDist * sign(latitude - lat0)/unit)
}

# ========

"geoDist" <-
function (lat1, lon1, lat2, lon2, NAOK = TRUE, DUP = TRUE) 
# Note de PL: Éliminer l'argument DUP, "For back-compatibility, accepted but ignored"
{
    n <- unique(c(length(lat1), length(lon1), length(lat2), length(lon2)))
    nok <- n
    if (length(n) > 1) 
        stop("Need all arguments of the same length:  got ", 
            paste(n, collapse = ", "))
    if (n < 1) 
        return(numeric())
    nas <- is.na(lat1) | is.na(lat2) | is.na(lon1) | is.na(lon2)
    if (NAOK) {
        if (any(nas)) {
            ok <- !nas
            lat1 <- lat1[ok]
            lon1 <- lon1[ok]
            lat2 <- lat2[ok]
            lon2 <- lon2[ok]
            nok <- sum(ok)
        }
    }
    else if (any(nas)) 
        stop("NA values found but not allowed")
    res <- .Fortran("GEODISTV", as.double(lat1), as.double(lon1), 
        as.double(lat2), as.double(lon2), dist = double(nok), 
        as.integer(nok), PACKAGE = "SoDA")$dist
    if (NAOK && any(nas)) {
        value <- rep(NA, n)
        value[ok] <- res
        value
    }
    else res
}

# ========