clean.natdb <- function(x){
    # Argument handling
    if(!inherits(x, "nasdb"))
        stop("'", deparse(substitute(x)), "' must be of type 'nasdb'")

    # Clean up any obvious weirdnesses with the site names

    # Clean up any obvious weirdnesses with the species names

    # Clean up any obvious weirdnesses with the abundances/counts/pres-abs
    return(x)
}

# This does more thorough cleaning of species names. Will likely
# require the addition of some sort of cache, as there will be
# *thousands* of species names that need adding in here
clean.nasdb.names <- function(x, thresh, ...){
    # Argument handling
    if(!inherits(x, "natdb"))
        stop("'", deparse(substitute(x)), "' must be of type 'natdb'")

    # This code doesn't work on a nasdb object, probably, but the general structure will
    spp <- unique(c(unique(x$numeric$species), unique(x$categorical$species)))
    dwn.spp <- gnr_resolve(spp)
    dwn.spp <- dwn.spp[!duplicated(dwn.spp$user_supplied_name),]
    dwn.spp$matched_name <- tolower(sapply(strsplit(dwn.spp$matched_name, " "), function(x) paste(x[1:2],collapse="_")))
    
    if(!missing(thresh))
        dwn.spp <- dwn.spp[dwn.spp$score >= thresh,]
    lookup <- with(dwn.spp, setNames(matched_name, user_supplied_name))
    
    x$numeric$species <- lookup[x$numeric$species]
    x$categorical$species <- lookup[x$categorical$species]
    return(x)
}
