#' Builds a soil database
#'
#' The key function of the nasdb package. When run with defaults, it
#' will download and build a database of species' traits from all the
#' manuscript sources in the package. This totals XXX
#' manuscripts/databases, XXX species, and XXX traits. Please note
#' that all parameters are interactive; thus specifying \code{species}
#' and \code{traits} constraints will constraint according to both,
#' for example. Please also note that specifying any kind of
#' constraints makes use of the package's built-in cache of what
#' species and traits information are available in each database;
#' making use of this on the GitHub (developer) build of this package
#' is not advisable, and (further) it is impossible for us to verify
#' whether the datasets NATDB searches have been updated since the
#' package was last built.
#' 
#' @param datasets Character vector of datasets to be searched for
#'     trait data. If not specified (the default) all trait datasets
#'     will be downloaded and returned.
#' @param species Character vector of species to be searched for trait
#'     data. If not specified (the default) data for all species will
#'     be downloaded and returned.
#' @param traits Character vector of traits to be searched for
#'     data. If not specified (the default) data for all traits will
#'     be downloaded and returned.
#' @return nasdb.data object. XXX
#' @author Will Pearse; Clint; Raimi; Ethan; etc.
#' #@examples
#' # Limit the scope of these as they have to work online on servers!...
#' #@seealso 
#' @export

nasdb <- function(cache, datasets, delay=5){
    #Check datasets
    if(missing(datasets)){
        datasets <- Filter(Negate(is.function), ls(pattern="^\\.[a-z]*\\.[0-9]+", name="package:nasdb", all.names=TRUE))
    } else {
        datasets <- paste0(".", tolower(datasets))
        datasets <- gsub("..", ".", datasets, fixed=TRUE)
    }
    if(!all(datasets %in% datasets)){
        missing <- setdiff(datasets, ls.funs())
        stop("Error: ", paste(missing, collapse=", "), "not in nasdb")
    }
    
    #Do data loads
    output <- vector("list", length(datasets))
    for(i in seq_along(datasets)){
        prog.bar(i, length(datasets))
        if(!missing(cache)){
            if(!file.exists(cache))
                stop("Cache directory does not exist")
            path <- file.path(cache,paste0(datasets[i], ".RDS"))
        } else path <- NA
        if(!is.na(path) && file.exists(path)){
            output[[i]] <- readRDS(path)
            next()
        }
        
        output[[i]] <- eval(as.name(datasets[i]))()
        
        output[[i]]$data$study <- datasets[i]
        output[[i]]$spp.metadata$study <- datasets[i]
        output[[i]]$site.metadata$study <- datasets[i]
        output[[i]]$study.metadata$study <- datasets[i]
        output[[i]]$data$site.id <- paste0(output[[i]]$data$site.id,datasets[i])
        output[[i]]$site.metadata$id <- paste0(output[[i]]$site.metadata$id,datasets[i])
        
        if(!is.na(path))
            saveRDS(output[[i]], path)
        Sys.sleep(delay)
    }
    
    # Merge data and return
    output <- list(
        data=do.call(rbind, lapply(output, function(x) x$data)),
        spp.metadata=do.call(rbind, lapply(output, function(x) x$spp.metadata)),
        site.metadata=do.call(rbind, lapply(output, function(x) x$site.metadata)),
        study.metadata=do.call(rbind, lapply(output, function(x) x$study.metadata))
    )
    class(output) <- "nasdb"
    return(output)
}

print.nasdb <- function(x, ...){
    # Argument handling
    if(!inherits(x, "nasdb"))
        stop("'", deparse(substitute(x)), "' must be of type 'nasdb'")
    
    # Create a simple summary matrix of species and sites in x
    n.species <- length(unique(species(x)))
    n.sites <- length(unique(sites(x)))
    n.total <- nrow(x$data)
    
    # Print it to screen
    cat("\nA Community DataBase containing:\nSpecies  : ", n.species, "\nSites    : ", n.sites, "\nTotal    : ", n.total,"\n")
    invisible(setNames(c(n.species,n.sites), c("n.species","n.sites")))
}

summary.nasdb <- function(x, ...){
    print.nasdb(x, ...)
}

"[.nasdb" <- function(x, sites, spp){
    # Argument handling
    if(!inherits(x, "nasdb"))
        stop("'", deparse(substitute(x)), "' must be of type 'nasdb'")

    # Setup null output in case of no match
    null <- list(
        data=data.frame(species=NA,site.id=NA,value=NA),
        study.metadata=data.frame(units=NA,other=NA),
        site.metadata=data.frame(id=NA,year=NA,name=NA,lat=NA,long=NA,address=NA,other=NA),
        spp.metadata=data.frame(species=NA, taxonomy=NA, other=NA)
    )
    class(null) <- "nasdb"

    # Site subsetting
    if(!missing(sites)){
        if(any(x$site.metadata$id %in% sites)){
            x$data <- x$data[x$data$site.id %in% sites,]
            x$spp.metadata <- x$spp.metadata[x$spp.metadata %in% x$data$species,]
            x$site.metadata <- x$site.metadata[x$site.metadata$id %in% sites,]
            x$study.metadata <- x$study.metadata[x$study.metadata %in% x$data$study,]
        } else {
            return(null)
        }
    }
    
    # Species subsetting
    if(!missing(spp)){
        if(any(x$spp.metadata$species %in% spp)){
            x$data <- x$data[x$data$species %in% spp,]
            x$spp.metadata <- x$spp.metadata[x$spp.metadata %in% spp,]
            x$site.metadata <- x$site.metadata[x$site.metadata %in% x$data$site,]
            x$study.metadata <- x$study.metadata[x$study.metadata %in% x$data$study,]
        } else {
            return(null)
        }
    }

    # Return (already checked for null case)
    return(output)
}

species <- function(x, ...){
    if(!inherits(x, "nasdb"))
        stop("'", deparse(substitute(x)), "' must be of type 'nasdb'")
    return(unique(x$spp.metadata$species))
    # Return a vector of the sites in nasdb (?)
}

sites <- function(x, ...){
    if(!inherits(x, "nasdb"))
        stop("'", deparse(substitute(x)), "' must be of type 'nasdb'")
    return(unique(x$site.metadata$id))
}

citations <- function(x){
    if(!inherits(x, "nasdb"))
        stop("'", deparse(substitute(x)), "' must be of type 'nasdb'")
    
    data(nasdb_citations)
    datasets <- Filter(Negate(is.function), ls(pattern="^\\.[a-z]*\\.[0-9]+[a-d]?", name="package:nasdb", all.names=TRUE))
    nasdb.citations$Name <- with(nasdb.citations, paste0(".", tolower(Author), ".", Year))

    return(as.character(nasdb.citations$BibTeX.citation[match(datasets, nasdb.citations$Name)]))
}
