#' Takes a matrix of data for a species, checks if its numeric, then puts the table into a long-format dataframe
#'
#' @param x a matrix of data, generally species in the columns and sites in the row
#' @param row.metadata metadata for the sites; in long format, it will be stored in each row with with the site pertaining to the data
#' @param col.metadata metadata for the species; will be stored in every 'n'th row, where 'n' is the number of rows in the original table
#' @param total.metadata metadata for table; will include publishing information
#' @importFrom reshape2 melt
#' @return data set in long format, with all metadata included
.matrix.melt <- function(x, study.metadata=data.frame(units=NA, other=NA),
                         site.metadata=data.frame(id=NA,year=NA,name=NA,lat=NA,long=NA,address=NA,area=NA,other=NA),
                         species.metadata=data.frame(species=NA, taxonomy=NA, other=NA)){

    #######################
    # Argument handling ###
    #######################
    for(i in seq_along(names(study.metadata)))
        if(is.factor(study.metadata[,i]))
            study.metadata[,i] <- as.character(study.metadata[,i])
    for(i in seq_along(names(site.metadata)))
        if(is.factor(site.metadata[,i]))
            site.metadata[,i] <- as.character(site.metadata[,i])
    for(i in seq_along(names(species.metadata)))
        if(is.factor(species.metadata[,i]))
            species.metadata[,i] <- as.character(species.metadata[,i])
    if(!is.numeric(x))
        stop("'value' is not numeric")
    if(!is.matrix(x))
        stop("'x' is not a matrix")
    if(length(dim(x))!=2)
        stop("'x' is not a two-dimensional matrix")
    if(!identical(rownames(x), site.metadata$id))
        stop("Mismatch between site (names?) and site meta-data")
    if(!identical(colnames(x), species.metadata$species))
        stop("Mismatch between species (names?) and species meta-data")

    ######################
    # Dispatch    ########
    # to .df.melt ########
    # and return  ########
    ######################
    site.id <- rownames(x)[as.numeric(row(x))]
    species <- colnames(x)[as.numeric(col(x))]
    value <- as.numeric(x)
    return(.df.melt(species, site.id, value, study.metadata, site.metadata, species.metadata))
}

.df.melt <- function(species, site.id, value,
                     study.metadata=data.frame(units=NA, other=NA),
                     site.metadata=data.frame(id=NA,year=NA,name=NA,lat=NA,long=NA,address=NA,area=NA,other=NA),
                     species.metadata=data.frame(species=NA, taxonomy=NA, other=NA)){
    #######################
    # Argument handling ###
    #######################
    if(!is.numeric(value))
        stop("'value' is not numeric")
    if(any(is.na(value)))
        stop("No NAs in 'value'")
    if(any(is.na(species)))
        stop("No NAs in 'species'")
    if(any(is.na(site.id)))
        stop("No NAs in 'site.id'")
    species <- as.character(species)
    site.id <- as.character(site.id)

    ######################
    # Meta-data ##########
    ######################
    .create.other <- function(metadata, columns){
        if(!all(columns %in% names(metadata))){
            other <- metadata[,!names(metadata) %in% columns, drop=FALSE]
            metadata <- metadata[,names(metadata) %in% columns, drop=FALSE]        
            other <- sapply(seq_along(names(other)), function(y) paste(names(other)[y],other[,y],sep=":"))
            if(nrow(metadata) > 1)
                other <- paste(other, collapse=";") else other <- apply(other, 1, paste, collapse=";")
            metadata$other <- other
        } else {
            metadata$other <- NA
        }
        return(metadata[,c(columns,"other")])
    }
    # Study
    if(nrow(study.metadata) > 1)
        stop("Only one row of meta-data per study")
    if(!all("units" %in% names(study.metadata)))
        stop("Incorrectly formatted study meta-data")
    if(is.na(study.metadata$units))
        stop("Study must have units of measurement")
    study.metadata <- .create.other(study.metadata, "units")
    # Site
    if(!all(c("id","year","name","lat","long","address","area") %in% names(site.metadata)))
        stop("Incorrectly formatted site meta-data")
    if(length(intersect(unique(site.id), site.metadata$id)) != nrow(site.metadata))
        stop("Site meta-data must contain information about all sites")
    if(length(intersect(site.metadata$id,unique(site.id))) != nrow(site.metadata))
        stop("Site meta-data must only contain information about present sites")
    site.metadata <- .create.other(site.metadata, c("id","year","name","lat","long","address","area"))
    # Species
    if(!all(c("species","taxonomy") %in% names(species.metadata)))
        stop("Incorrectly formatted species meta-data")
    if(length(intersect(unique(species), species.metadata$species)) != nrow(species.metadata))
        stop("Species meta-data must contain information about all species")
    if(length(intersect(species.metadata$species,unique(species))) != nrow(species.metadata))
        stop("Species meta-data must only contain information about present species")
    species.metadata <- .create.other(species.metadata, c("species","taxonomy"))
    
    ######################
    # Format and return ##
    ######################
    # Reformat data
    output <- list(
        data=data.frame(site.id, species, value),
        spp.metadata=species.metadata,
        site.metadata=site.metadata,
        study.metadata=study.metadata
    )
    for(i in seq_along(output))
        for(j in seq_len(ncol(output[[i]])))
            if(is.factor(output[[i]][,j]))
                output[[i]][,j] <- as.character(output[[i]][,j])
    class(output) <- "nacdb"
    return(output)
}

# Takes a data already in long format that will be converted to a string of metadata. Each row will be a single string, and the
# function will return the list of these strings
#
# @param data a dataframe exclusively containing the columns of metadata
# @return a list of metadata strings
.make.metadata <- function(data){
  sapply(1:nrow(data), function(y) {
    char.list <- c(rbind(colnames(data), "=", as.character(data[y,]), ", "))
    char.list <- head(char.list, -1)
    metadata <- paste(char.list, collapse="")
    return(metadata)
  })
}

# Unzips a file from a downloaded zip file
# param file name of file to be extracted from zip
# param zip location and name of zip file (e.g.,
#     ~/Downlaods/a_file.zip)
# param to.save.dir directory to save resulting file (DEFAULT: a new
#     temporary directory will be used)
# param to.save.name name to save the file as (DEFAULT: it will be
#     named paste(zip,file, sep='_'))
# return Complete path to unzipped file
#' @importFrom utils unzip download.file
#' @importFrom reshape2 melt
#' @importFrom httr GET
#' @importFrom stats setNames
.unzip <- function(file, zip, to.save.dir, to.save.name){
    if(missing(to.save.dir))
        to.save.dir <- tempdir()
    if(missing(to.save.name))
        to.save.name <- file
    
    files <- unzip(zip, list=TRUE)
    if(!file %in% files$Name)
        stop("Required file not in zipfile ", zip)

    file <- unzip(zip, file)
    file.rename(file, file.path(to.save.dir, to.save.name))
    return(file.path(to.save.dir, to.save.name))
}

.fac.sim <- function(x){
    x <- Filter(Negate(is.na), x)
    x <- x[x != "" & x != " "]
    x <- unique(x)
    return(paste(x,collapse="_"))
}

#' @importFrom stats model.matrix
.expand.factor <- function(factor_to_expand, name){
    names <- rep(name, length(unique(factor_to_expand)))
    output <- model.matrix(~factor_to_expand-1)
    colnames(output) <- paste(names, gsub("factor_to_expand", "", colnames(output)), sep="_")
    return(as.data.frame(output))
}

.download <- function(url, dir, save.name, cache=TRUE){
    destination <- file.path(dir, save.name)
    suffix <- .file.suffix(url, 4)

    if(cache==TRUE & file.exists(destination)){
        if(!is.na(suffix))
            attr(destination, "suffix") <- suffix
        return(destination)
    }

    result <- download.file(url, destination, quiet=TRUE)
    if(result != 0)
        stop("Error code", result, " downloading file; file may not exist")

    if(!is.na(suffix))
        attr(destination, "suffix") <- suffix
    return(destination)
}

.save.name <- function(doi, save.name, file){
    if(is.na(save.name)){
        save.name <- paste(doi,file, sep="_")
        save.name <- gsub(.Platform$file.sep, "_", save.name, fixed=TRUE)
    }
    return(save.name)
}

.grep.url <- function(url, regexp, which=1){
    html <- as.character(GET(url))
    return(.grep.text(html, regexp, which))
}

.grep.text <- function(text, regexp, which=1){
    links <- gregexpr(regexp, text)
    if(which > length(links[[1]]))
        stop("SI number '", which, "' greater than number of detected SIs (", length(links[[1]]), ")")
    pos <- as.numeric(links[[1]][which])
    return(substr(text, pos, pos+attr(links[[1]], "match.length")[which]-1))
}

.file.suffix <- function(text, max.length=4){
    suffix <- .grep.text(text, "[a-zA-Z]+$")
    if(nchar(suffix) <= max.length & nchar(suffix) > 0)
        return(suffix)
    return(NA)
}

prog.bar <- function(x, y){
    if(y < 100){
        cat(".")} else {
            z <- Filter(function(z) z>=0, seq(1,y,length.out=100)-x)
            if(length(z) > 0)
                tryCatch(if(z[1] < 1) if((length(z) %% 10)==0) cat("|") else cat("."), error=function(z) cat("."))
        }
}    
