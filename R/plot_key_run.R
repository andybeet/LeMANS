#' Plots output from the key_run script
#'
#' Plots catch, N, M1, M2, SSB, recruits output from the model.
#' This is a rudimentary tool used for troubleshooting the model.
#' The user is expected to develop their own set of plotting tools for publication purposes
#'
#'@param dataSet The data set to be plotted
#'@param ylabel The yaxis label for the figure
#'@param is.aggregated A logical value. Determines if aggregation is required over sizeClass
#'@param speciesNames A data frame in the format of \code{\link{rochet_GB_species}}
#'@param scales Text indicating is facet plots are displayed using same yaxis range ("fixed") or not ("free")
#'@param speciesSuitabilityPlot Either species number or species name. Only required if passing suitability or size preference data. Default = 1.
#'
#'
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 "aes" "geom_line" "facet_wrap" "ylab" "xlab"
#'
#'
#' @export

plot_key_run <- function(dataSet,ylabel="ylabel",is.aggregated=T,speciesNames,scales="fixed",speciesSuitability=1) {

  speciesSuitability <- error_check(speciesSuitability,speciesNames)

  dimOfData <- dim(dataSet)
  nSpecies <- dim(speciesNames)[1]
  df <- reshape::melt(dataSet)

  if ((length(dimOfData) > 2) & (dimOfData[2]==nSpecies)) { # output with size structure, by species over time (catch, N, M2)
    # transform data into a 2 array for use with ggplot
    names(df) <- c("sizeClass","speciesNumber","time","dataF")
    df <- dplyr::inner_join(df,speciesNames,by="speciesNumber")
    df <- dplyr::select(df,sizeClass,commonName,time,dataF)
    df$commonName <- factor(df$commonName, levels = unique(df$commonName))

    if (is.aggregated) {
      # sum over sizeclass to get species totals for each time step
      newdf <- dplyr::group_by(df,commonName,time) %>% dplyr::summarize(aggData = sum(dataF,na.rm=T))

      p <- ggplot2::ggplot(data=newdf) +
        geom_line(mapping =  aes(x = time, y = aggData)) + ylab(ylabel) +
        facet_wrap( ~ commonName, scales = scales)
      print(p)
    } else { # plot sizeclasses on same figure
      options(warn =-1)
      p <- ggplot2::ggplot(data=df) +
        geom_line(mapping =  aes(x = time, y = dataF, group = sizeClass)) + ylab(ylabel) +
        facet_wrap( ~ commonName, scales = scales)
      print(p)
      options(warn=0)

    }
  } else if ((length(dimOfData) > 2) & (dimOfData[1] != nSpecies)) { #  size preference and suitability
    # we plot each prey species and plot
    nSizeClasses <- dim(dataSet)[2]
    predSc <- data.frame(species = rep(c(1:nSpecies),each  = nSizeClasses),sizeClass=rep(1:nSizeClasses))
    addOnMAt <- apply(predSc,2,rep,nSpecies*nSizeClasses)
    # map first column (predator/size class to df)
    df <- cbind(predSc,df)
    names(df) <- c("speciesNumber","predSizeClass","index","preySizeClass","preySpecies","dataF")
    df <- dplyr::inner_join(df,speciesNames,"speciesNumber")
    df <- dplyr::select(df,commonName,predSizeClass,preySizeClass,preySpecies,dataF)
    df$commonName <- factor(df$commonName, levels = unique(df$commonName))

    # need to plot a facet plot for each species
    # see species suitability input
    dfTemp <- dplyr::filter(df,preySpecies==speciesSuitability)
       p <- ggplot2::ggplot(data=dfTemp) +
          geom_line(mapping =  aes(x = preySizeClass, y = dataF, group = predSizeClass)) +
          ylab(ylabel) + xlab(speciesNames$commonName[speciesSuitability]) +
          facet_wrap( ~ commonName, scales = scales)
        print(p)
       # options(warn=0)
       #

  } else if (length(dimOfData) == 2) { # output by species over time (SSB & recruits) or M1)
    if (dimOfData[1] == dim(speciesNames)[1]) { # nSpeces x  nYears
      if (is.aggregated) warning("Aggregregation has no efect for 2D data sets")
      names(df) <- c("speciesNumber","year","dataF")
      df <- dplyr::inner_join(df,speciesNames,"speciesNumber")
      df <- dplyr::select(df,commonName,year,dataF)
      df$commonName <- factor(df$commonName, levels = unique(df$commonName))


      p <- ggplot2::ggplot(data=df) +
        geom_line(mapping =  aes(x = year, y = dataF)) + ylab(ylabel) +
        facet_wrap( ~ commonName, scales = scales)
      print(p)
    } else { # M1
      names(df) <- c("sizeClass","speciesNumber","dataF")
      df <- dplyr::inner_join(df,speciesNames,"speciesNumber")
      df <- dplyr::select(df,commonName,sizeClass,dataF)
      df$commonName <- factor(df$commonName, levels = unique(df$commonName))
      df[df==0] <- NA
      options(warn=-1)
      p <- ggplot2::ggplot(data=df) +
        geom_line(mapping =  aes(x = sizeClass, y = dataF)) + ylab(ylabel) +
        facet_wrap( ~ commonName, scales = scales)
      print(p)
      options(warn=0)


    }



  } else {
    stop("Not coded for structures other than 2D and 3D arrays")
  }

}

# internal function - not exported

error_check <- function(speciesID,speciesData) {


  nSpecies <- dim(speciesData)[1]
  if (is.numeric(speciesID) ) {
    if (((speciesID %% 1) == 0)  &  ((speciesID > 0) & (speciesID<=nSpecies))) { # integer between 1 and nSpecies
      return(speciesID) # all is well
    } else {
      stop("speciesSuitability value not valid. Please input either a species name or a species number")
    }
  } else if (is.character(speciesID)) {
    if (any(tolower(speciesID)  %in%  tolower(speciesData$commonName))) {
      speciesID <- which(tolower(speciesID) == tolower(speciesData$commonName))
      return(speciesID)
    } else {
      stop("speciesSuitability value not valid. Please input either a species name or a species number")
    }

  } else {
    stop("speciesSuitability value not valid. Please input either a species name or a species number")
  }

}
