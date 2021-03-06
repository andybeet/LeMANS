#' Plots output from the key_run script
#'
#' Plots catch, N, M1, M2, SSB, recruits, suitability, growth efficiency, ration, maturity, and growth output from the model run.
#' See \code{\link{key_run}} for details regarding obtaining output. This is a rudimentary tool used for troubleshooting the model.
#' The user is expected to develop their own set of plotting tools for publication purposes
#'
#'@param dataSet The data set to be plotted. See \code{\link{key_run}}
#'@param ylabel The yaxis label for the figure. Match this with the data.
#'@param is.aggregated A logical value. Determines if aggregation is required over sizeClass. Only valid for N, M2 and catch.
#'@param speciesNames A data frame in the format of \code{\link{rochet_GB_species}}
#'@param scales Text indicating whether to display yaxis range fixed for all species ("fixed") or not ("free").
#'@param speciesSuitability Either species number or species name. Only required if passing suitability or size preference data. Default = 1.
#'@param predOrPrey character string indicating whether to plot suitabiliy for species as the prey or predator. Values "prey" or "pred". Default = "Prey"
#'
#'@section Notes:
#'
#'Not all arguments are required for plotting.
#'\code{speciesSuitability} and \code{predOrPrey} are only required for plotting suitabilities (and size preferences).
#'If you are plotting any other output data then you can ignore these arguments.
#'In some cases \code{is.aggregated} is irrelevant. A warning will let you know that this is so.
#'
#'The suitability plots (similar for size preference plots) are interpreted in the following way:
#'If \code{predOrPrey} = "prey" and \code{speciesSuitability} = "Atlantic cod" then each facet represents a predator of cod.
#'Each line drawn represents a size class of the predator and displays its suitability for each size class (x axis) of Cod (prey)
#'
#'The units of the output variable are explained in \code{\link{key_run}}.
#'
#'
#''@examples
#'\dontrun{
#'# runs the model with bundled data from Rochet et al (2011).
#' output <- key_run(Ffull=.4,nYrs=50,rochet_GB_modelSetup,rochet_GB_parameterValues,rochet_GB_initialValues,rochet_GB_foodweb,rochet_GB_species)
#'
#'# plots suitability of cod as a prey species to all other species in the model.
#'plot_key_run(output$suitability,ylabel="Suitability",speciesNames=rochet_GB_species,scales="free",speciesSuitability = "atLANTIC coD",predOrPrey="prey")
#'
#'# plots catch (in millions of individuals) aggregated over size class for each species in the model.
#'plot_key_run(output$catch/1E6,ylabel="Catch (millions of individuals)",is.aggregated = T,speciesNames=rochet_GB_species,scales="free")
#'
#'}
#'
#'
#'
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 "aes" "geom_line" "facet_wrap" "ylab" "xlab"
#'
#'
#' @export

plot_key_run <- function(dataSet,ylabel="change this label",is.aggregated=F,speciesNames,scales="fixed",speciesSuitability=1,predOrPrey="Prey") {

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
    if (is.aggregated) warning("Aggregregation has no efect this data set")
        # we plot each prey species and plot
      nSizeClasses <- dim(dataSet)[2]
      predSc <- data.frame(species = rep(c(1:nSpecies),each  = nSizeClasses),sizeClass=rep(1:nSizeClasses))
      addOnMAt <- apply(predSc,2,rep,nSpecies*nSizeClasses)
      # map first column (predator/size class to df)
      df <- cbind(predSc,df)

      if (tolower(predOrPrey) == "prey") {
        names(df) <- c("speciesNumber","predSizeClass","index","preySizeClass","preySpecies","dataF")
        df <- dplyr::inner_join(df,speciesNames,"speciesNumber")
        df <- dplyr::select(df,commonName,predSizeClass,preySizeClass,preySpecies,dataF)
        df$commonName <- factor(df$commonName, levels = unique(df$commonName))

        # see species suitability input. Facet plot
        df <- dplyr::filter(df,preySpecies==speciesSuitability)
        p <- ggplot2::ggplot(data=df) +
          geom_line(mapping =  aes(x = preySizeClass, y = dataF, group = predSizeClass)) +
          ylab(ylabel) + xlab(paste0("Prey (sizeClass) = ",speciesNames$commonName[speciesSuitability])) +
          facet_wrap( ~ commonName, scales = scales)
        print(p)
      } else if (tolower(predOrPrey) == "pred") {
        names(df) <- c("predSpecies","predSizeClass","index","preySizeClass","speciesNumber","dataF")
        df <- dplyr::inner_join(df,speciesNames,"speciesNumber")
        df <- dplyr::select(df,predSpecies,predSizeClass,preySizeClass,commonName,dataF)
        df$commonName <- factor(df$commonName, levels = unique(df$commonName))

        # see species suitability input. Facet plot
        df <- dplyr::filter(df,predSpecies==speciesSuitability)
        p <- ggplot2::ggplot(data=df) +
          geom_line(mapping =  aes(x = predSizeClass, y = dataF, group = preySizeClass)) +
          ylab(ylabel) + xlab(paste0("Predator (sizeClass) = ",speciesNames$commonName[speciesSuitability])) +
          facet_wrap( ~ commonName, scales = scales)
        print(p)

     }

  } else if (length(dimOfData) == 2) { # output by species over time (SSB & recruits) or M1)
    if (is.aggregated) warning("Aggregregation has no efect for 2D data sets")
    if (dimOfData[1] == dim(speciesNames)[1]) { # nSpeces x  nYears
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
# check to make sure species entered for suitability conforms to either an interger or the species name
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
