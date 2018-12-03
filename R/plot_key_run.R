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
#'
#'
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 "aes" "geom_line" "facet_wrap" "ylab"
#'
#' @export

plot_key_run <- function(dataSet,ylabel="ylabel",is.aggregated=T,speciesNames,scales="fixed") {

  dimOfData <- dim(dataSet)

  if (length(dimOfData) > 2) { # output with size structure, by species over time (catch, N, M2)
    df <- reshape::melt(dataSet)
    # transform data into a 2 array for use with ggplot
    names(df) <- c("sizeClass","speciesNumber","time","dataF")
    df <- dplyr::inner_join(df,speciesNames,"speciesNumber")
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


  } else if (length(dimOfData) == 2) { # output by species over time (SSB & recruits) or M1)
    if (dimOfData[1] == dim(speciesNames)[1]) { # nSpeces x  nYears
      if (is.aggregated) warning("Aggregregation has no efect for 2D data sets")
      df <- reshape::melt(dataSet)
      names(df) <- c("speciesNumber","year","dataF")
      df <- dplyr::inner_join(df,speciesNames,"speciesNumber")
      df <- dplyr::select(df,commonName,year,dataF)
      df$commonName <- factor(df$commonName, levels = unique(df$commonName))


      p <- ggplot2::ggplot(data=df) +
        geom_line(mapping =  aes(x = year, y = dataF)) + ylab(ylabel) +
        facet_wrap( ~ commonName, scales = scales)
      print(p)
    } else { # M1
      df <- reshape::melt(dataSet)
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
