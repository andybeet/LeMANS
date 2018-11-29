#' Plots output from the key_run script
#'
#'@param dataSet The data set to be plotted
#'@param title The title of the figure
#'
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 "aes" "geom_line" "facet_wrap" "ylab"
#'
#' @export

plot_key_run <- function(dataSet,title,is.aggregated=T) {

  dimOfData <- dim(dataSet)
  if (length(dimOfData) > 2) {
    # transform data into a 2 array for use with ggplot
    df <- reshape::melt(dataSet)
    names(df) <- c("sizeClass","species","time","data")

    newdf <- dplyr::group_by(df,species,time) %>% dplyr::summarize(aggData = sum(data,na.rm=T))

    p <- ggplot2::ggplot(data=newdf) +
      geom_line(mapping =  aes(x = time, y = aggData)) + ylab("Catch (units?)") +
      facet_wrap( ~ species)
    print(p)

  }



#return(newdf)

  if (length(dimOfData) == 3) {
    nSizes <- dim(dataSet)[1]
    nSpecies <- dim(dataSet)[2]
    nTimeSteps <- dim(dataSet)[3]
  }
  if (is.aggregated) { # aggregate over size class
    aggData2 <- apply(dataSet,c(2,3),sum,na.rm=T)
  }

  par(mfrow=c(4,6))
  par(mar=c(2,2,3,2)+0.1)
  par(oma=c(2,4,4,0))

  for (isp in 1:nSpecies) {
    plot(aggData2[isp,],main = title)
  }

  # plot(catchOverSize, col = "darkred", xlab = "Time", nc = 6,
  #      panel = function(x) { grid(col = "lightgray");lines(x) })


}
