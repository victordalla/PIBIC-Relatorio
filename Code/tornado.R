## Tornado Diagram
tornado <- function(lb, ub, base, lowlabels, uplabels, titleplot="Tornado Diagram",
                    xlim=c(min(lb), max(ub)), mval=11, dsa.var=dsavar, ncol=1, space=.2){
  ## This function plots a tornad diagram used to present the results of deterministic sensitivity analysis
  ## parameteres:
  ## lb and ub = lower and upper ICER, respectively
  ## base = ICER for the base case scenario
  ## lowlabels and uplabels = labels to be plotted in the lower end upper end of the bars
  ## mval = variable that controls de width of the left margin, so that you can fit the names of the variables
  ## dsa.var = names of the variables
  ## ncol = from 1 to 8 to select a color for the bars
  ## space = space between bars

  library(RColorBrewer)
  ## mycol <- brewer.pal(8,"Dark2")[ncol]
  ## mycol <- gray(.7) ## color for grey scale plots
  mycol <- rgb(110, 141, 130, max=255) ## Same color used in the model diagram

  range <- ub-lb
  ix <- sort(range, index.return=TRUE)$ix

  par(mar=c(5, mval, 4, 2) + .1, mgp=c(2, .5, 0))
  bp1 <- barplot(lb[ix] - base, main=titleplot, space=space, horiz=TRUE, names.arg=dsa.var[ix],
                 cex.names=1, xlab=quote(bold("ICER in thousands USD")), offset=base, xlim=xlim,
                 col=mycol, las=1, xaxt="n")
  abline(v=base, col=grey(.5))
  bp2 <- barplot(ub[ix] - base, space=space, horiz=TRUE, add=TRUE, col=mycol, offset=base, las=1, xaxt="n")

  ticks <- seq(0, xlim[2], by=20000)
  axis(1, at=ticks, labels=as.character(ticks/1000), col=grey(.5))
  text(lb[ix], bp1, lowlabels[ix], cex=.9, pos=2)
  text(ub[ix], bp2, uplabels[ix], cex=.9, pos=4)
  ## return(bp1)
}
