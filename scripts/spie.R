spiechart <- function(expected, observed, ...) {
# Plot a Spie chart.

   if (any(expected == 0)) {
      warning("0 in 'expected': cannot plot spie chart")
      return
   }

   # Calls to 'tapply()', which a relatively frequent in this
   # context return an 'array', causing trouble downstream.
   if (is.array(expected)) expected <- as.vector(expected)
   if (is.array(observed)) observed <- as.vector(observed)

   # A bit of maths.
   angle <- 2*pi * diffinv(expected / sum(expected))
   area <- pi * observed / sum(observed)
   radius <- sqrt(2 * area / diff(angle))
   R <- max(radius)

   # Dispatch dot arguments.
   dotargs = list(...)
   valid_plot_args <- c("main")
   plotargs <- dotargs[names(dotargs) %in% valid_plot_args]
   dotargs <- dotargs[!names(dotargs) %in% valid_plot_args]
   plotargs[["type"]] <- "n"
   plotargs[["bty"]] <- "n"
   plotargs[["xaxt"]] <- "n"
   plotargs[["yaxt"]] <- "n"
   plotargs[["xlab"]] <- ""
   plotargs[["ylab"]] <- ""
   plotargs[["x"]] <- c(-R,R)
   plotargs[["y"]] <- c(-R,R)

   # Do the plot.
   do.call(what=plot, args=plotargs);

   # Manually recycle the parameters passed to 'polygon()'.
   # NB: 'polargs.i' is a list of lists of arguments.
   polargs <- rep(
      do.call(
         what = mapply,
         args = c(list(FUN=list, SIMPLIFY=FALSE), dotargs)
      ),
      length.out = length(radius)
   )

   for (i in 1:length(radius)) {
      nsteps <- (angle[i+1] - angle[i]) / .01
      steps <- seq(from=angle[i], to=angle[i+1], length.out=nsteps)
      x <- c(0, radius[i]*cos(steps))
      y <- c(0, radius[i]*sin(steps))
      do.call(what=polygon, args=c(list(x=x, y=y), polargs[[i]]))
   }
   steps <- seq(0, 2*pi, length.out=628)
   lines(cos(steps), sin(steps), col="grey50")
}
