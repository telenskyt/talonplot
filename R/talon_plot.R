
#' The 'talon plot': a plot of contributions of survival and recruitment to population growth (lambda)
#'
#' @description
#' Makes a so called "talon plot", a plot showing contributions of adult survival and recruitment
#' to population growth. It is a 2D plot with three axes, population growth being a sum of adult
#' survival and recruitment. It can reveal important patterns in population dynamics,
#' see Telenský et al. (2023) <doi:>.
#'
#' @param recr_pts,surv_pts vector of mean values of recruitment and survival, respectively,
#'	for each data point (data point corresponding e.g. to a temporal occasion)
#' @param recr_samples,surv_samples optional. Matrices - rows correspond to MCMC samples,
#'	columns to the same data points as the *_pts variables. Used for plotting and/or calculating 95\% confidence ellipses.
#' @param recr_CI_low,recr_CI_high,surv_CI_low,surv_CI_high 95\% confidence limits for *_pts variables
#' @param CI_type "ellipse", "cross", or "none"
#' @param CI_transform only for ellipses now, says whether they should be for bivariate normal
#'	on the scale [recr,surv] for CI_transform = FALSE,
#'	or on the scale [log(recr),logit(surv)] for CI_transform = TRUE
#' @param ... other graphical parameters
#'
#' @return
#' @export
#' @references Telenský, T., Storch, D., Klvaňa, P., Reif, J. (2023). Extension of Pradel capture-recapture survival-recruitment model accounting for transients. Methods in Ecology and Evolution. In press.
#'
#' @examples
#' @importFrom car dataEllipse
#' @importFrom grDevices col2rgb dev.cur dev.new dev.off dev.set rgb dev.size
#' @importFrom graphics arrows lines par points segments text
#' @importFrom gridGraphics grid.echo
#' @importFrom grid grid.grab grid.newpage grid.text unit gpar pushViewport grid.draw viewport
talon_plot <- function (recr_pts, surv_pts, recr_samples = NULL, surv_samples = NULL, recr_CI_low, recr_CI_high, surv_CI_low, surv_CI_high, col.pts = "black", col.samples = alpha("lightgrey", 128),
	plot.samples = TRUE,
	CI_type = c("ellipse", "cross", "none"),
	CI_transform = FALSE,
	col.CI = "grey",
	lim_include_zero_one = TRUE, main, #labels = NULL,
	big_font,
	ylab_offset = if (big_font) 0.25 else 0.15,
	cex = if (big_font) 1.5 else 1,
	cex.pts = 1,  # cex.pts = cex
	cex.main = 1.2 * cex,
	cex.lab = cex,
	cex.axis = cex,
	cex.samples = 1,
	font.main = 2,
	mar = if (big_font) c(5.1, 5.1, 5.1, 5.1) else c(4.1, 4.1, 4.1, 4.1),
	... # other graphical parameters
	)
{
	CI_type <- match.arg(CI_type)
	if (lim_include_zero_one) {
		if (CI_type == "none") {
			xylim <- c(0, 1.05*max(1, recr_pts, surv_pts))
		} else {
			xylim <- c(0, 1.05*max(1, recr_CI_high, surv_CI_high))
		}
	} else {
		if (CI_type == "none") {
			xylim <- expand_range(range(recr_pts, surv_pts), 0.05)
		} else {
			xylim <- expand_range(range(recr_CI_high, recr_CI_low, surv_CI_high, surv_CI_low), 0.05)
		}
	}
	stopifnot(length(recr_pts) == length(surv_pts))
	n.occasions <- length(recr_pts)
	dev <- dev.cur()
	dev.new(width = 7, height = 7, noRStudioGD = TRUE, ...) # create auxilliary device, for the plot that will be rotated 45° degrees
	old.mar <- par("mar")
	par(mar = mar)
	old.cex.axis <- par("cex.axis")
	par(cex.axis = cex.axis)
	plot(recr_pts, surv_pts, col = col.pts, pch = 16, xlim = xylim, ylim = xylim,
		ylab = "", xlab = "recruitment", bty = "l", las = 1, yaxs="i", xaxs="i",
		cex.main = cex.main, cex.lab = cex.lab, cex = cex.pts) # asp = 1 doesn't work, it breaks the axes
	bbox <- par("usr")


	text(bbox[1] + -ylab_offset*(bbox[2]-bbox[1]), mean(bbox[3:4]), "adult survival", srt = -90, xpd = TRUE, pos = 1, cex = cex.lab)
	arrows(bbox[1], bbox[3], bbox[2], bbox[4], length = 0.15) # plot the diagonal axis
	arrows(bbox[1], bbox[3], bbox[1], bbox[4], length = 0.15, xpd = TRUE)
	arrows(bbox[1], bbox[3], bbox[2], bbox[3], length = 0.15, xpd = TRUE)


	if (plot.samples)
		points(recr_samples, surv_samples, pch = ".", col = col.samples, cex = cex.samples)
	if (CI_type == "ellipse") { # confidence ellipses!
		library(car)
		for (i in 1:(n.occasions - 1)) {
			if (!CI_transform) { # bivariate normal on the scale [recr,surv]:
				ell <- dataEllipse(x = recr_samples[,i], y = surv_samples[,i], levels = 0.95, draw = FALSE)
			}
			else { # bivariate normal on the scale [log(recr),logit(surv)]! And we scale it back! # but in the end it is less corresponding to normal distribution than the untransformed; so we prefer the untransformed
				# (see "Transformovat ci netransformovat CI elipsy v paratech.docx" for more detail)
				ell <- dataEllipse(x = log(recr_samples[,i]), y = logit(surv_samples[,i]), levels = 0.95, draw = FALSE)
				ell[,1] <- exp(ell[,1]) # transform back
				ell[,2] <- inv.logit(ell[,2])
			}
			lines(ell, col = col.CI, lwd = 1)
		}
	} else if (CI_type == "cross") {
		segments(recr_CI_low, surv_pts, recr_CI_high, surv_pts, col = col.CI)
		segments(recr_pts, surv_CI_low, recr_pts, surv_CI_high, col = col.CI)
	}
	points(recr_pts, surv_pts, pch = 16, cex = cex.pts, col = col.pts)

	#abline(1, -1, lty = 2, lwd = 2, col = "green")
	segments(1, 0, 0, 1, lty = 2, lwd = 1) # col = "green") # horizontal line for lambda = 1
	#text(bbox[1] + 0.9*(bbox[2]-bbox[1]), bbox[3] + 0.9*(bbox[4]-bbox[3]), " population  growth     ", srt = -45, xpd = TRUE, pos = 3, cex = cex)
	text(bbox[1] + 0.9*(bbox[2]-bbox[1]), bbox[3] + 0.9*(bbox[4]-bbox[3]), "population ", srt = -45, xpd = TRUE, adj = c(1, NA), cex = cex)
	text(bbox[1] + 0.9*(bbox[2]-bbox[1]), bbox[3] + 0.9*(bbox[4]-bbox[3]), " growth", srt = -45, xpd = TRUE, adj = c(0, NA), cex = cex)

	par(mar = old.mar)
	par(cex.axis = old.cex.axis)

	## !! NOW DO THE 45° DEGREE ROTATION!

	library(gridGraphics) # solution from https://stackoverflow.com/a/32198964/684229
	grid.echo() # redraws it without reason, no matter what I try (invisible(), x = NULL, newpage = FALSE)
	g <- grid.grab() # again redraws it without reason... even with invisible()


	dev.off() # close auxilliry device
	#if (dev != "null device")
	dev.set(dev) # need to come back to the device that was active before, because dev.off() could skip it (https://stackoverflow.com/q/77426629/684229)
	grid.newpage() #recording = FALSE)
	grid.text(label = main, x = unit(0.5, "npc"), y = unit(1, "npc") - unit(1, "line"), just = c("centre", "top"), gp = gpar(cex = cex.main, font = font.main)) # have to use grid.text() and not title(), see https://stackoverflow.com/q/77427520/684229
	pushViewport(viewport(width = unit(sqrt(2)/2*min(dev.size("in")), "in"), height = unit(sqrt(2)/2*min(dev.size("in")), "in"), angle=45, clip = "off"))
	grid.draw(g)
}

r.signif <- function (r.info) sign(r.info['2.5%']) == sign(r.info['97.5%'])

alpha<-function(colname,alpha){r<-col2rgb(colname);rgb(r[1],r[2],r[3],alpha,maxColorValue =255)}

# will expand given range by adding a given percentage of the range on each side
expand_range <- function(range, perc)
{
	(matrix(c(1+perc, -perc, -perc, 1+perc), nrow = 2) %*% range(range))[,1]
}

logit <- function (x) log(x/(1-x))
inv.logit <- plogis

