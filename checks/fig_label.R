# A label with the panel behind it masked, for labels that have to sit over a curve or a bar.
#
# WHY. label_ink_check.R found twenty-one labels across nine generators printed over by the
# figure's own ink: readable in the R window at 1300 pixels, unreadable in the built PDF at
# journal width. Moving each one is not always possible, because the place that carries the
# meaning is often the place the ink is. Masking the panel behind the text is, and it is what
# the rest of the release already does by hand in one figure.
#
#   source("checks/fig_label.R")
#   lab_on(x, y, "text", bg = "#fbfbfa", col = mark, cex = 0.8, adj = 0)
#
# bg must be the panel's own background so the patch is invisible. pad is a fraction of the
# string's height, added on every side.
lab_on <- function(x, y, labels, bg = "white", pad = 0.18, ...) {
  o <- list(...)
  cex <- if (!is.null(o$cex)) o$cex else 1
  adj <- if (!is.null(o$adj)) o$adj[1] else 0.5
  for (i in seq_along(labels)) {
    xi <- x[[min(i, length(x))]]; yi <- y[[min(i, length(y))]]
    w <- abs(strwidth(labels[i], cex = cex)); h <- abs(strheight(labels[i], cex = cex))
    x0 <- xi - adj*w
    graphics::rect(min(x0, x0 + w) - pad*h, yi - (0.5 + pad)*h,
                   max(x0, x0 + w) + pad*h, yi + (0.5 + pad)*h,
                   col = bg, border = NA)
  }
  graphics::text(x, y, labels, ...)
}
