## Renders man/figures/model-structure.png, the schematic README.md embeds, from its SVG with
## macOS's qlmanage. dev/make-doc-figures.R makes the other figures, and calls this when rsvg is
## not installed. Run from the package root:
##   Rscript dev/make-readme-figures.R

## ---- model-structure.png, from the SVG --------------------------------------------------
## Both documents embed the schematic; MODEL_DESCRIPTION uses the SVG and README the PNG, so
## the PNG has to be kept in step by hand unless it is generated. qlmanage fits an SVG to a
## SQUARE canvas and crops the overflow, which silently loses the right-hand panel - so the
## viewBox is padded to square first, then the whitespace is trimmed back off.
svg <- "man/figures/model-structure.svg"
png <- "man/figures/model-structure.png"
tmp <- file.path(tempdir(), "square.svg")
s <- readLines(svg, warn = FALSE)
s <- sub('viewBox="0 0 1600 1010" width="1600" height="1010"',
         'viewBox="0 -295 1600 1600" width="1600" height="1600"', s, fixed = TRUE)
writeLines(s, tmp)
out <- file.path(tempdir(), "render"); unlink(out, recursive = TRUE); dir.create(out)
system2("qlmanage", c("-t", "-s", "3000", "-o", shQuote(out), shQuote(tmp)),
        stdout = FALSE, stderr = FALSE)
rendered <- list.files(out, pattern = "[.]png$", full.names = TRUE)[1]
stopifnot(!is.na(rendered))

if (requireNamespace("png", quietly = TRUE)) {
  im <- png::readPNG(rendered)
  ## composite onto white, then trim the padding qlmanage added
  if (dim(im)[3] == 4) {
    a <- im[, , 4]
    for (k in 1:3) im[, , k] <- im[, , k] * a + (1 - a)
    im <- im[, , 1:3]
  }
  ink <- apply(im, c(1, 2), min) < 0.98
  rows <- which(apply(ink, 1, any)); cols <- which(apply(ink, 2, any))
  im <- im[min(rows):max(rows), min(cols):max(cols), ]
  png::writePNG(im, png)
  cat(sprintf("wrote %s (%d x %d, aspect %.2f; the SVG is 1.58)\n",
              png, dim(im)[2], dim(im)[1], dim(im)[2] / dim(im)[1]))
  if (abs(dim(im)[2] / dim(im)[1] - 1.58) > 0.2)
    warning("aspect is wrong - qlmanage probably cropped the width; check the render")
} else {
  stop("needs the 'png' package to trim the render: install.packages('png')")
}
