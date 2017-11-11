.onLoad <- function(...) {
  packageStartupMessage("Note: this is a pre-release development version of microsynth. Please send any observed bugs, functionality requests, or comments to the developers at mrobbins@rand.org or sdavenpo@rand.org.")
}

makemicrosynth <- function(x) {
  class(x) <- c("microsynth")
  x
}
