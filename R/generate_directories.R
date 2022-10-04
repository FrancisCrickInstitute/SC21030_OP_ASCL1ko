#' Generate directory if it does not exist
#'
#' @param dir_path path to the new directory
#'
#' @return
#' @export
#'
#' @examples
#' 
#' generate_directories("./new/dir")

generate_directories <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}