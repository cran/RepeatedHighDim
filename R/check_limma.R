#' Check for package dependency
#'
#' @title Check for 'limma' availability
#' @description checks if the 'limma' package is installed. If not already installed,
#' limma will be installed automatically.
#' @author Sergej Ruff
#' @importFrom utils install.packages menu
#' @export
#' @keywords internal
check_limma <- function() # Returns TRUE if available, FALSE otherwise
{
  if(requireNamespace("limma", quietly=TRUE)) return(TRUE)
  if(!interactive()) return(FALSE)
  inst <- menu(c("Yes", "No"), title="Package {limma} required but not installed.\nDo you want to install it now?")
  if(inst != 1)
  {
    message("To run this example, first install {limma} following the directions at 'https://bioconductor.org/packages/limma'")
    return(FALSE)
  }
  # the following could be wrapped in try and conditionally return TRUE / FALSE
  if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", quiet=TRUE)
  BiocManager::install("limma", update=FALSE, ask=FALSE, quiet=TRUE)
  return(TRUE)
}
