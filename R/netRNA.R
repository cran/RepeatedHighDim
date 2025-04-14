#' @title netRNA:Network meta-analysis for gene expression data
#'
#' @description This function conducts network meta-analysis using gene expression data to make indirect comparisons between different groups. It computes the p values for each gene and the fold changes, and provides a dataframe  containing these results.
#'
#' @param TE A list containing log fold changes from two individual studies.
#' Index names of the list should be the gene names; otherwise,
#' each value of the 'name' column in the output dataframe will correspond to
#' the position in the list, rather than gene identifiers.
#' @param seTE A list containing standard errors of log fold changes from two individual studies.
#' @param treat1 A vector with Label/Number for first treatment.
#' @param treat2 A vector with Label/Number for second treatment.
#' @param studlab A vector containing study labels
#'
#'
#' @details
#' The function supports a simple network with three nodes, where one node represents
#' a control group and the two other nodes represent treatment (or diseased) groups.
#' While the user provides fold changes and their standard errors of each treatment versus control as input,
#' the function calculates the fold changes for the indirect comparison between the two treatments.
#'
#' @return A list containing the p values for each gene, the fold changes,
#' the upper and lower bounds for the 95\% CI of the log fold changes, and a summary dataframe with results for each gene.
#'
#' @author Klaus Jung, Sergej Ruff
#' @references
#'
#' Winter, C., Kosch, R., Ludlow, M. et al. Network meta-analysis correlates with analysis of merged independent transcriptome expression data.
#' \emph{BMC Bioinformatics} \strong{20}, 144 (2019). \doi{https://doi.org/10.1186/s12859-019-2705-9}
#'
#' Rücker G. (2012). Network meta-analysis, electrical networks and graph theory.
#' \emph{Research synthesis methods}, \strong{3(4)}, 312–324. \doi{https://doi.org/10.1002/jrsm.1058}
#'
#' @importFrom netmeta netmeta
#' @importFrom progress progress_bar
#' @export
#' @examples
#'
#'\dontrun{
#' #'#######################
#'### Data generation ###
#'#######################
#'n = 100 ### Sample size per group
#'G = 100 ### Number of genes
#'
#'### Basic expression, fold change, batch effects and error
#'alpha.1 = rnorm(G, 0, 1)
#'alpha.2 = rnorm(G, 0, 1)
#'beta.1 = rnorm(G, 0, 1)
#'beta.2 = rnorm(G, 0, 1)
#'gamma.1 = rnorm(G, 0, 1)
#'gamma.2 = rnorm(G, 2, 1)
#'delta.1 = sqrt(invgamma::rinvgamma(G, 1, 1))
#'delta.2 = sqrt(invgamma::rinvgamma(G, 1, 2))
#'sigma.g = rep(1, G)
#'
#'# Generate gene names
#'gene_names <- paste("Gene", 1:G, sep = "")
#'
#'### Data matrices of control and treatment (disease) groups
#'C.1 = matrix(NA, G, n)
#'C.2 = matrix(NA, G, n)
#'T.1 = matrix(NA, G, n)
#'T.2 = matrix(NA, G, n)
#'
#'for (j in 1:n) {
#'  C.1[,j] = alpha.1 + (0 * beta.1) + gamma.1 + (delta.1 * rnorm(1, 0, sigma.g))
#'  C.2[,j] = alpha.1 + (0 * beta.2) + gamma.2 + (delta.2 * rnorm(1, 0, sigma.g))
#'  T.1[,j] = alpha.2 + (1 * beta.1) + gamma.1 + (delta.1 * rnorm(1, 0, sigma.g))
#'  T.2[,j] = alpha.2 + (1 * beta.2) + gamma.2 + (delta.2 * rnorm(1, 0, sigma.g))
#'}
#'
#'study1 = cbind(C.1, T.1)
#'study2 = cbind(C.2, T.2)
#'
#'# Assign gene names to row names
#'#rownames(study1) <- gene_names
#'#rownames(study2) <- gene_names
#'#############################
#'### Differential Analysis ###
#'#############################
#'
#'if(check_limma()){
#'### study1: treatment A versus control
#'group = gl(2, n)
#'M = model.matrix(~ group)
#'fit = limma::lmFit(study1, M)
#'fit = limma::eBayes(fit)
#'p.S1 = fit$p.value[,2]
#'fc.S1 = fit$coefficients[,2]
#'fce.S1 = sqrt(fit$s2.post) * sqrt(fit$cov.coefficients[2,2])
#'
#'### study2: treatment B versus control
#'group = gl(2, n)
#'M = model.matrix(~ group)
#'fit = limma::lmFit(study2, M)
#'fit = limma::eBayes(fit)
#'p.S2 = fit$p.value[,2]
#'fc.S2 = fit$coefficients[,2]
#'fce.S2 = sqrt(fit$s2.post) * sqrt(fit$cov.coefficients[2,2])
#'
#'
#'
#'#############################
#'### Network meta-analysis ###
#'#############################
#'p.net = rep(NA, G)
#'fc.net = rep(NA, G)
#'treat1 = c("uninfected", "uninfected")
#'treat2 = c("ZIKA", "HSV1")
#'studlab = c("experiment1", "experiment2")
#'fc.true = beta.2 - beta.1
#'
#'TEs <- list(fc.S1, fc.S2)
#'seTEs <- list(fce.S1, fce.S2)
#'}
#'
#'# Example usage:
#'test <- netRNA(TE = TEs, seTE = seTEs, treat1 = treat1, treat2 = treat2, studlab = studlab)
#'}
#'@seealso
#'For more information, please refer to the package's documentation and the tutorial: \url{https://software.klausjung-lab.de/}.
netRNA <- function(TE, seTE, treat1, treat2, studlab) {




  # Check if TE is a list
  if (!is.list(TE)) {
    stop("Error: TE should be a list containing log fold changes.Currently TE is not a list")
  }
  # Check if seTE is a list
  if (!is.list(seTE)) {
    stop("Error: seTE should be a list containing standard errors of treatment effects.Currently seTE is not a list")
  }
  # Check if there are exactly three unique treatments
  if (length(unique(c(treat1, treat2))) != 3) {
    stop("Error: netRNA currently only supports the comparison of two distinct treatments.\n",
         "Ensure that treat1 contains only one common control group with the same label,\n",
         "and treat2 contains two distinct labels for the two groups intended for comparison.")
  }



  # Check log-fold values
  if (length(TE) == 0) {
    stop("Error: TE has zero values.")
  }
  # Check if all elements in TE are numeric
  if (!all(sapply(TE, is.numeric))) {
    stop("Error: Not all elements in TE are numeric. Make sure they are numeric.")
  }

  # Check standard error values
  if (length(seTE) == 0) {
    stop("Error: seTE has zero values.")
  }
  # Check if all elements in seTE are numeric
  if (!all(sapply(seTE, is.numeric))) {
    stop("Error: Not all elements in seTE are numeric.Make sure they are numeric.")
  }
  # Check if lengths of seTE and TE are the same
  if (length(seTE) != length(TE)) {
    stop("Error: Lengths of seTE and TE are not the same. Please ensure that each log fold change has a corresponding standard error.")
  }
  # Check treatment groups for emptiness
  if (length(treat1) == 0 || length(treat2) == 0) {
    stop("Error: The treatment groups cannot be empty. Please provide both control groups (treat1) and two groups intended for comparison (treat2).")
  }
  if (any(treat1 == treat2)) {
    stop("Error: Treatments must be different (arguments 'treat1' and 'treat2'). Make sure that elements do not repeat between treat1 and treat2")

  }




  G <- min(lengths(c(TE, seTE)))

  p.net <- vector("numeric", length = G)
  fc.net <- vector("numeric", length = G)
  fc.upper <- vector("numeric", length = G)
  fc.lower <- vector("numeric", length = G)
  names_list <- vector("list", length = G)  # Initialize a list to store names



  pb <- progress_bar$new(total = G, format = "[:current/:total] genes completed [:bar] :elapsed")

  for (j in 1:G) {
    pb$tick()

    # Extracting names and storing in the list
    if (is.null(names(TE[[1]][j]))) {
      names_list[[j]] <- j
    } else {
      names_list[[j]] <- names(TE[[1]][j])
    }


    TE_j <- unlist(lapply(TE, "[[", j))  # Fold changes from individual studies
    seTE_j <- unlist(lapply(seTE, "[[", j))  # Error of fold changes from individual studies

    exampledata <- data.frame(TE = TE_j, seTE = seTE_j, treat1, treat2, studlab)
    net1 <- netmeta(TE_j, seTE_j, treat1, treat2, studlab, data = exampledata, sm = "MD")
    S <- summary(net1)

    p.net[[j]] <- S$fix$p[3,1]
    fc.net[[j]] <- S$fixed$TE[3,1]
    fc.upper[[j]] <- S$fix$upper[3,1]
    fc.lower[[j]] <- S$fix$lower[3,1]

  }
  pb$terminate()

  # Create dataframe
  df <- data.frame(name = unlist(names_list),p.net = p.net, fc.net = fc.net,
                   fc.lower=fc.lower,fc.upper=fc.upper)



  return(list(p.net=p.net,fc.net=fc.net,fc.lower=fc.lower,fc.upper=fc.upper
              ,df=df))


}



