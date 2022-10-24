#' Print method for Dunn test results
#'
#' @param x a `dunn.htest` object
#' @param digits significant digits
#' @param ...
#'
#' @return `x`, silently
#' @export
print.dunn.htest <- function(x, digits = getOption("digits"), prefix = "\t", ...) {

    cat("\n")
    cat(strwrap(x$method, prefix = prefix), sep = "\n")
    cat("\n")
    cat("data:  ", x$data.name, "\n", sep = "")

    preTestLine <- function(x) {
        out <- character()
        if(!is.null(x$statistic))
            out <- c(out, paste(names(x$statistic), "=",
                                format(signif(x$statistic, max(1L, digits - 2L)))))
        if(!is.null(x$parameter))
            out <- c(out, paste(names(x$parameter), "=",
                                format(signif(x$parameter, max(1L, digits - 2L)))))
        if(!is.null(x$p.value)) {
            fp <- format.pval(x$p.value, digits = max(1L, digits - 3L))
            out <- c(out, paste("p-value",
                                if(substr(fp, 1L, 1L) == "<") fp else paste("=",fp)))
        }
        cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
    }

    if(!x$pass.pre.test) {
        preTestLine(x)
        cat(paste0("Because P >= ", signif(x$alpha, max(1L, digits - 3L)), ", the post-hoc test was not run.\n"))
        return(invisible(x))
    }

    preTestLine(x$pre.test)

    sigL0 <- ifelse(x$alpha > 0.05, x$alpha, NA)
    sigL1 <- min(x$alpha, 0.05)
    sigL2 <- ifelse(sigL1<=0.01, ifelse(sigL1<=0.001, NA, 0.001), 0.01)
    sigL3 <- ifelse((sigL1<=0.001)||(sigL2<=0.001), NA, 0.001)
    getSignifCode <- function(p) {
        if(!is.na(sigL3) && p <= sigL3)
            return("***")
        if(!is.na(sigL2) && p <= sigL2)
            return("**")
        if(!is.na(sigL1) && p <= sigL1)
            return("*")
        if(!is.na(sigL0) && p <= sigL0)
            return(".")
        return("")
    }

    cat("\n-----------------------------------------")
    cat("\nPost-hoc comparisons between groups")
    cat(paste0("\nP value adjustment method: ", x$p.adjust.method))
    cat("\n-----------------------------------------\n")
    tabForm <- x$contrast.table
    tabForm$Z <- format(signif(tabForm$Z, max(1L, digits - 2L)))
    tabForm$p.unadjust <- format.pval(tabForm$p.unadjust, digits = max(1L, digits - 3L))
    tabForm$p.value <- format.pval(tabForm$p.value, digits = max(1L, digits - 3L))
    tabForm$Sig. <- sapply(x$contrast.table$p.value, getSignifCode)
    print(tabForm, quote = FALSE)
    cat("-----------------------------------------\n")

    if(!is.na(sigL0)) cat("., P<", sigL0, "; ", sep = "")
    if(!is.na(sigL1)) cat("*, P<", sigL1, ifelse(is.na(sigL2), "", "; "), sep = "")
    if(!is.na(sigL2)) cat("**, P<", sigL2, ifelse(is.na(sigL3), "", "; "), sep = "")
    if(!is.na(sigL3)) cat("***, P<", sigL3, sep = "")
    cat("\n")

    invisible(x)
}
