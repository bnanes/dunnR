#' Dunn test for multiple comparisons using rank sums
#'
#' Non-parametric testing procedure for identifying differences between multiple
#' groups in a population. Implementation of Dunn JO (1964). Multiple
#' Comparisons Using Rank Sums, Technometrics, 6:3, 241-252. doi: (10.1080/00401706.1964.10490181)[https://doi.org/10.1080/00401706.1964.10490181]
#'
#' @param x a numeric vector of data values, or a list of numeric data vectors.
#'   Non-numeric elements of a list will be coerced, with a warning.
#' @param g a vector or factor object giving the group for the corresponding
#'   elements of x. Ignored with a warning if x is a list.
#' @param formula a formula of the form response ~ group where response gives
#'   the data values and group a vector or factor of the corresponding groups.
#' @param compare.all if `TRUE` (default), will consider all possible
#'   between-group comparisons
#' @param baseline.group if `compare.all = FALSE`, will compare this group to
#'   all others.
#' @param contrast if `compare.all = FALSE` and `baseline.group = NULL`, matrix
#'   specifying groups to compare. See Details.
#' @param p.adj.method method for adjusting P values to account for multiple
#'   comparisons. Default is the Bonferroni correction, following the original
#'   manuscript. See [p.adjust()].
#' @param alpha level at which the null hypothesis is rejected. See Details.
#'
#' @return Object of class `dunn.htest`. If the ensemble test is passed, the
#'   result will also have class `pairwise.htest`. If the ensemble test is not
#'   passed and the between-group comparisons are not run, the result will also
#'   have class `htest`.
#'
#'   This design provides for a specially-formatted `print` method ([print.dunn.htest()]), as
#'   well as downstream compatibility with methods expecting `htest` or `pairwise.htest`
#'   methods as input (such as [broom::tidy()]). It also holds additional information
#'   not provided for by `pairwise.htest`.
#'
#'   `dunn.htest` is a list containing the following elements:
#'   * `method` Name of the procedure
#'   * `data.name` Input data
#'   * `p.value` Matrix of adjusted P values. `NxN`, where `N` is the number of groups. Or a single P value if only the ensemble test is run.
#'   * `p.adjust.method` Adjustment method for multiple comparisons
#'   * `p.unadjust` Matrix of unadjusted P values
#'   * `statistic` Matrix of Z scores
#'   * `statistic.name` By analogy to the `htest` class
#'   * `contrast.table` `data.frame` containing the above information in a potentially
#'     more useful format, with each row reporesenting one between-group comparison.
#'   * `alpha` Level for rejecting the null hypothesis
#'   * `pass.pre.test` `TRUE` if the ensemble test was passed
#'   * `pre.test` An `htest` object containing the results for the Kruskal-Wallis ensemble test.
#'     Included only if the ensemble test is passed. If the ensemble test is not passed it is returned directly.
#'
#' @details This package implements Dunn's procedure for identifying differences
#'   between multiple groups within a population. It is a two-step procedure,
#'   first deciding if all groups come from identically-distributed populations,
#'   then, if not, proceeding to post-hoc comparisons between individual groups.
#'   The procedure compares the average ranks of values in each group. The
#'   first-step ensemble test is the Kruskal-Wallis rank sum test
#'   ([kruskal.test()]). The post-hoc tests are similar to Wilcoxan rank sum
#'   tests ([pairwise.wilcox.test()]), but the ranks of values within the larger sample
#'   are maintained. The procedure includes an adjustment for ties.
#'
#'   Essentially, the Dunn test is a non-parametric analog of a one-way ANOVA
#'   ensemble test followed by post-hoc tests such as Tukey. It has become
#'   commonly used in biomedical research, likely because it is included in
#'   popular statistics and graphing software packages (ex.:
#'   [https://www.graphpad.com/guides/prism/latest/statistics/how_the_kruskal-wallis_test_works.htm]).
#'
#'   In principle, any subset of possible between-group comparisons may be
#'   considered. By default, all possible comparisons are made. Alternatively, a
#'   baseline group may be compared to all others by setting `compare.all =
#'   FALSE` and indicating the `baseline.group`. To specify
#'   an arbitrary set of comparisons, also set `baseline.group = NULL` and
#'   identify desired between-group comparisons as an `NxN` `logical` matrix,
#'   where `N` is the number of groups. Matrix row and column names should be
#'   set to the group levels. A full set of comparisons is indicated by the
#'   lower triangle ([lower.tri()]); the upper triangle gives the same result,
#'   but ordered differently. Note that comparing a group to itself (diagonal
#'   elements) or double-counting comparisons will adversely affect the results.
#'
#'   Dunn applies the Bonferroni correction to control the family-wise error
#'   rate, and that is the default here. This approach can be overridden using
#'   the `p.adj.method` parameter to specify any of the methods supported by
#'   [p.adjust()].
#'
#'   The post-hoc between-group inferences depend on the determination of the
#'   preceding ensemble test. Therefore, if a Kruskal-Wallis fails to reject the
#'   null hypothesis that all groups are drawn from identically-distributed
#'   populations (`P >= alpha`), between-group comparisons will not run, and
#'   only the ensemble test is reported. To override this behavior and force the
#'   post-hoc tests to run, increase `alpha` above its default value of .05.
#'
#' @examples
#' dunn.test(mpg ~ cyl, data = mtcars)
#'
#' @references
#' Dunn JO (1964). Multiple Comparisons Using Rank Sums, Technometrics, 6:3, 241-252.  doi: (10.1080/00401706.1964.10490181)[https://doi.org/10.1080/00401706.1964.10490181]
#'
#' @export
dunn.test <- function(x, ...) UseMethod("dunn.test")

#' @rdname dunn.test
#' @export
dunn.test.default <- function(x, g,
                        compare.all = TRUE, baseline.group = ifelse(missing(g), NULL, levels(factor(g))[1]), contrast = NULL,
                        p.adj.method = p.adjust.methods, alpha = 0.05, ...) {
    if (is.list(x)) {
        if (length(x) < 2L)
            stop("'x' must be a list with at least 2 elements")
        if (!missing(g))
            warning("'x' is a list, so ignoring argument 'g'")
        dataName <- deparse(substitute(x))
        if (!all(sapply(x, is.numeric)))
            warning("some elements of 'x' are not numeric and will be coerced to numeric")
        k <- length(x)
        l <- sapply(x, "length")
        if (any(l == 0L))
            stop("all groups must contain data")
        g <- factor(rep.int(seq_len(k), l))
        x <- unlist(x)
    } else {
        if (length(x) != length(g))
            stop("'x' and 'g' must have the same length")
        if (!is.factor(g))
            g <- factor(g)
        if (!all(is.finite(g)))
            stop("all group levels must be finite")
        dataName <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    }

    # Default to Bonferroni, which is the method used in the original publication
    if(length(p.adj.method>1) || is.null(p.adj.method)) p.adj.method <- "bonferroni"
    p.adj.method <- match.arg(p.adj.method)

    # Default to all between-group comparisons
    # Then compare baseline.group to all others
    # Then use the provided contrast matrix
    if(compare.all)
        contrast = lower.tri(matrix(nrow=nlevels(g),ncol=nlevels(g)))
    else if(!is.null(baseline.group)) {
        contrast = matrix(FALSE, nrow=nlevels(g),ncol=nlevels(g))
        contrast[,baseline.group] <- TRUE
        contrast[baseline.group,baseline.group] <- FALSE
    }

    .dunn.test(x, g, dataName, contrast = contrast, p.adj.method=p.adj.method, alpha=alpha)
}

#' @rdname dunn.test
#' @export
dunn.test.formula <- function(formula, data, ...) {
    if(missing(formula) || (length(formula) != 3L))
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m[[1L]] <- quote(stats::model.frame)
    mf <- model.frame(formula, data)
    if(length(mf) > 2L)
        stop("'formula' should be of the form response ~ group")
    dataName <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    y <- dunn.test(unlist(mf[1]), unlist(mf[2]), ...)
    y$data.name <- dataName
    y
}

.dunn.test <- function(values, groups, dataName, contrast=NULL, p.adj.method="bonferroni", alpha=0.05, ...) {
    # values: vector with values to test
    # groups: indicator of the group to which each value belongs
    # contrast: square logical matrix of group comparisons to test, if null defaults to lower triangle (all comparisons)
    # p.adj.method: multiple comparison correction

    if(!is.factor(groups))
        groups <- factor(groups)
    ngroups = nlevels(groups)
    if(ngroups<2)
        stop("Only one group")

    kwPretest <- stats::kruskal.test(values, groups)
    if(kwPretest$p.value > alpha) {
        kwPretest$alpha <- alpha
        kwPretest$pass.pre.test <- FALSE
        class(kwPretest) <- c("dunn.htest", "htest")
        return(kwPretest)
    }

    if(is.null(contrast))
        contrast = lower.tri(matrix(nrow=nlevels(groups),ncol=nlevels(groups)))

    nas <- is.na(values)
    values <- values[!nas]
    groups <- groups[!nas]
    dups <- factor(values[duplicated(values)])
    dupGroupLengths <- sapply(levels(dups), function(l) sum(dups==l)+1)
    tieCor <- ifelse(
        length(dupGroupLengths)>0,
        sum(sapply(dupGroupLengths, function(t) (t^3)-t)),
        0)
    N <- length(values)

    values <- rank(values, ties.method="average")

    Qstat <- function(x, y) { # x and y are sets of *ranked* data points to compare
        abs(mean(x)-mean(y)) /
            sqrt(
                ((N*(N+1)/12) - (tieCor/(12*(N-1)))) *
                    ((1/length(x))+(1/length(y)))
            )
    }

    if(nrow(contrast) != ncol(contrast) || nrow(contrast) != ngroups)
        stop("`contrast` must be a square matrix with width and height equal to the number of groups")

    compInd <- which(contrast>0, arr.ind=TRUE)
    compTab <- NULL
    Q <- matrix(nrow=ngroups, ncol=ngroups)
    P.unadj <- matrix(nrow=ngroups, ncol=ngroups)
    for(i in 1:nrow(compInd)) {
        xi <- values[groups==levels(groups)[compInd[i,2]]]
        xj <- values[groups==levels(groups)[compInd[i,1]]]
        Qthis <- Qstat(xi,xj)
        Q[compInd[i,1], compInd[i,2]] <- Qthis
        P.unadj[compInd[i,1], compInd[i,2]] <- pnorm(Qthis, lower.tail=FALSE)*2
        compTab <- rbind(compTab, data.frame( # Potentially more convient data frame output
            Group.A=levels(groups)[compInd[i,2]],
            Group.B=levels(groups)[compInd[i,1]],
            Z=Qthis,
            p.unadjust=P.unadj[compInd[i,1], compInd[i,2]]))
    }
    P.adj <- p.adjust(P.unadj, method=p.adj.method)
    P.adj <- matrix(P.adj,nrow=nrow(P.unadj),ncol=ncol(P.unadj))
    compTab$p.value <- p.adjust(compTab$p.unadjust, method=p.adj.method)

    rObj <- list(
        method = "Dunn test for multiple comparisons using rank sums",
        data.name = dataName,
        p.value = P.adj,
        p.adjust.method = ifelse(p.adj.method=="bonferroni", "classic Dunn test (Bonferroni)", p.adj.method),
        p.unadjust = P.unadj,
        statistic = Q
    )
    dimnames(rObj$p.value) <- list(levels(groups), levels(groups))
    dimnames(rObj$p.unadjust) <- list(levels(groups), levels(groups))
    dimnames(rObj$statistic) <- list(levels(groups), levels(groups))
    rObj$statistic.name <- "Z"
    rObj$contrast.table <- compTab
    rObj$alpha <- alpha
    rObj$pass.pre.test <- TRUE
    rObj$pre.test <- kwPretest

    if(sum(!is.na(rObj$p.value[1,])) == 0 && sum(!is.na(rObj$p.value[,ncol(rObj$p.value)])) == 0) { # To try to match other `pairwise.htest` functions, trim the p.value matrices if using a lower triangle contrast matrix
        rObj$p.value <- rObj$p.value[-1,-1*ncol(rObj$p.value)]
        rObj$p.unadjust <- rObj$p.unadjust[-1,-1*ncol(rObj$p.unadjust)]
    }

    class(rObj) <- c("dunn.htest", "pairwise.htest")

    rObj
}
