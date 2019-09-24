#' @title fit_models

#' @description fits several ecological models to test the effect of environmental variables on species traits across ecological gradients
#'

#' @param data The dataframe that contains the data
#' @param dep_var numeric. The name of the column in `data` that holds the trait for which to fit the predictive models
#' @param ind_var numeric. The name of the column in `data` that holds the numeric environmental variable (light, nutrients)
#' that will act as main predictor of the response variable `dep_var`
#' @param modelname name of the object that contains the scientific function to be fit. User must provide also a list containing
#' the initial values of the parameters of the model, the highest and the lowest possible values for those parameters. The names of these
#' lists must be exactly the same as `modelname` followed by the suffixes `.par`, `.par.hi` and `par.lo`, respectively.
#' @param species factor. Name of the species for which the model is gonna be fitted.
#'  If `species = NULL` then the model is fit for all the dataframe together.
#' @param size_var numeric. Name of the column in `data` that contains the size of plants. If `size = NULL` no effect of size is included in
#' the model fitted.
#' @param sites factor. Name of the column in `data` that contains the name of the different sites along the gradient.
#' This can be different elevations, different latitude, or different climate, for example. If `sites = NULL` the function will fit one model for
#' for all sites pooled together. Otherwise, the parameters of the model indicated in `rep_vectors` will be allowed to take a different value
#' for each site.
#' @param rep_vectors vector. In the parameter list (provided by the user), position of the parameters that will be allowed to vary across sites.
#' If sites != NULL and rep_vectors = NULL, all parameters in the function will be allowed to vary across sites.
#' @param max_iter maximum number of iterations allowed for the `anneal` function. By default `max_iter = 15000`.
#'

#'
#' @return This function fit a user-given scientific function to assess the effect of an environmental variable (`ind_var`) on
#' a given plant trait (`dep_var`) across an environmental gradient of a second variable (`sites`). The model uses the function ànneal`in the
#' `likelihood`package to fit the model, and returns a list containing:
#' - the best estimates of the parameters in the function
#' - the upper limits
#' - the lower limits
#' - LL
#' - AICc
#' - slope
#' -R2
#' - k
#'
#' (please see the documentation for the `anneal`function in the `likelihood` package for more info in the output of this function.

#'
#'
#'
#' @examples
#' data(ecophysio)
#' rdpi(ecophysio,sp,SB, Piso)
#'
#' @export


fit_models <- function (data, dep_var, ind_var, modelname, species = NULL, size_var = NULL,
                        sites = NULL, rep_vectors, max_iter = 15000) {

    # Create a working dataset for a species that has adequate sample sizes for all of those subzones
    if (!is.null(species)) {
        working_data <- dplyr::filter(data, Sp == species)  %>%
            mutate(Site = droplevels(Site))
    } else {
        working_data <- data
    }

    # drop any observations containing missing values for the independent or dependent variables
    working_data <- working_data %>%
        drop_na(!!dep_var, !!ind_var, !!size_var, !!sites) %>%
        dplyr::select(!!dep_var, !!ind_var, !!size_var, !!sites)


    ###########################################
    #  ANNEALING ALGORITHM
    ###########################################

    model <- get(modelname)

    # Define the parameters to use

    par <- get(paste0(modelname,".par"))
    par_lo <- get(paste0(modelname,".par.lo"))
    par_hi <- get(paste0(modelname,".par.hi"))

    if(is.null(sites)) {
        numsites <- 1
    } else {
        numsites <- length(unique(working_data[,sites]))
    }

    par[rep_vectors] <- map(par[rep_vectors],.f = rep, times = numsites )
    par_lo[rep_vectors] <- map(par_lo[rep_vectors],.f = rep, times = numsites)
    par_hi[rep_vectors] <- map(par_hi[rep_vectors],.f = rep, times = numsites)

    var <- list(X1 = ind_var, X2 = size_var,
                site = sites)


    var$x <- dep_var
    var$mean <- "predicted"
    var$log <- TRUE


    #now call the annealing algorithm, specifying the proper model
    results <- anneal(model, par, var, working_data, par_lo,
                      par_hi, pdf.gauss.prop,
                      dep_var, hessian = TRUE, max_iter = max_iter)


    results <- list(best_pars = results$best_pars,
                    upper_limits = results$upper_limits,
                    lower_limits = results$lower_limits,
                    LL = results$max_likeli,
                    AICc = results$aic_corr,
                    slope = results$slope,
                    R2 = results$R2)


}
#
# first_test <- run_models (moist, dep_var = "ApiG", ind_var = "PPFD",
#             modelname = "site.lin.size.model",size_var = "D10", rep_vectors = 1:2,sites = "Site", max_iter = 2000)



