#' @title RDPI

#' @description Function to compute the RDPI (Relative Distance Plasticity Index,
#' Valladares et al, (2006) Quantitative estimation of phenotypic plasticity:
#'  bridging the gap between the evolutionary concept and its ecological applications,
#'  Journal of Ecology, 94(6):1103-1116.
#'

#' @param data The dataframe that contains the data
#' @param sp The bare (unquoted) name of the column whose values will be used as independent variable.
#'The function will compare RDPI values among values of this variable. It can be species, provenances, etc.
#' @param trait The bare (unquoted) name of the column that holds the trait for which to calculate RDPI. Must be numeric
#' @param factor the bare (unquoted) name of the column that holds the environmental factor for which we will calculate RDPI. By definition, RDPI computes distances between pairs of
#' observations that are at different levels of this factor.
#'
#'
#' @return This function computes RDPI to the environmental factor for each species (or any other variable defined in 'sp')
#' of the dataset
#' Then it makes an ANOVA or t-test of the values of RDPI across species
#' and plots the boxplot
#'
#' @examples
#' data(ecophysio)
#' rdpi(ecophysio,sp,SB, Piso)
#'
#' @export


rdpi <- function(data, sp, trait, factor) {

  # Load the required libraries

  library(agricolae)
  library(psych)
  library(dplyr)
  library(sciplot)

  # Modify the parameters to add quotes (")
  sp <- deparse(substitute(sp))
  trait <- deparse(substitute(trait))
  factor <- deparse(substitute(factor))

  # Create the object (an empty data frame) that will store the results of RDPI
  RDPI <- data.frame(species=character(0), value=numeric(0))

  # Since we need to compute everything per species, we make a 'for' loop

  for (a in 1:nlevels(data[[sp]])) {

    # subset the data for a given species
    data<-data[data[[sp]]== levels(data[[sp]])[a],]


    # NOTE: RDPI is based on pairwise distances (Canberra distance) between individuals that
    # belong to different levels of an environmental variable. We perform this in three steps:

    # Step1: Compute pairwise Canberra distance (aka RDPI) for all individuals in the dataset
    RDPI_temp<-as.matrix(dist(x=data[[trait]], method="canberra"))

    # Step 2: Generate a matrix where value is "TRUE" only if observation i and observation j
    # belong to different levels of the factor
    filter_frame <- data.frame(matrix(NA,nrow(data),nrow(data)))
    for (i in 1:nrow(filter_frame)) {
      for (j in 1:ncol(filter_frame)){
        ifelse((data[[factor]][i]==data[[factor]][j]),
               filter_frame[i,j] <- FALSE,
               filter_frame[i,j] <- TRUE)
      }
    }
    filter_frame[upper.tri(filter_frame,diag = T)] <- FALSE         #only keep lower triangle


    # Step 3: Subset RDPI so that it only includes comparisons between individuals that
    # belong to different levels of an environmental variable
    RDPI_temp <- RDPI_temp[filter_frame==TRUE]
    if(length(RDPI_temp)>0) { ########### beginning of addition
      RDPI_sp <- data.frame(species=levels(data[[sp]])[a],
                          value=RDPI_temp)

      RDPI<- rbind(RDPI, RDPI_sp)
    } ########### end of addition
  }
  summary <- RDPI %>%
    group_by(species) %>%
    summarise(mean=mean(value,na.rm=T),
              sd=sd(value,na.rm=T),
              se= se(value, na.rm=T))
  print(summary)

  boxplot(RDPI$value~RDPI$species)


  if (nlevels(RDPI$species) < 3) {
    fit <- t.test(RDPI$value ~ RDPI$species)
    print("t-test")
    print(fit)
  } else {
    fit<-aov(RDPI$value ~ RDPI$species)

    print("ANOVA test")
    print(summary(fit))
    Tuk<-HSD.test (fit, trt='RDPI$sp')
    print("Tukey HSD test for differences accross groups")
    print(Tuk$groups)
  }

  return(RDPI)

}
