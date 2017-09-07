#' Paleo Time Series Constructor
#'
#' This function provides a standard time series plot using ggplot's geom_line.
#'
#' @param data A dataframe
#' @param write_folder Folder in which to save results
#' @param write_file Name of file to be saved
#'
#' @return p Plot
#'
#'
#' @export

paleo.ts <- function(ts, time_scale, site_prefix="", wy_first_month=NULL, ...){
	### Put in a bunch of if statements making sure it is in the right format
	### use stop command to return error
	
	x <- list(ts=ts, time_scale=time_scale, site_prefix=site_prefix, wy_first_month=wy_first_month, ...)

	class(x) <- "paleo.ts"
	return(x)
}

#' Paleo Time Series Print Function
#'
#' This function provides a standard time series plot using ggplot's geom_line.
#'
#' @param data A dataframe
#' @param write_folder Folder in which to save results
#' @param write_file Name of file to be saved
#'
#' @return p Plot
#'
#'
#' @export

print.paleo.ts <- function(x, ...) {
  n_nonna <- sum(!is.na(x$ts$flow)) 
  
  cat("Object of class", class(x), "with", n_nonna, "records at", x$time_scale, "scale\nSite prefix is", x$site_prefix, "\n\n")
  print(head(x$ts))
  cat("...\n")
  print(tail(x$ts))
# print(unclass(head(x)))
}




#' Generic No documentation
#'
#' This function needs documentation.
#'
#' @param data A dataframe
#' @param write_folder Folder in which to save results
#' @param write_file Name of file to be saved
#'
#' @return p Plot
#'
#'
#' @export
### Object to hold normalization
paleo.norm <- function(fit, distr, gof, prefix, ref_period, time_scale){
	### Put in a bunch of if statements making sure it is in the right format
	### use stop command to return error
	x <- list(fit=fit, distr=distr, gof=gof, prefix=prefix, ref_period=ref_period, time_scale=time_scale)

	class(x) <- "paleo.norm"
	return(x)
}


#' Generic No documentation
#'
#' This function needs documentation.
#'
#' @param data A dataframe
#' @param write_folder Folder in which to save results
#' @param write_file Name of file to be saved
#'
#' @return p Plot
#'
#'
#' @export
print.paleo.norm <- function(x, ...) {
  cat("Object of class", class(x), "at", x$time_scale, "scale\n\n")
  print(x$coef)
  cat("\n")
  print(x$gof)
# print(unclass(head(x)))
}

#' Generic No documentation
#'
#' This function needs documentation.
#'
#' @param data A dataframe
#' @param write_folder Folder in which to save results
#' @param write_file Name of file to be saved
#'
#' @return p Plot
#'
#'
#' @export
plot.paleo.norm <- function(x, ...) {
  if (x$time_scale == "annual") {
  	plot(x$fit)
  } else if (x$time_scale == "monthly") {
  	for (j in seq(1,length(x$fit))){
  		plot(x$fit[[j]])
  		invisible(readline(prompt= paste0("Month ",j, "\nPress [enter] to continue")))
  	}
  
  }
}


#' Generic No documentation
#'
#' This function needs documentation.
#'
#' @param data A dataframe
#' @param write_folder Folder in which to save results
#' @param write_file Name of file to be saved
#'
#' @return p Plot
#'
#'
#' @export
paleo.fit <- function(method, reconst_data, annual_norm=NULL, monthly_norm=NULL, first_month_wy=10, reg_eq=NULL, pred_ts=NULL){
	x <- list(method=method, reconst_data=reconst_data, first_month_wy=first_month_wy)
	
	if(method == "mf"){
		model_list <- list()	
	} else if (method == "ap") {
		model_list <- list(annual_norm=annual_norm, monthly_norm=monthly_norm)
	} else if (method == "apr") {
		model_list <- list(reg_eq=reg_eq, annual_norm=annual_norm, monthly_norm=monthly_norm)
	}
	
	x <- c(x, model_list)
	
	class(x) <- "paleo.fit"
	return(x)
}


#' Generic No documentation
#'
#' This function needs documentation.
#'
#' @param data A dataframe
#' @param write_folder Folder in which to save results
#' @param write_file Name of file to be saved
#'
#' @return p Plot
#'
#'
#' @export
coef.paleo.fit <- function(fit, ...) {
	if (fit$method == "ap") {
		coef_mat <-  matrix(rep(1,12),1, 12)
		rownames(coef_mat) <- "annual_norm" 
	} else if (fit$method == "apr") {
		coef_mat <- sapply(fit$reg_model, function(x) coef(x))
		coef_mat <- as.matrix(coef_mat)
	}
	
	colnames(coef_mat) <- paste0("M_", seq(1,12))
	return(coef_mat)
}




