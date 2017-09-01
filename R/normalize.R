#' Fit Parametric Distribution and test
#'
#' This function fits a distribution to a time series of monthly or annual flows for normalization purposes. It alsos run several tests to verify the goodness of fit, using the function fit_perc.
#'
#' @param flow_data Dataframe with columns "year" and "flow". Column "month" should also be included if running monthly.
#' @param time_scale Character of either "monthly" or "annual". 
#' @param distr Character name of the probability distribution function
#' @param ref_period Vector of 2 years to subset the data. Defaults to NULL (full datasets)
#' @param save_prefix Character prefix for saved files. Typically a site name or location.
#' @param save_plots TRUE or FALSE Plots should be saved to disk? Defaults to TRUE
#' @param save_dir Character for save directory. Defaults to working directory
#'
#' @return monthly_ts a dataframe with results of the null model reconstruction
#'
#' @author James Stagge, \email{james.stagge@usu.edu}
#' @seealso \code{\link{fit_perc}}
#' @keywords normalization, percentile
#'
#'
#' @export

fit_norm_dist <- function(flow_series, distr, ref_period=NULL, save_prefix="", save_plots=TRUE, save_dir=getwd()){
 
### Load packages
require(testthat)
require(fitdistrplus)
require(goftest)

### Process object
time_scale <- flow_series$time_scale
flow_ts <- flow_series$ts
save_prefix <- flow_series$site_prefix

### Change column name if water year
if(!is.null(flow_series$wy_first_month)){
	names(flow_ts)[which(names(flow_ts) == "water_year")] <- "year"
}

### Verify that there is a year column
expect_true(all(c("year") %in% names(flow_ts)), label='flow_ts has column named "year"')

### Extract to years in the reference period if reference period is included
if (!is.null(ref_period)){
	year_test <- flow_ts$year >= ref_period[1] & flow_ts$year <= ref_period[2]
	flow_ts <- flow_ts[year_test,]
}

if(time_scale=="monthly") {
################################################
### Fit distribution for normalization of monthly flows and test goodness of fit
#################################################
### Check that flow_ts has the correct column names
expect_true(all(c("year", "month", "flow") %in% names(flow_ts)), label='flow_ts has columns named "year", "month", and "flow"')

### If only one distribution given, repeat it 12 times
if(length(distr) == 1){
distr <- rep(distr, 12)
}

### Loop through 12 months
	for (j in 1:12) {
		### Create a test for the month and extract flows from observed for this month
		month_test <- flow_ts$month == j
		month_flows <- flow_ts$flow[month_test]
	
		### Run the distribution fit
		month_name <- paste0("month_",j)
		month_fit_temp <- fit_dist_generic(flows=month_flows, distr=distr[j], save_prefix=paste0(save_prefix,"_",month_name), save_plots=save_plots, save_dir=save_dir)
		
		if (j ==1) {
			fit_param <- month_fit_temp
			
			### Trick to prepare fit as a list
			fit_1 <- fit_param$fit
			fit_param$fit <- list()
			fit_param$fit[[1]] <- fit_1
		} else {
			fit_param$fit[[j]] <- month_fit_temp$fit
			fit_param$gof <- rbind(fit_param$gof, month_fit_temp$gof)
		}
	}	
} else if (time_scale=="annual") {
################################################
### Fit distribution for normalization of annual flows and test goodness of fit
#################################################
### Check that flow_ts has the correct column names
expect_true(all(c("year", "flow") %in% names(flow_ts)), label='flow_ts has columns named "year" and "flow"')
### Check that there is only 1 measurement per year
expect_true(all(table(flow_ts$year)==1), label='flow_ts has 1 measurement per year')

### Run perc_fit on annual flows
fit_param <- fit_dist_generic(flows=flow_ts$flow, distr=distr, save_prefix=paste0(save_prefix,"_annual"), save_plots=save_plots, save_dir=save_dir)

}

### Convert to paleo.norm object and return
fit_param <- paleo.norm(fit=fit_param$fit, distr=distr, gof=fit_param$gof, prefix=fit_param$prefix, ref_period=ref_period, time_scale=time_scale)

return(fit_param)
}


#' Distribution Fit (generic function)
#'
#' This function is called by fit_norm_dist.  It fits a distribution to a time series and runs several tests to verify the goodness of fit.
#'
#' @param flows Vector of numeric flows to be fit.
#' @param distr Character name of the probability distribution function
#' @param save_plots TRUE or FALSE Plots should be saved to disk? Defaults to TRUE
#' @param save_dir Character for save directory. Defaults to working directory
#'
#' @return monthly_ts a dataframe with results of the null model reconstruction
#'
#' @author James Stagge, \email{james.stagge@usu.edu}
#' @seealso \code{\link{flow_perc_fit}}
#' @keywords normalization, percentile
#'
#'
#' @export
fit_dist_generic <- function(flows, distr, save_prefix, save_plots=TRUE, save_dir) {
require(staggefuncs)
require(fitdistrplus)

	if(save_plots==TRUE){
		### Create folder for output
		dir.create(file.path(save_dir, "pdf"), showWarnings = FALSE)
		dir.create(file.path(save_dir, "png"), showWarnings = FALSE)
	}
	
	### extract the non NA values
	flows <- flows[!is.na(flows)]
	
	####################################################
	### Plot and save a Cullen-Fry diagram for flows
	####################################################
	### Create Cullen-Fry with bootstrap
	descdist(flows, boot=500)
	
	if(save_plots==TRUE){
		#### Save to both pdf and png
		pdf_location <- file.path(save_dir, paste0("pdf/",save_prefix,"_",distr,"_cullenfrey.pdf"))
		d_pdf = dev.copy(pdf,pdf_location,width=6, height=6)
    	dev.off(d_pdf)
	
		png_location <- file.path(save_dir, paste0("png/",save_prefix,"_",distr,"_cullenfrey.png"))
		d_png = dev.copy(png,png_location,width=6, height=6,units="in",res=600)
    	dev.off(d_png)
	}
	####################################################
	### Fit the distribution to flows
	####################################################
	### This function fits a distribution to flows
	### Could come back and do a bootstrap to estimate the uncertainty and carry
	### through the calculations.
	flow_fit <- try(fitdist(flows, distr))

	### Can play with this if the future
	#a_kernel_fit <- spdfit(month_flows, type="pwm")
	#pspd(c(1,2,3,4,5), a_kernel_fit)

	####################################################
	### Calculate Goodness of fit and return fit results
	####################################################
	### If the fit runs, collect results and goodness of fit statistics
	if (class(flow_fit)!= "try-error") {
		### Create a dataframe with results
		fit_coef <- data.frame(distr=distr,t(as.matrix(flow_fit$estimate)))
		rownames(fit_coef) <- save_prefix
		
		### Add goodness of fit (AIC) and bootstrapped p-values from KS, AD, and CVM tests
		fit_gof <- data.frame(aic = gofstat(flow_fit)$aic)
		rownames(fit_gof) <- save_prefix
		p_tests <- gof_bootstrap(flow_fit, n_sims=5e3, parallel=TRUE)
		fit_gof$ks <- p_tests$ks_p
		fit_gof$ad <- p_tests$ad_p
		fit_gof$cvm <- p_tests$cvm_p
		
		fit_result <- list(fit=flow_fit, coef=fit_coef, gof=fit_gof, prefix=save_prefix)
	}
	return(fit_result)	
}






#' Convert flow series to normalized percentile
#'
#' @export

flow_to_norm <- function(flow_series, dist_object){

### Extract flow time series
flow_ts <- flow_series$ts
### Create object to hold result
perc_est <- rep(NA, dim(flow_ts)[1])

### Run for monthly
if(flow_series$time_scale == "monthly"){
	for (j in seq(1,12)) {
		### Extract only this month
		month_test <- flow_ts$month == j
		### Convert distribution parameters to a list
		param_list <- list(q = flow_ts$flow[month_test])
		### Save to parameter list
		param_name <- names(dist_object$fit[[j]]$estimate)
		param_list[param_name] <- dist_object$fit[[j]]$estimate
		### Create a p function for the given annual distribution and run to obtain
		### Annual time series of percentiles
		p_distr <- match.fun(paste("p",dist_object$fit[[j]]$distname,sep=""))
		perc_est[month_test] <- do.call(p_distr, param_list)
	}
} else {
	### Convert distribution parameters to a list
	param_list <- list(q = flow_ts$flow)
	### Save to parameter list
	param_name <- names(dist_object$fit$estimate)
	param_list[param_name] <- dist_object$fit$estimate
	### Create a p function for the given annual distribution and run to obtain
	### Annual time series of percentiles
	p_distr <- match.fun(paste("p",dist_object$fit$distname,sep=""))
	perc_est <- do.call(p_distr, param_list)

}

### Convert to normalized scale and return
norm_est <- qnorm(perc_est)
return(norm_est)
}





#' Convert normalized percentile to flow series
#'
#' @export

norm_to_flow <- function(norm_series, dist_object){

### Extract flow time series
norm_ts <- norm_series$ts
norm_ts$perc <- pnorm(norm_ts$norm)
### Create object to hold result
flow_est <- rep(NA, dim(norm_ts)[1])

### Run for monthly
if(norm_series$time_scale == "monthly"){
	for (j in seq(1,12)) {
		### Extract only this month
		month_test <- norm_ts$month == j
		### Convert distribution parameters to a list
		param_list <- list(p = norm_ts$perc[month_test])
		### Save to parameter list
		param_name <- names(dist_object$fit[[j]]$estimate)
		param_list[param_name] <- dist_object$fit[[j]]$estimate
		### Create a p function for the given annual distribution and run to obtain
		### Annual time series of percentiles
		q_distr <- match.fun(paste("q",dist_object$fit[[j]]$distname,sep=""))
		flow_est[month_test] <- do.call(q_distr, param_list)
	}
} else {
	### Convert distribution parameters to a list
	param_list <- list(q = norm_ts$perc)
	### Save to parameter list
	param_name <- names(dist_object$fit$estimate)
	param_list[param_name] <- dist_object$fit$estimate
	### Create a p function for the given annual distribution and run to obtain
	### Annual time series of percentiles
	q_distr <- match.fun(paste("q",dist_object$fit$distname,sep=""))
	flow_est <- do.call(q_distr, param_list)

}

return(flow_est)

}




