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
### Generic wrapper for fitting model
fit_model <- function(method, regmethod="lm", reconst_data, annual_norm=NULL, monthly_norm=NULL, pred_ts=NULL, reg_eq=NULL, monthly_obs=NULL,...){

	################################################
	### Process Reconstructed time series
	#################################################	
		flow_ts <- reconst_data$ts  
		### Test if the reconstruction uses calendar years
		are_cal_years <- is.null(reconst_data$wy_first_month)

		### Create a monthly time series of all combinations of months and year
		### For calendar years, need the range of years and 12 months
		if (are_cal_years){
			year_range <- c(min(flow_ts["year"], na.rm=TRUE), max(flow_ts["year"], na.rm=TRUE))
			monthly_ts <- expand.grid(month = seq(1,12), year = seq(year_range[1], year_range[2]))
		} else {
		### For water years, go one year before and after to make sure it works
		### Then calculate water years
			year_range <- c(min(flow_ts["water_year"], na.rm=TRUE)-1, max(flow_ts["water_year"], na.rm=TRUE)+1)
			monthly_ts <- expand.grid(month = seq(1,12), year = seq(year_range[1], year_range[2]))
			monthly_ts$water_year <- usgs_wateryear(monthly_ts$year, monthly_ts$month, reconst_data$wy_first_month)
		}

		### Set an order column to allow easier resorting at the end
		monthly_ts$t <- seq(1,dim(monthly_ts)[1])		

	################################################
	### Begin to fit model
	#################################################			
	if(method == "mf"){

		### Merge annual flows and month
		monthly_ts <- merge(monthly_ts, flow_ts, all.x=TRUE)
		### Rename flow column to annual recon
		names(monthly_ts)[which(names(monthly_ts) == "flow")] <- "annual_recon"
				
		### Re-sort to obtain time series
		monthly_ts <- monthly_ts[with(monthly_ts, order(t)), ]

		### Process observed flow
		obs_ts <- monthly_obs$ts  

		### Test if the reconstruction uses calendar years
		are_cal_years <- is.null(reconst_data$wy_first_month)

		### Run the MF model
		if(are_cal_years == TRUE) {
			mf_prop <- mf_fit(obs_ts=obs_ts, wy_first_month=1)
		} else {
			mf_prop <- mf_fit(obs_ts=obs_ts, wy_first_month=reconst_data$wy_first_month)
		}
			
		### Prepare to return fit object
		x <- paleo.fit(method=method, reconst_data=reconst_data, mf_prop=mf_prop)			
		x$reconst_data$ts <- monthly_ts
		x$reconst_data$time_scale <- "monthly"
	} 
	else if (method == "ap" | method == "apr") {
		### Check that annual reconstructed timeseries and normalization are from same data
		expect_that(annual_norm$prefix, is_identical_to(paste0(reconst_data$site_prefix, "_annual")))
		
 		### Calculate normalized annual flow
 		flow_ts$annual_norm <- flow_to_norm(flow_series=reconst_data, dist_object=annual_norm)	
		
		### Merge annual flows and month
		monthly_ts <- merge(monthly_ts, flow_ts, all.x=TRUE)
		### Rename flow column to annual recon
		names(monthly_ts)[which(names(monthly_ts) == "flow")] <- "annual_recon"
				
		### Re-sort to obtain time series
		monthly_ts <- monthly_ts[with(monthly_ts, order(t)), ]

	################################################
	### Separate methods for AP and APR
	#################################################	
		### For AP model
		if (method == "ap") {
			### Reformat column names	
			rownames(monthly_ts) <- monthly_ts$t
			monthly_ts <- subset(monthly_ts, select=-c(t))		
		
			### No regression model needed
			reg_model <- rep(list(NA),12)
		
		### For APR model
		} else if (method == "apr") {
			
			### Merge observed flows and month
			monthly_obs_ts_df <- monthly_obs$ts
			names(monthly_obs_ts_df)[which(names(monthly_obs_ts_df) == "flow")] <- "obs_flow"
 			monthly_obs_ts_df$obs_norm <- flow_to_norm(flow_series=monthly_obs, dist_object=monthly_norm)	

			### Merge annual flows and month
			monthly_ts <- merge(monthly_ts, monthly_obs_ts_df, all.x=TRUE)
		
			### Merge predictors
			for (j in seq(1,length(pred_ts))) {
				monthly_ts <- merge(monthly_ts, pred_ts[[j]]$ts, all.x=TRUE)
			}
			
			### Re-sort to obtain time series
			monthly_ts <- monthly_ts[with(monthly_ts, order(t)), ]
			
			### Run the APR model
			reg_model <- apr_fit(recon_data=monthly_ts, reg_eq=reg_eq, regmethod=regmethod)
			monthly_ts <- reg_model$pred_df
			reg_model <- reg_model$reg_model
		
		}
			
		### Prepare to return fit object
		x <- paleo.fit(method=method, reconst_data=reconst_data, annual_norm=annual_norm, monthly_norm=monthly_norm)
		x$reconst_data$ts <- monthly_ts
		x$reconst_data$time_scale <- "monthly"
		x$reg_model <- reg_model	
	}

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
flow_reconstr <- function(recon_model, post_proc=FALSE, interpolate=FALSE){
	
	### Process reconstruction object
	method <- recon_model$method
	site_prefix <- recon_model$reconst_data$site_prefix
	
	if(method == "mf"){
		recon_ts <- recon_model$reconst_data
		
		### Merge with MF proportion
		mf_ts <- merge(recon_ts$ts, recon_model$mf_prop, by="month", all.x=TRUE)
		### Multiply monthly fraction by annual reconstruction * 12
		mf_ts$flow_est <- mf_ts$prop_mean * mf_ts$annual_recon*12
		### Re-sort to obtain time series
		mf_ts <- mf_ts[with(mf_ts, order(t)), ]
		
		recon_ts$ts <- mf_ts
		
		### Create object to return
		monthly_rec_ts <- paleo.ts(ts=recon_ts, time_scale="monthly", site_prefix=site_prefix, wy_first_month=recon_model$first_month_wy, method=method)
			
	} else if (method == "ap" | method == "apr") {	
	
		if (method == "ap") {
			recon_ts <- recon_model$reconst_data
			recon_ts$ts$norm_est <- recon_ts$ts$annual_norm
		} else if (method == "apr") {
			recon_ts <- apr_reconst(recon_model=recon_model)
		}
		
		### Post-Processing option
		if (post_proc == TRUE){
			### Determine mean and std dev for predicted during reference period
			ref_period <- recon_model$monthly_norm$ref_period		
			ref_test <- recon_ts$ts$year >= ref_period[[1]] & recon_ts$ts$year <= ref_period[[2]]
						
			mean_ref <- rep(NA, 12)
			sd_ref <- rep(NA, 12)
									
			for (j in seq(1,12)){
				### Set reference period and monthly test
				month_test <- recon_ts$ts$month == j
				ref_month_test <- ref_test & month_test
				### Extract the mean and standard deviation for norm estimate of each month wihin the refernce period
				mean_ref[j] <- mean(recon_ts$ts$norm_est[ref_month_test], na.rm=TRUE)
				sd_ref[j] <- sd(recon_ts$ts$norm_est[ref_month_test], na.rm=TRUE)
				### Apply this to the entire time series within the month
				recon_ts$ts$norm_est[month_test] <- qnorm(pnorm(recon_ts$ts$norm_est[month_test], mean_ref[j], sd_ref[j]))
				recon_ts$post_proc <- list(mean=mean_ref, sd=sd_ref)
			}	
		}
		
		### Change column name to norm
		column_test <- names(recon_ts$ts) == "norm_est"
		names(recon_ts$ts)[column_test] <- "norm"
		
		### Apply distribution to convert back to flows
		flow_est <- norm_to_flow(norm_series=recon_ts, dist_object=recon_model$monthly_norm)
		recon_ts$ts$flow_est <- flow_est
		names(recon_ts$ts)[column_test] <- "norm_est"
		
		### Create object to return
		monthly_rec_ts <- paleo.ts(ts=recon_ts, time_scale="monthly", site_prefix=site_prefix, wy_first_month=recon_model$first_month_wy, method=method)
	}
	
return(monthly_rec_ts)
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
mf_fit <- function(obs_ts, wy_first_month=1) {
require(data.table)

### Add a water year column
obs_ts$water_year <- usgs_wateryear(obs_ts$year, obs_ts$month, wy_first_month)

################################################
### Calculate Annual Flow
#################################################
obs_ts <- data.table(obs_ts)

### Calculate flow by water year 
flow_obs_annual <- obs_ts[,list(annual_sum=sum(flow), annual_mean=mean(flow)),by="water_year"]
### Merge back with observed flow
obs_ts <-  merge(obs_ts,flow_obs_annual, by=c("water_year"), all.x = T)
	

################################################
### Calculate Monthly Fraction
#################################################
mf_ts <- obs_ts

### Only consider those water years with 12 months of data, so first sum months with flows greater than 0.  Negative flows will confuse the proportion method.
months_by_wy <- mf_ts[,list(count_months=sum(flow>0, na.rm=TRUE)), by=water_year]
	
### Determine full water years and subset null_dt to only full years 
full_wys <- months_by_wy$water_year[months_by_wy$count_months==12]
mf_ts <- mf_ts[mf_ts$water_year %in% full_wys,]

### Calculate monthly proportion 
mf_ts$prop_month <- mf_ts$flow/mf_ts$annual_sum

### Summarize null_dt by calculating the mean and median proportion across all years
null_prop <- mf_ts[,list(prop_mean=mean(prop_month, na.rm=TRUE)), by=month] 

### Rescale so that months always add up to 1
null_prop$prop_mean <- null_prop$prop_mean/sum(null_prop$prop_mean)
null_prop <- data.frame(null_prop)

### Re-sort by months
null_prop <- null_prop[order(null_prop$month),]

### Return monthly proportion
return(null_prop)

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
apr_fit <- function(recon_data, reg_eq, regmethod) {

	################################################
	### Process regression equation
	################################################# 
	reg_vars <- gsub("[[:space:]]", "", reg_eq)
	reg_vars <- c(sapply(strsplit(reg_vars, split='+', fixed=TRUE), `[`))

	### Determine which variables are lagged
	lag_test <- grepl("l(", reg_vars, fixed=TRUE)
	
	vars_plain <- reg_vars[!lag_test]
	vars_lagged <- reg_vars[lag_test]
	
	################################################
	### Process lagged predictors to create predictor data
	################################################# 
	### Remove lagging parentheses
	vars_lagged <- substr(vars_lagged, 3, nchar(vars_lagged)-1)
	### Split at first comma
	XX <- "SPLITTERTEXT"
	vars_lagged <- strsplit(sub(",\\s*", XX, vars_lagged), XX)
	vars_lagged_list <- c()
	
	### Loop through all combinations to create lagged predictor data
	### Requires shift function from staggefuncs
	for (j in seq(1,length(vars_lagged))) {
		### Calculate lags for variable	
		name_j <- vars_lagged[[j]][1]
		if(name_j == "flow") {name_j <- "annual_norm"}
		lags_j <- eval(parse(text=vars_lagged[[j]][2]))
			
		### create dataframe
		base_data <- recon_data[name_j]
		
		### Loop through all lags and create the lagged predictor dataframe
		for (k in seq(1, length(lags_j))){
			lag_k <- lags_j[k]
			lagged_temp <- shift(as.matrix(base_data), shift_by=lag_k*12)
			lagged_temp <- data.frame(lagged_temp)
			names(lagged_temp) <- paste0(name_j,".",lag_k)
			if(k == 1) {
				lagged_df <- lagged_temp
			} else {
				lagged_df <- cbind(lagged_df, lagged_temp)
			}
		}
		### Combine all variables
		if (j == 1){
			lagged_pred <- lagged_df
		} else {
			lagged_pred <- cbind(lagged_pred, lagged_df)
		}
		
	}

	################################################
	### Extract plain variables and combine with lagged to create final predictor set
	################################################# 
	### change flow to annual_perc
	vars_plain[vars_plain %in% "flow"] <- "annual_norm"
	### Extract these variables
	plain_pred <- recon_data[vars_plain]
	### Create predictor dataframe
	pred_df <- cbind(recon_data[c("year", "month")],plain_pred, lagged_pred)


	################################################
	### Extract monthly percentiles
	#################################################  
	monthly_obs <- recon_data[c("year", "month", "obs_norm")]


   	################################################
	### Fit monthly model
	################################################# 
	### If elastic net regression
	if (regmethod == "elastic"){	  

	require(caret)
	require(doParallel)
	### Set elastic grid search parameters
	lambda_grid <- 10^seq(2,-6,length=100)
	#lambda_grid <- seq(0,1,length=10)
	alpha_grid <- seq(0,1,length=10)
	eGrid <- expand.grid(.alpha = alpha_grid, .lambda = lambda_grid)

	### Set control for runs, 10-fold repeated measures, 4 iterations
	Control <- trainControl(method = "repeatedcv",repeats = 8, number=10, verboseIter =FALSE)
	
	### Create parallel clusters
	cores <- detectCores()
	cl <- makePSOCKcluster(cores)
	registerDoParallel(cl)
	
	for (j in seq(1,12)) {
		### Extract to months
		pred_j <- pred_df[pred_df$month == j, ]
		norm_j <- monthly_obs[monthly_obs$month == j,][c("year", "month", "obs_norm")]
		
		### Merge to ensure complete cases
		pred_and_norm <- merge(norm_j, pred_j)
		pred_and_norm <- pred_and_norm[complete.cases(pred_and_norm), !(names(pred_and_norm) %in% c("year", "month"))]
			
		### Extract the values
		x <- as.matrix(pred_and_norm[,!(names(pred_and_norm) %in% c("year", "month", "obs_norm"))])
		y <- pred_and_norm$obs_norm
					
		### Fit the model
		netFit <- train(x =x, y = y,
          method = "glmnet",
          tuneGrid = eGrid,
          trControl = Control)

		### Fit the model
		netFit <- train(obs_norm ~ . , data = pred_and_norm,
          method = "glmnet",
          tuneGrid = eGrid,
          trControl = Control)
    	
		### Extract the best tuning parameters
		my_glmnet_model <- netFit$finalModel
		
		### Extract coefficients
		model_coef <- coef(my_glmnet_model, s = netFit$bestTune$lambda)		
		model_results <- matrix(model_coef)
		rownames(model_results) <- dimnames(model_coef)[[1]]
	
		### Combine coefficients to a single matrix
		if (j == 1) {
			final_coef <- model_results
		} else {
			final_coef <- cbind(final_coef, model_results)
		}
	}
	
	### Stop clusters
	stopCluster(cl)
	
	} else if (regmethod == "lm") {
	
		for (j in seq(1,12)) {
		### Extract to months
		pred_j <- pred_df[pred_df$month == j, ]
		norm_j <- monthly_obs[monthly_obs$month == j,][c("year", "month", "obs_norm")]
		
		### Merge to ensure complete cases
		pred_and_norm <- merge(norm_j, pred_j)
		pred_and_norm <- pred_and_norm[complete.cases(pred_and_norm), !(names(pred_and_norm) %in% c("year", "month"))]
			
		### Extract the values
		x <- as.matrix(pred_and_norm[,!(names(pred_and_norm) %in% c("year", "month", "obs_norm"))])
		y <- pred_and_norm$monthly_norm

		### Fit the model
		netFit <- lm(obs_norm ~ . , data = pred_and_norm)

		### Extract the best tuning parameters
		my_glmnet_model <- netFit
		
		### Save regression model
		if (j ==1) {reg_model_month <- list()}
		reg_model_month[[j]] <- my_glmnet_model
				
		}
	}
	
	### return results
	return(list(reg_model=reg_model_month, pred_df=pred_df))
  
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
apr_reconst <- function(recon_model) {

	### Process object
	monthly_ts <- recon_model$reconst_data$ts  
	monthly_ts$t <- rownames(monthly_ts)
	
	### Create column to hold prediction
	monthly_ts$norm_est <- NA
	
	### Loop through months and predict value 
	for (j in seq(1,12)) {
		month_test <- monthly_ts$month == j
		monthly_ts$norm_est[month_test] <- predict(recon_model$reg_model[[j]], monthly_ts[month_test,])
	}
	
	################################################
	### Reformat results to return
	#################################################
	### Re-sort
	monthly_ts <- monthly_ts[with(monthly_ts, order(year, month)), ]
	monthly_ts$t <- seq(1, dim(monthly_ts)[1])
	rownames(monthly_ts) <- monthly_ts$t
	### Cut to columns
	monthly_ts <- monthly_ts[c("month", "year", "t", "norm_est")]
	
	### Create reconstructed object
	recon_result_ts <- recon_model$reconst_data
	recon_result_ts$ts <- monthly_ts	
	
  return(recon_result_ts)
}

