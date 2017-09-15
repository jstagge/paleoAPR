#' Shift Indices
#'
#' No description
#'
#' @param annual_rec a dataframe with columns "water_year" and "annual_flow"
#' @param monthly_prop a dataframe with columns "month" and "prop"
#' @param first_month_wy a numeric variable with the month that signifies the start of the water year, usually 10
#'
#' @return monthly_ts a dataframe with results of the null model reconstruction
#'
#' @examples
#' mf_model()
#'
#' @export


shift_indices <- function(data, name, lags) {
	for (j in seq(1,length(lags))) {
		shift_by <- lags[j]
		data_temp <- shift_bylag(data, shift_by=shift_by)
	
		if (j == 1){
			data_final <- data_temp
		} else {
			data_final <- cbind(data_final, data_temp)
		}
	}

	colnames(data_final) <- paste0(name, "_", lags)
	return(data_final)
}


#' Shift Indices Dataframe
#'
#' No description
#'
#' @param annual_rec a dataframe with columns "water_year" and "annual_flow"
#' @param monthly_prop a dataframe with columns "month" and "prop"
#' @param first_month_wy a numeric variable with the month that signifies the start of the water year, usually 10
#'
#' @return monthly_ts a dataframe with results of the null model reconstruction
#'
#' @examples
#' mf_model()
#'
#' @export

shift_indices_df <- function(data, lags, date_firstcol=TRUE) {

### Extract dates if the first column contains these
if(date_firstcol==TRUE){
dates <- data[,1]
dates_name <- colnames(data)[[1]]
data <- data[,seq(2,dim(data)[2])]
}

### Loop and apply lag to all other columns
for (j in seq(1,dim(data)[2])) {
	data_temp <- shift_indices(data[,j], colnames(data)[j], lags)
	if (j == 1) {
		data_lag <- data_temp
	} else {
		data_lag <- cbind(data_lag, data_temp)
	}
}

### Recombine if first column contained dates
if(date_firstcol==TRUE){
data_lag <- data.frame(dates, data_lag)
### Rename first column
colnames(data_lag)[[1]] <- dates_name

}

return(data_lag)
}


