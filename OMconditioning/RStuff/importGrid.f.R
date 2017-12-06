
importGrid.f <- function(
    gridList=gridY2List,
    gridDir=gridDir,
    covar=F)
{
  for(i in 1:length(gridList)){
    rm(tmp)
    tmp <- SS_output(dir=gridDir %&% gridList[i], repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=covar,forecast=F)
    converged <- file.exists(gridDir %&% gridList[i] %&% "\\" %&% "ss3_opt.std")
    M <- tmp$endgrowth$M
    tmp <- list(tmp$derived_quants, tmp$recruit, tmp$maximum_gradient_component,tmp$likelihoods_used,tmp$likelihoods_raw_by_fleet,tmp$cpue, tmp$Length_comp_Eff_N_tuning_check, tmp$timeseries[,1:8], tmp$parameters[,1:9],M,
         tmp$equil_yield, converged )
    names(tmp) <- c('derived_quants', 'recruit', 'maximum_gradient_component','likelihoods_used',"likelihoods_raw_by_fleet",'cpue', 'Length_comp_Eff_N_tuning_check', 'timeseries', 'parameters','M',
          'equil_yield', 'converged')

    assign(gridList[i],tmp, envir = .GlobalEnv)
  }
  return(999)
}
