#make an assessment grid for BET based on Langley et al (2013)


makeGridB1.f <- function (path='H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\BET\\gridB1', makeGrid=T, splitMasterBatch=T, doHess=T,
  sp.val  = c( "R1P1"), #1 regions, 1 population, no mixing
  h.val   = c( "h70","h80","h90"), #SR steepness
  M.val   = c( "M10", "M08","M06"),#mort
  q.val   = c( "q0", "q1")) # longline catchability
  # t.val =c("t10","t01","t00"))   #tag weight - no tags in aggregate model
  # one population in this grid 
  
{                 
    gridList<- NULL
    nGridElements <- 0

    nGRidElements <- 0
    for (sp in sp.val){
     for (h in h.val){
      for (M in M.val){
       for (iq in q.val){
           gridList <- c(gridList, sp %&% "_" %&% h %&% "_" %&% M %&% "_" %&% iq )
           nGridElements <- nGridElements + 1
        }
       }
      }
     }


      #reminder of options are in the template
        
        
      #define master Batch File  
      masterBatchFName <- 'gridB1'
      if(makeGrid) mf <- file(  path %&% "\\" %&% masterBatchFName %&% ".bat", open="w" )

      if(splitMasterBatch){
        if(makeGrid) mf1 <- file(  path %&% "\\" %&% masterBatchFName %&% "1.bat", open="w" )
        if(makeGrid) mf2 <- file(  path %&% "\\" %&% masterBatchFName %&% "2.bat", open="w" )
        if(makeGrid) mf3 <- file(  path %&% "\\" %&% masterBatchFName %&% "3.bat", open="w" )
        if(makeGrid) mf4 <- file(  path %&% "\\" %&% masterBatchFName %&% "4.bat", open="w" )
      }

      gridElement <- 0
      #loop over all grid factors 
        for( sp in 1:length(sp.val) )
        {
        for( h in 1:length(h.val) )
        {
        for( M in 1:length(M.val) )
        {
        for( iq in 1:length(q.val) )
        {
          gridElement <- gridElement + 1
          # create folder structure and batch files and edit template.DATA.SS and template. CONTROL.SS  

          #define gridName as unique scenario identifier
          gridName <- paste( sp.val[sp],h.val[h], M.val[M], q.val[iq], sep="_" )

          if(makeGrid){ #make the batch files and dirs, otherwise just return the list   
            #create directories and copy requisite files          
            mkdir(c(path %&% '\\' %&% gridName))
            copyString <- " copy "  %&% path %&% "\\gridTemplate " %&% path %&% "\\" %&% gridName 
            system(paste(Sys.getenv("COMSPEC"),"/c ", copyString))


            #create .ctl  file from template files by removing # from switches
            template <- scan(path %&% "\\gridTemplate" %&% "\\templateBET.ctl", what="", sep='\n')
            write ('# ' %&% gridName, file=path %&% "\\" %&% gridName %&% "\\BET.ctl", append=F)

            for(i in 1:length(template)){
              str <- template[i]

              str <- sub('# xxx ' %&%  sp.val[sp], '', str)
              str <- sub('# xxx ' %&%  h.val[h],   '', str)
              str <- sub('# xxx ' %&%  M.val[M],   '', str)
              str <- sub('# xxx ' %&%  q.val[iq], '', str)

              write (str, file=path %&% "\\" %&% gridName %&% "\\BET.ctl", append=T)
            } 

            #create  .dat file from template files by removing # from switches
            template <- scan(path %&% "\\gridTemplate" %&% "\\templateBET.dat", what="", sep='\n')
            write ('# ' %&% gridName, file=path %&% "\\" %&% gridName %&% "\\BET.dat", append=F)

            for(i in 1:length(template)){
              str <- template[i]

              str <- sub('# xxx ' %&%  sp.val[sp], '', str)
              str <- sub('# xxx ' %&%  h.val[h], '', str)
              str <- sub('# xxx ' %&%  M.val[M], '', str)
              str <- sub('# xxx ' %&%  q.val[iq], '', str)

              write (str, file=path %&% "\\" %&% gridName %&% "\\BET.dat", append=T)
            }

            if(doHess) projCall <- "\ncall projBat"
            if(doHess==F) projCall <- "\ncall projBatNoHess"

            
            #write the master Batch file call 
            cat("\ncd " %&% gridName, file=mf, sep='\n' )
            cat(projCall, file=mf, sep='\n' )
            cat("cd ..", file=mf, sep='\n' )

            #make 4 optional Batch files for splitting grid
            if(splitMasterBatch){
              if(gridElement <= 0.25*nGridElements) cat("\ncd " %&% gridName %&% projCall %&% "\ncd .. \n", file=mf1, sep='\n' )
              if(gridElement  > 0.25*nGridElements & gridElement <= 0.50*nGridElements) cat("\ncd " %&% gridName %&% projCall %&% "\ncd .. \n", file=mf2, sep='\n' )
              if(gridElement  > 0.50*nGridElements & gridElement <= 0.75*nGridElements) cat("\ncd " %&% gridName %&% projCall %&% "\ncd .. \n", file=mf3, sep='\n' )
              if(gridElement  > 0.75*nGridElements) cat("\ncd " %&% gridName %&% projCall %&% "\ncd .. \n", file=mf4, sep='\n' )
            }

          } #if makeGrid
        } # end grid factor loops
        }
        }
        }

      if(makeGrid){
        close(mf)
        if(splitMasterBatch){
          close(mf1)
          close(mf2)
          close(mf3)
          close(mf4)
        }
      }         
      #closeAllConnections()
      #graphics.off()

    return(gridList)
}

makeGridB1.f(doHess=F, makeGrid=F)