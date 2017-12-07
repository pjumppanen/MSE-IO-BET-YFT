#  robustness case OM checking the impact of sigmaR

makeGridY17.2noTag.f <- function (path='H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.2noTag', makeGrid=T, splitMasterBatch=10, doHess=T,
#makeGridY17.2noTag.f <- function (path='C:\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.2noTag', makeGrid=T, splitMasterBatch=T, doHess=T,
  sp.val   = c( "R4MvEst"),          #Spatial&Population structure - assessment default single option
  #h.val    = c( "h80"),  #SR steepness
  h.val   = c( "h70","h80","h90"),  #SR steepness
  sigR.val = c( "sr4", "sr6", "sr8"),  #SR steepness
  M.val    = c( "M10", "M08","M06"), #mort
  t.val    = c("t0001"),   #tag weight
  q.val    = c("q0","q1"),           #CPUE q trend % per y
  mix.val  = c("x3"),            #tag mixing period (qtrs); irrelevant oif tags not used
#  ind.val  = c("iH","i10H","iC","i10C"),            #CPUE standardization method: cluster vs HBF  sigma CPUE = 0.3 or 0.1
  ind.val  = c("iH","i10H","iC","i10C"),            #CPUE standardization method: cluster vs HBF  sigma CPUE = 0.3 or 0.1
  sel.val  = c("SS"),           # LL selectivity - stationary of 15 years of devs
#  sel.val  = c("SS","LLd"),    # LL selectivity - stationary of 15 years of devs
  ess.val  = c("CLRW","ess5"))   #size comp weighting; 1 iteration of re-weighting or all fisheries at ess=5

{
    gridList<- NULL

    nGridElements <- 0
    for (sp in sp.val){
     for (h in h.val){
     for (isigR in sigR.val){
      for (M in M.val){
       for (t in t.val){
        for (iq in q.val){
         for (imix in mix.val){
          for (iind in ind.val){
            for (isel in sel.val){
              for (iess in ess.val){
                gridList <- c(gridList, sp %&% "_" %&% h %&% "_" %&% isigR  %&% "_" %&% M %&% "_" %&% t %&% "_" %&% iq %&% "_" %&% imix %&% "_" %&% iind %&% "_" %&% isel %&% "_" %&% iess)
                nGridElements <- nGridElements + 1
              }
            }
          }
         }
        }
       }
      }
      }
     }
    }
      #define master Batch File
      masterBatchFName <- 'gridY17.2noTag'
      if(makeGrid) mf <- file(  path %&% "\\" %&% masterBatchFName %&% ".bat", open="w" )

      if(splitMasterBatch){
        if(makeGrid) mf1 <- file(  path %&% "\\" %&% masterBatchFName %&% "1.bat", open="w" )
        if(makeGrid) mf2 <- file(  path %&% "\\" %&% masterBatchFName %&% "2.bat", open="w" )
        if(makeGrid) mf3 <- file(  path %&% "\\" %&% masterBatchFName %&% "3.bat", open="w" )
        if(makeGrid) mf4 <- file(  path %&% "\\" %&% masterBatchFName %&% "4.bat", open="w" )
        if(makeGrid) mf5 <- file(  path %&% "\\" %&% masterBatchFName %&% "5.bat", open="w" )
        if(makeGrid) mf6 <- file(  path %&% "\\" %&% masterBatchFName %&% "6.bat", open="w" )
        if(makeGrid) mf7 <- file(  path %&% "\\" %&% masterBatchFName %&% "7.bat", open="w" )
        if(makeGrid) mf8 <- file(  path %&% "\\" %&% masterBatchFName %&% "8.bat", open="w" )
      }

ctltemplate <- scan(path %&% "\\gridTemplate" %&% "\\templateYFT.ctl", what="", sep='\n')
dattemplate <- scan(path %&% "\\gridTemplate" %&% "\\templateYFT.dat", what="", sep='\n')

      gridElement <- 0
      #loop over all grid factors
        for( sp in 1:length(sp.val) )
        {
        for( h in 1:length(h.val) )
        {
        for( isigR in 1:length(sigR.val) )
        {
        for( M in 1:length(M.val) )
        {
        for( t in 1:length(t.val) )
        {
        for( iq in 1:length(q.val) )
        {
        for( imix in 1:length(mix.val) )
        {
        for( iind in 1:length(ind.val) )
        {
        for( isel in 1:length(sel.val) )
        {
        for( iess in 1:length(ess.val) )
        {
          gridElement <- gridElement + 1
          # create folder structure and batch files and edit template.DATA.SS and template. CONTROL.SS

          #define gridName as unique scenario identifier
          gridName <- paste( sp.val[sp],h.val[h], sigR.val[isigR] , M.val[M], t.val[t], q.val[iq], mix.val[imix], ind.val[iind], sel.val[isel], ess.val[iess],  sep="_" )

          if(makeGrid){ #make the batch files and dirs, otherwise just return the list
            #create directories and copy requisite files
            #mkdir(c(path %&% '\\' %&% gridName))  #mvbutils versions
            newDir <- path %&% '\\' %&% gridName
            dir.create(newDir)

            copyString <- " copy "  %&% path %&% "\\gridTemplate " %&% path %&% "\\" %&% gridName
            system(paste(Sys.getenv("COMSPEC"),"/c ", copyString))


            #create CONTROL.SS from template by removing comments on appropriate switches
            ctlFile <- file(path %&% "\\" %&% gridName %&% "\\YFT.ctl", open="w")
            cat('# ' %&% gridName, file=ctlFile, sep='\n')
            #write ('# ' %&% gridName, file=path %&% "\\" %&% gridName %&% "\\YFT.ctl", append=F)

            for(i in 1:length(ctltemplate)){
              str <- ctltemplate[i]

              str <- sub('# xxx ' %&%  sp.val[sp],    '', str)
              str <- sub('# xxx ' %&%  h.val[h],      '', str)
              str <- sub('# xxx ' %&%  sigR.val[isigR],      '', str)
              str <- sub('# xxx ' %&%  M.val[M],      '', str)
              str <- sub('# xxx ' %&%  t.val[t],      '', str)
              str <- sub('# xxx ' %&%  q.val[iq],     '', str)
              str <- sub('# xxx ' %&%  mix.val[imix], '', str)
              str <- sub('# xxx ' %&%  ind.val[iind], '', str)
              str <- sub('# xxx ' %&%  sel.val[isel], '', str)

              #write (str, file=path %&% "\\" %&% gridName %&% "\\YFT.ctl", append=T)
              cat(str, file=ctlFile, sep='\n')

            }
            close(ctlFile)

#            template <- scan(path %&% "\\gridTemplate" %&% "\\templateYFT.dat", what="", sep='\n')

            datFile <- file(path %&% "\\" %&% gridName %&% "\\YFT.dat", open="w")
            cat('# ' %&% gridName, file=datFile, sep='\n')
            #write ('# ' %&% gridName, file=path %&% "\\" %&% gridName %&% "\\YFT.dat", append=F)

            for(i in 1:length(dattemplate)){
              str <- dattemplate[i]

            #  str <- sub('# xxx ' %&%  sp.val[sp], '', str)
            #  str <- sub('# xxx ' %&%  h.val[h],   '', str)
            #  str <- sub('# xxx ' %&%  M.val[M],   '', str)
            #  str <- sub('# xxx ' %&%  t.val[t], '', str)
              str <- sub('# xxx ' %&%  q.val[iq], '', str)
              str <- sub('# xxx ' %&%  mix.val[imix], '', str)
              str <- sub('# xxx ' %&%  ind.val[iind], '', str)
              str <- sub('# xxx ' %&%  ess.val[iess], '', str)

              cat(str, file=datFile, sep='\n')
              #write (str, file=path %&% "\\" %&% gridName %&% "\\YFT.dat", append=T)
            }
            close(datFile)

            #create DATA.SS from template
            ### loop over oldtxt and replace switches...
#            template <- scan(path %&% "\\gridTemplate" %&% "\\templateYFT.dat", what="", sep='\n')
#            write ('# ' %&% gridName, file=path %&% "\\" %&% gridName %&% "\\YFT.dat", append=F)
#
#            for(i in 1:length(template)){
#              str <- template[i]
#
#              str <- sub('# xxx ' %&%  sp.val[sp], '', str)
#              str <- sub('# xxx ' %&%  h.val[h], '', str)
#              str <- sub('# xxx ' %&%  M.val[M], '', str)
#              str <- sub('# xxx ' %&%  t.val[t], '', str)
#
#              write (str, file=path %&% "\\" %&% gridName %&% "\\YFT.dat", append=T)
#            }

            if(doHess) projCall <- "\ncall projBat"
            if(doHess==F) projCall <- "\ncall projBatNoHess"


            #write the master Batch file call
            cat("\ncd " %&% gridName, file=mf, sep='\n' )
            cat(projCall, file=mf, sep='\n' )
            cat("cd ..", file=mf, sep='\n' )

            #make nbf optional Batch files for splitting grid
            if(splitMasterBatch){
              if(gridElement <= 0.125*nGridElements) cat("\ncd " %&% gridName %&% projCall %&% "\ncd .. \n", file=mf1, sep='\n' )
              if(gridElement  > 0.125*nGridElements & gridElement <= 0.25*nGridElements) cat("\ncd " %&% gridName %&% projCall %&% "\ncd .. \n", file=mf2, sep='\n' )
              if(gridElement  > 0.25*nGridElements & gridElement <= 0.375*nGridElements) cat("\ncd " %&% gridName %&% projCall %&% "\ncd .. \n", file=mf3, sep='\n' )
              if(gridElement  > 0.375*nGridElements & gridElement <= 0.5*nGridElements) cat("\ncd " %&% gridName %&% projCall %&% "\ncd .. \n", file=mf4, sep='\n' )
              if(gridElement  > 0.50*nGridElements & gridElement <= 0.625*nGridElements) cat("\ncd " %&% gridName %&% projCall %&% "\ncd .. \n", file=mf5, sep='\n' )
              if(gridElement  > 0.625*nGridElements & gridElement <= 0.75*nGridElements) cat("\ncd " %&% gridName %&% projCall %&% "\ncd .. \n", file=mf6, sep='\n' )
              if(gridElement  > 0.75*nGridElements & gridElement <= 0.875*nGridElements) cat("\ncd " %&% gridName %&% projCall %&% "\ncd .. \n", file=mf7, sep='\n' )
              if(gridElement  > 0.875*nGridElements) cat("\ncd " %&% gridName %&% projCall %&% "\ncd .. \n", file=mf8, sep='\n' )
            }

          } #if makeGrid
        } # end grid factor loops
        }
        }
        }
        }
        }
        }
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
          close(mf5)
          close(mf6)
          close(mf7)
          close(mf8)
        }
      }
      #closeAllConnections()
      #graphics.off()

    return(gridList)
}

makeGridY17.2noTag.f(doHess=F, makeGrid=F)