#### MMoL project: pooled AKM or sector AKM with switch
rm(list = ls())

library(haven)
library(data.table)
library(fixest)
############################################################
# TEMPORARY CHECK: inspect raw variable names in the .dta file
# Remove after checking
############################################################

vNamesRaw <- names(read_dta("Veneto_panel_1984_2001.dta", n_max = 0))
print(vNamesRaw)
############################################################
# 0. User settings
#
# sAKMMode:
# - "pooled_then_sector"
# - "sector_by_sector"
#
# Interpretation:
# pooled_then_sector:
#   estimate one AKM on the pooled LCC sample, then describe
#   decomposition inside each sector using the pooled FE.
#
# sector_by_sector:
#   estimate one AKM separately inside each sector, using the
#   pooled LCC sample as the starting sample.
#
# This lets you compare the two empirical strategies directly.
############################################################

sDataPath <- "Veneto_panel_1984_2001.dta"
sAKMMode <- "pooled_then_sector" # Here we can change if we use AKM on whole panel or by sector
iMinObsSector <- 10000

############################################################
# 1. Helper: build cleaned AKM sample
#
# We repeat the exact same cleaning logic as in script 1.
# This keeps the script self-contained.
############################################################

fBuildAKMSample <- function(sPath) {
  
  dtRaw <- as.data.table(
    read_dta(
      sPath,
      col_select = c(
        id,
        year,
        firmid,
        retrib03,
        gior_r,
        log_dailywages,
        age,
        sec1
      )
    )
  )
  
  setnames(
    dtRaw,
    old = c(
      "id",
      "year",
      "firmid",
      "retrib03",
      "gior_r",
      "log_dailywages",
      "age",
      "sec1"
    ),
    new = c(
      "iWorker",
      "iYear",
      "iFirm",
      "dEarn",
      "dDays",
      "dLogWage",
      "iAge",
      "iSector"
    )
  )
  
  print(names(dtRaw))
  str(dtRaw)
  
  dtOut <- dtRaw[
    !is.na(iWorker) &
      !is.na(iYear) &
      !is.na(iFirm) &
      !is.na(dEarn) &
      !is.na(dDays) &
      !is.na(dLogWage) &
      !is.na(iAge) &
      !is.na(iSector) &
      dEarn > 0 &
      dDays > 0 &
      iAge >= 20 &
      iAge <= 60 &
      !(iSector %in% c(9, 10))
  ]
  
  setorder(dtOut, iWorker, iYear, -dEarn)
  dtOut <- unique(dtOut, by = c("iWorker", "iYear"))
  
  dtSectorMap <- data.table(
    iSector = 0:8,
    sSectorName = c(
      "Agriculture",
      "Energy_gas_water",
      "Mining_extraction",
      "Manufacturing_mechanical",
      "Manufacturing_food_clothes_furniture",
      "Construction",
      "Commerce_hospitality",
      "Telecom_transport",
      "Finance_insurance"
    )
  )
  
  dtOut <- merge(
    dtOut,
    dtSectorMap,
    by = "iSector",
    all.x = TRUE,
    all.y = FALSE
  )
  
  return(dtOut)
}

############################################################
# 2. Helper: restrict to pooled LCC from script 1
############################################################

fRestrictToPooledLCC <- function(dtIn) {
  
  if (!file.exists("AKM_pooled_LCC_worker_ids.csv.gz")) {
    stop("Run script 1 first: AKM_pooled_LCC_worker_ids.csv.gz not found.")
  }
  
  if (!file.exists("AKM_pooled_LCC_firm_ids.csv.gz")) {
    stop("Run script 1 first: AKM_pooled_LCC_firm_ids.csv.gz not found.")
  }
  
  dtWorkersLCC <- fread("AKM_pooled_LCC_worker_ids.csv.gz")
  dtFirmsLCC   <- fread("AKM_pooled_LCC_firm_ids.csv.gz")
  
  dtOut <- dtIn[
    iWorker %in% dtWorkersLCC$iWorker &
      iFirm %in% dtFirmsLCC$iFirm
  ]
  
  return(dtOut)
}

############################################################
# 3. Helper: compute AKM decomposition from an estimated model
#
# Theory:
#   log wage = worker FE + firm FE + year FE + age^2 + residual
#
# We work with the variance of net log wages:
#   y_net = log wage - observable component
#
# Important practical choice:
# We exclude the linear age term. This is in line with your
# current implementation and avoids the near-collinearity
# problem with worker FE + year FE.
############################################################

fAKMDecompFromModel <- function(dtUsed, mAKM, sModelType, iSectorValue, sSectorLabel) {
  
  lFE <- fixef(
    mAKM,
    fixef.iter = 50000,
    fixef.tol = 1e-4
  )
  
  dtWork <- copy(dtUsed)
  
  dtWork[, dWorkerFE := unname(lFE$iWorker[as.character(iWorker)])]
  dtWork[, dFirmFE   := unname(lFE$iFirm[as.character(iFirm)])]
  dtWork[, dYearFE   := unname(lFE$iYear[as.character(iYear)])]
  
  vCoef <- coef(mAKM)
  dBetaAge2 <- unname(vCoef[grepl("iAge\\^2", names(vCoef))])
  
  if (length(dBetaAge2) != 1) {
    stop("Could not recover age-squared coefficient.")
  }
  
  dtWork[, dXObs := dYearFE + dBetaAge2 * (iAge^2)]
  dtWork[, dYNet := dLogWage - dXObs]
  if (!("dResid" %in% names(dtWork))) {
    
    vResid <- resid(mAKM)
    
    if (length(vResid) != nrow(dtWork)) {
      stop("Residual length does not match dtWork. In pooled_then_sector mode, attach pooled residuals before splitting by sector.")
    }
    
    dtWork[, dResid := vResid]
  }
  
  dVarYNet   <- var(dtWork$dYNet, na.rm = TRUE)
  dVarWorker <- var(dtWork$dWorkerFE, na.rm = TRUE)
  dVarFirm   <- var(dtWork$dFirmFE, na.rm = TRUE)
  dCovWF     <- 2 * cov(dtWork$dWorkerFE, dtWork$dFirmFE, use = "complete.obs")
  dVarResid  <- var(dtWork$dResid, na.rm = TRUE)
  
  dVarExplained <- dVarWorker + dVarFirm + dCovWF + dVarResid
  
  dtOut <- data.table(
    sModelType = sModelType,
    iSector = iSectorValue,
    sSectorName = sSectorLabel,
    iNObsUsed = nrow(dtWork),
    iNWorkersUsed = uniqueN(dtWork$iWorker),
    iNFirmsUsed = uniqueN(dtWork$iFirm),
    dVarYNet = dVarYNet,
    dVarWorker = dVarWorker,
    dVarFirm = dVarFirm,
    dTwoCovWorkerFirm = dCovWF,
    dVarResidual = dVarResid,
    dShareWorker = dVarWorker / dVarYNet,
    dShareFirm = dVarFirm / dVarYNet,
    dShareSorting = dCovWF / dVarYNet,
    dShareResidual = dVarResid / dVarYNet,
    dShareTotal = dVarExplained / dVarYNet,
    dDiff = dVarYNet - dVarExplained
  )
  
  return(dtOut)
}

############################################################
# 4. Build pooled LCC sample
############################################################

dtAKM <- fBuildAKMSample(sDataPath)
dtLCC <- fRestrictToPooledLCC(dtAKM)

cat("\nSample in pooled LCC:\n")
print(
  data.table(
    iNObs = nrow(dtLCC),
    iNWorkers = uniqueN(dtLCC$iWorker),
    iNFirms = uniqueN(dtLCC$iFirm)
  )
)

cat("\nCleaned sample before LCC:\n")
print(data.table(
  iNObs = nrow(dtAKM),
  iNWorkers = uniqueN(dtAKM$iWorker),
  iNFirms = uniqueN(dtAKM$iFirm),
  iYears = uniqueN(dtAKM$iYear)
))

cat("\nPooled LCC sample:\n")
print(data.table(
  iNObs = nrow(dtLCC),
  iNWorkers = uniqueN(dtLCC$iWorker),
  iNFirms = uniqueN(dtLCC$iFirm),
  iYears = uniqueN(dtLCC$iYear)
))

cat("\nYears covered:\n")
print(range(dtAKM$iYear, na.rm = TRUE))
############################################################
# 5. Run AKM according to chosen mode
############################################################

lResults <- list()

if (sAKMMode == "pooled_then_sector") {
  
  ##########################################################
  # 5A. One pooled AKM
  ##########################################################
  
  cat("\nRunning pooled AKM ...\n")
  
  dtPool <- dtLCC[, .(iWorker, iFirm, iYear, iAge, dLogWage, iSector, sSectorName)]
  
  mPool <- feols(
    dLogWage ~ I(iAge^2) | iYear + iWorker + iFirm,
    data = dtPool,
    fixef.iter = 50000,
    mem.clean = TRUE
  )
  
  dtPoolUsed <- dtPool[obs(mPool)]
  dtPoolUsed[, dResid := resid(mPool)]
  
  
  
  ############################################################
  # Two-sector summaries aligned with the structural notebook
  #
  # Goal:
  #   Collapse the pooled AKM sample into:
  #   - Energy
  #   - Non-energy
  #
  # This gives us:
  #   1) mean log wages by the two sectors
  #   2) pooled-FE AKM decomposition by the same two sectors
  ############################################################
  
  dtTwoSector <- copy(dtPoolUsed)
  
  # Energy is sec1 == 1 in the cleaned private-sector sample
  dtTwoSector[, iSector2 := fifelse(iSector == 1L, 1L, 0L)]
  dtTwoSector[, sSector2Name := fifelse(iSector == 1L, "Energy", "Non_energy")]
  
  ############################################################
  # 1. Mean log wage by two-sector grouping
  ############################################################
  
  dtMeanLogWage2 <- dtTwoSector[
    , .(
      iNObs = .N,
      dMeanLogWage = mean(dLogWage, na.rm = TRUE),
      dSDLogWage = sd(dLogWage, na.rm = TRUE)
    ),
    by = .(iSector2, sSector2Name)
  ][order(iSector2)]
  
  cat("\nMean log wages: Energy vs Non-energy\n")
  print(dtMeanLogWage2)
  
  ############################################################
  # 2. Pooled-FE AKM decomposition by two-sector grouping
  #
  # Important:
  #   This matches the notebook logic:
  #   estimate AKM on the pooled sample first,
  #   then describe the decomposition inside Energy vs Non-energy
  #   using the pooled FE.
  ############################################################
  
  lResults2 <- list()
  
  for (iS in sort(unique(dtTwoSector$iSector2))) {
    
    dtS <- dtTwoSector[iSector2 == iS]
    
    lResults2[[as.character(iS)]] <- fAKMDecompFromModel(
      dtUsed = dtS,
      mAKM = mPool,
      sModelType = "pooled_then_two_sector",
      iSectorValue = iS,
      sSectorLabel = dtS[1, sSector2Name]
    )
  }
  
  dtAKM2 <- rbindlist(lResults2, fill = TRUE)
  
  cat("\nPooled-FE AKM decomposition: Energy vs Non-energy\n")
  print(dtAKM2)
  
  ############################################################
  # 3. Save
  ############################################################
  
  fwrite(
    dtMeanLogWage2,
    "AKM_mean_log_wage_energy_vs_nonenergy_poolAKM.csv"
  )
  
  fwrite(
    dtAKM2,
    "AKM_results_energy_vs_nonenergy_poolAKM.csv"
  )
  ##########################################################
  # 5A(i). Overall pooled decomposition
  ##########################################################
  
  lResults[["pooled_all"]] <- fAKMDecompFromModel(
    dtUsed = dtPoolUsed,
    mAKM = mPool,
    sModelType = "pooled_then_sector",
    iSectorValue = -1,
    sSectorLabel = "All_sectors_pooled"
  )
  
  ##########################################################
  # 5A(ii). Sector decomposition using pooled FE
  #
  # Important interpretation:
  # These are sector-specific decompositions based on the pooled
  # wage model, not sector-specific AKMs.
  ##########################################################
  
  vSectors <- sort(unique(dtPoolUsed$iSector))
  
  for (iS in vSectors) {
    
    dtS <- dtPoolUsed[iSector == iS]
    
    if (nrow(dtS) < iMinObsSector) next
    
    lResults[[paste0("sector_", iS)]] <- fAKMDecompFromModel(
      dtUsed = dtS,
      mAKM = mPool,
      sModelType = "pooled_then_sector",
      iSectorValue = iS,
      sSectorLabel = dtS[1, sSectorName]
    )
  }
}

if (sAKMMode == "sector_by_sector") {
  
  ##########################################################
  # 5B. Separate AKM inside each sector
  #
  # Important interpretation:
  # Here each sector gets its own worker FE / firm FE / year FE
  # estimation, but still within the pooled LCC sample.
  ##########################################################
  
  vSectors <- sort(unique(dtLCC$iSector))
  
  for (iS in vSectors) {
    
    cat("\nRunning sector AKM:", iS, "\n")
    
    dtS <- dtLCC[
      iSector == iS,
      .(iWorker, iFirm, iYear, iAge, dLogWage, iSector, sSectorName)
    ]
    
    if (nrow(dtS) < iMinObsSector) {
      cat("Skipped: too few observations\n")
      next
    }
    
    mS <- feols(
      dLogWage ~ I(iAge^2) | iYear + iWorker + iFirm,
      data = dtS,
      fixef.iter = 50000,
      mem.clean = TRUE
    )
    
    dtSUsed <- dtS[obs(mS)]
    
    lResults[[paste0("sector_", iS)]] <- fAKMDecompFromModel(
      dtUsed = dtSUsed,
      mAKM = mS,
      sModelType = "sector_by_sector",
      iSectorValue = iS,
      sSectorLabel = dtSUsed[1, sSectorName]
    )
    
    rm(dtS, mS, dtSUsed)
    gc()
  }
}

############################################################
# 6. Collect and save results
############################################################

dtResults <- rbindlist(lResults, fill = TRUE)

print(dtResults)

if (sAKMMode == "pooled_then_sector") {
  fwrite(dtResults, "AKM_results_pooled_then_sector.csv")
}

if (sAKMMode == "sector_by_sector") {
  fwrite(dtResults, "AKM_results_sector_by_sector.csv")
}

