
############################################################
# 0. Check that both result files exist
############################################################

if (!file.exists("AKM_results_pooled_then_sector.csv")) {
  stop("Missing file: AKM_results_pooled_then_sector.csv")
}

if (!file.exists("AKM_results_sector_by_sector.csv")) {
  stop("Missing file: AKM_results_sector_by_sector.csv")
}

############################################################
# 1. Load results
############################################################

dtPool <- fread("AKM_results_pooled_then_sector.csv")
dtSec  <- fread("AKM_results_sector_by_sector.csv")

############################################################
# 2. Keep sector rows only
#
# In pooled_then_sector mode, iSector = -1 is the full pooled
# decomposition. We exclude that row here because the comparison
# is sector-by-sector.
############################################################

dtPool <- dtPool[iSector >= 0]
dtSec  <- dtSec[iSector >= 0]

############################################################
# 3. Merge the two result tables
############################################################

dtCompare <- merge(
  dtPool,
  dtSec,
  by = c("iSector", "sSectorName"),
  suffixes = c("_pool", "_sector")
)

############################################################
# 4. Compute differences
#
# Positive difference means:
#   sector-by-sector estimate > pooled-then-sector estimate
############################################################

dtCompare[
  ,
  `:=`(
    dDiffVarYNet = dVarYNet_sector - dVarYNet_pool,
    dDiffVarWorker = dVarWorker_sector - dVarWorker_pool,
    dDiffVarFirm = dVarFirm_sector - dVarFirm_pool,
    dDiffVarSorting = dTwoCovWorkerFirm_sector - dTwoCovWorkerFirm_pool,
    dDiffVarResidual = dVarResidual_sector - dVarResidual_pool,
    dDiffShareWorker = dShareWorker_sector - dShareWorker_pool,
    dDiffShareFirm = dShareFirm_sector - dShareFirm_pool,
    dDiffShareSorting = dShareSorting_sector - dShareSorting_pool,
    dDiffShareResidual = dShareResidual_sector - dShareResidual_pool
  )
]

############################################################
# 5. Keep a compact comparison table
############################################################

dtCompareOut <- dtCompare[
  ,
  .(
    iSector,
    sSectorName,
    
    iNObsUsed_pool,
    iNObsUsed_sector,
    iNWorkersUsed_pool,
    iNWorkersUsed_sector,
    iNFirmsUsed_pool,
    iNFirmsUsed_sector,
    
    dVarYNet_pool,
    dVarYNet_sector,
    dDiffVarYNet,
    
    dShareWorker_pool,
    dShareWorker_sector,
    dDiffShareWorker,
    
    dShareFirm_pool,
    dShareFirm_sector,
    dDiffShareFirm,
    
    dShareSorting_pool,
    dShareSorting_sector,
    dDiffShareSorting,
    
    dShareResidual_pool,
    dShareResidual_sector,
    dDiffShareResidual
  )
]

############################################################
# 6. Print and save
############################################################

cat("\nComparison of pooled-then-sector vs sector-by-sector AKM:\n")
print(dtCompareOut)

fwrite(dtCompareOut, "AKM_comparison_pooled_vs_sector.csv")
