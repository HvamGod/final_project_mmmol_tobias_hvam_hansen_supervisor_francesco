#### MMoL project: pooled largest connected component
rm(list = ls())

library(haven)
library(data.table)
library(igraph)

############################################################
# 0. User settings
############################################################

sDataPath <- "Veneto_panel_1984_2001.dta"

############################################################
# 1. Load only variables needed for:
#    - sample construction
#    - sector description
#    - worker-firm mobility graph
#
# Important empirical choice:
# We drop sectors 9 and 10 now, because they will not be part
# of the project going forward.
############################################################

dtRaw <- as.data.table(
  read_dta(
    sDataPath,
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

############################################################
# 2. Clean sample
#
# We keep the same baseline AKM sample restrictions as in your
# current code, so that the connected set is computed on the
# exact sample used later in AKM.
############################################################

dtAKM <- dtRaw[
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

rm(dtRaw)
gc()

############################################################
# 3. Keep main job per worker-year
#
# Important modeling choice:
# The raw data are job-level. For AKM we want one employer
# per worker-year. We keep the job with the highest annual
# earnings in each worker-year.
############################################################

setorder(dtAKM, iWorker, iYear, -dEarn)
dtAKM <- unique(dtAKM, by = c("iWorker", "iYear"))

gc()

############################################################
# 4. Add sector labels
############################################################

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

dtAKM <- merge(
  dtAKM,
  dtSectorMap,
  by = "iSector",
  all.x = TRUE,
  all.y = FALSE
)

############################################################
# 5. Build worker-firm mobility graph
#
# Standard AKM logic:
# firm effects are only comparable across firms that are linked
# directly or indirectly through worker mobility.
#
# We therefore build a bipartite graph with:
# - worker nodes
# - firm nodes
# - one edge for each observed worker-firm relationship
#
# Then we keep the largest connected component.
############################################################

dtEdges <- unique(
  dtAKM[, .(iWorker, iFirm)]
)

dtEdges[, sWorkerNode := paste0("w_", iWorker)]
dtEdges[, sFirmNode   := paste0("f_", iFirm)]

gMobility <- graph_from_data_frame(
  d = dtEdges[, .(sWorkerNode, sFirmNode)],
  directed = FALSE
)

dtComp <- data.table(
  sNode = names(components(gMobility)$membership),
  iComp = as.integer(components(gMobility)$membership)
)

dtCompSize <- dtComp[, .N, by = iComp]
iLargestComp <- dtCompSize[which.max(N), iComp]

dtCompLargest <- dtComp[iComp == iLargestComp]

vWorkerNodesLCC <- dtCompLargest[grepl("^w_", sNode), sNode]
vFirmNodesLCC   <- dtCompLargest[grepl("^f_", sNode), sNode]

dtWorkersLCC <- data.table(
  iWorker = as.numeric(sub("^w_", "", vWorkerNodesLCC))
)

dtFirmsLCC <- data.table(
  iFirm = as.numeric(sub("^f_", "", vFirmNodesLCC))
)

############################################################
# 6. Restrict sample to pooled largest connected component
############################################################

dtAKM[
  ,
  bWorkerInLCC := iWorker %in% dtWorkersLCC$iWorker
]

dtAKM[
  ,
  bFirmInLCC := iFirm %in% dtFirmsLCC$iFirm
]

dtAKM[
  ,
  bInLCC := bWorkerInLCC & bFirmInLCC
]

dtLCC <- dtAKM[bInLCC == TRUE]

############################################################
# 7. Descriptive summary: overall
############################################################

dtOverall <- data.table(
  sSample = c("Before_LCC", "After_LCC"),
  iNObs = c(nrow(dtAKM), nrow(dtLCC)),
  iNWorkers = c(uniqueN(dtAKM$iWorker), uniqueN(dtLCC$iWorker)),
  iNFirms = c(uniqueN(dtAKM$iFirm), uniqueN(dtLCC$iFirm))
)

dtOverall[
  ,
  dShareObsKept := iNObs / dtOverall[sSample == "Before_LCC", iNObs]
]

dtOverall[
  ,
  dShareWorkersKept := iNWorkers / dtOverall[sSample == "Before_LCC", iNWorkers]
]

dtOverall[
  ,
  dShareFirmsKept := iNFirms / dtOverall[sSample == "Before_LCC", iNFirms]
]

############################################################
# 8. Descriptive summary: by sector
#
# We report retained shares of:
# - observations
# - workers
# - firms
#
# Note:
# This is done after pooled-LCC construction.
# Alternative approach:
# build sector-specific LCCs inside each sector separately.
#
# Pros of sector-specific LCCs:
# - tailored to sector AKMs
# - can retain more observations in some sectors
#
# Cons:
# - each sector then lives on a different mobility graph
# - comparability across sectors becomes weaker
#
# Since your immediate task is whole-sample AKM first, pooled
# LCC is the clean default.
############################################################

dtSectorBefore <- dtAKM[
  ,
  .(
    iNObs_before = .N,
    iNWorkers_before = uniqueN(iWorker),
    iNFirms_before = uniqueN(iFirm)
  ),
  by = .(iSector, sSectorName)
]

dtSectorAfter <- dtLCC[
  ,
  .(
    iNObs_after = .N,
    iNWorkers_after = uniqueN(iWorker),
    iNFirms_after = uniqueN(iFirm)
  ),
  by = .(iSector, sSectorName)
]

dtSectorSummary <- merge(
  dtSectorBefore,
  dtSectorAfter,
  by = c("iSector", "sSectorName"),
  all.x = TRUE,
  all.y = FALSE
)

for (sCol in c("iNObs_after", "iNWorkers_after", "iNFirms_after")) {
  set(
    dtSectorSummary,
    i = which(is.na(dtSectorSummary[[sCol]])),
    j = sCol,
    value = 0L
  )
}

dtSectorSummary[
  ,
  dShareObsKept := iNObs_after / iNObs_before
]

dtSectorSummary[
  ,
  dShareWorkersKept := iNWorkers_after / iNWorkers_before
]

dtSectorSummary[
  ,
  dShareFirmsKept := iNFirms_after / iNFirms_before
]

setorder(dtSectorSummary, iSector)

############################################################
# 9. Print and save results
#
# We save:
# - overall LCC description
# - sector-level LCC description
# - worker IDs in pooled LCC
# - firm IDs in pooled LCC
#
# These ID files are then used by script 2.
############################################################

cat("\nOverall pooled LCC summary:\n")
print(dtOverall)

cat("\nSector-level pooled LCC summary:\n")
print(dtSectorSummary)

fwrite(dtOverall, "AKM_pooled_LCC_summary_overall.csv")
fwrite(dtSectorSummary, "AKM_pooled_LCC_summary_by_sector.csv")

fwrite(dtWorkersLCC, "AKM_pooled_LCC_worker_ids.csv.gz")
fwrite(dtFirmsLCC, "AKM_pooled_LCC_firm_ids.csv.gz")
