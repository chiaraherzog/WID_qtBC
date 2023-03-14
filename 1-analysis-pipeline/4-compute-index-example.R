# Compute index in datasets for evaluation: -----------------------
# indices are computed for each dataset using the following code:
source("0-source-code/index_mQTL.R")
wid_qtbc <- index_mQTL(beta)

# The following datasets are used:
# EGAD00010002079 and EGAD00010002081 (Validation set)
# EGAD00010002073 and GSE133985 (type, "tissue:ch1" = ipsilateral-normal (normal adjacent)), (Breast tissue methylation set)
# 
# 

