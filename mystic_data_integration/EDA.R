# Run main_pipeline_data.R first! This is not intended as an automated process.

len(unique(clinical_db$subjid))
len(unique(clinical_db_no_SF$subjid))

ntable(clinical_db$actual_arm)
ntable(clinical_db_no_SF$actual_arm)
