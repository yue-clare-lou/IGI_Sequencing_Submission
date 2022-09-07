# IGI_Sequencing_Submission
#### This python script helps generate the template form that is needed for submission

```
usage: IGISeqCoreSubmissionTable.py [-h] [-m MOLARITY] [-db ORIGINAL_QPCR_CT_DIR] [-pn PROJECT_NAME] [-ip INDEX_PLATE] [-l LIBRARY_TYPE] [-s SAMPLE_STATUS] [-idb INDEX_DB]
                                    [-n NEBINDEX_INDEXES_SHEET] [-o OUTFILE] [-pl | --pooling_only | --no-pooling_only]

Generating the table required by the IGI core sequencing facility

GENERAL ARGUMENTS:
  -h                    show this help message and exit
  -m MOLARITY, --molarity MOLARITY
                        molarity for pooling libraries (default: 25)
  -db ORIGINAL_QPCR_CT_DIR, --original_qPCR_Ct_dir ORIGINAL_QPCR_CT_DIR
                        the dataframe that contains your library concentrations in nM (default: )
  -pn PROJECT_NAME, --project_name PROJECT_NAME
                        project name (required by IGI core sequencer) (default: )
  -ip INDEX_PLATE, --index_plate INDEX_PLATE
                        the plate where i5/i7 were taken from (default: NEBNext E6441)
  -l LIBRARY_TYPE, --library_type LIBRARY_TYPE
                        the type of libraries (required by IGI core sequencer) (default: Metagenomic)
  -s SAMPLE_STATUS, --sample_status SAMPLE_STATUS
                        sample status (required by IGI core sequencer) (default: Pooled Library)
  -idb INDEX_DB, --index_db INDEX_DB
                        2-column sample and indices table (column names: User Sample Name and Well Location) (default: sheets/sample2indices.csv)
  -n NEBINDEX_INDEXES_SHEET, --nebindex_indexes_sheet NEBINDEX_INDEXES_SHEET
                        NEBNext indexes sheet (default: sheets/NEBNEXT_indexes_Sheet1.csv)
  -o OUTFILE, --outfile OUTFILE
                        Output file (default: ./NGS_Library_Import_Form.csv)
  -pl, --pooling_only, --no-pooling_only
                        If True, then skip generating the submission table (default: None)
                        
```
