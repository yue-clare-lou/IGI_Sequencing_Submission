#! /usr/bin/env python3

import argparse
import os
import glob
from subprocess import call
import numpy as np
import pandas as pd
from collections import defaultdict
from qPCR_calculations import qPCR_calculator



def IGI_Seq_Table_generator(molarity,
                            original_qPCR_Ct_dir,
                            project_name,
                            index_plate_name,
                            lib_type,
                            sample_status,
                            index_db,
                            nebindex_indexes_sheet,
                            output):

    Cdb = qPCR_calculator(original_qPCR_Ct_dir)
    # print(Cdb)

    Vdb=Cdb.copy()
    Vdb['volume for pooling (uL)']=[molarity/x for x in Vdb['Concentration  of undiluted library (nM)']]
    Vdb['total volume']=round(sum(Vdb['volume for pooling (uL)']),2)
    Vdb['nM (Final)']=[round(molarity*len(Vdb)/x,2) for x in Vdb['total volume']]

    Sdb = Cdb[['Sample']]
    Sdb.insert(0, 'Project name/ID', project_name)
    Sdb.insert(1, 'Plate number/name', index_plate_name)
    # Sdb.insert(2, 'Well Location', 'X')
    Sdb.insert(3, 'User Sample Name', ["_".join([x, str(y)]) for x, y in zip(Sdb['Project name/ID'], Sdb['Sample'])])
    Sdb = Sdb.drop(columns=['Sample'])
    Sdb['Type of Library'] = lib_type
    Sdb['Sample Status'] = sample_status
    Sdb['i7 name'] = 'X'
    Sdb['i7 Sequence'] = 'X'
    Sdb['i5 name'] = 'X'
    Sdb['i5 Sequence'] = 'X'


    Idb = pd.read_csv(index_db)
    Sdb2=Sdb.merge(Idb,on='User Sample Name',how='left')
    # print(Sdb2)
    Sdb2.insert(2,'Well Location',Sdb2.pop('Well Location'))
    Sdb2.insert(3, 'Well number', [int(x[1:]) for x in Sdb2['Well Location']])
    Sdb2.insert(3, 'Well Letter', [x[0] for x in Sdb2['Well Location']])
    Sdb2 = Sdb2.sort_values(['Well Letter', 'Well number']).reset_index(drop=True)

    IndexDB = pd.read_csv(nebindex_indexes_sheet)
    IndexDB2 = IndexDB.rename(columns={IndexDB.filter(regex='Forward').columns[0]: 'i5 Sequence FWD', \
                                       IndexDB.filter(regex='Reverse').columns[0]: 'i5 Sequence REV'})

    DIdb = Sdb2.drop(columns=['i7 name', 'i7 Sequence', 'i5 name', 'i5 Sequence']).merge(IndexDB2, on=['Well Letter',
                                                                                                       'Well number'], \
                                                                                         how='left')
    DIdb = DIdb.rename(columns={'I7_Index_ID': 'i7 name', 'i7_index': 'i7 Sequence', \
                                'I5_Index_ID': 'i5 name'})
    DIdb.insert(6, 'User Sample Name num', [int(x.split("_")[-1]) for x in DIdb['User Sample Name']])
    DIdb = DIdb.sort_values('User Sample Name num').reset_index(drop=True)

    DIdb_for_submission = DIdb.drop(columns=['Well Letter', 'Well number', 'User Sample Name num'])
    DIdb_for_submission = DIdb_for_submission.rename(columns={'i5 Sequence FWD': 'i5 Sequence'})
    DIdb_for_submission = DIdb_for_submission.drop(columns=['i5 Sequence REV'])
    DIdb_for_submission.to_csv(output, index=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,add_help=False,
                                     description='Generating the table required by the IGI core sequencing facility')

    # General Arguments
    GenArgs = parser.add_argument_group('GENERAL ARGUMENTS')
    GenArgs.add_argument('-h', action="help",help="show this help message and exit")
    GenArgs.add_argument('-m','--molarity', default=25, help='molarity for pooling libraries')
    GenArgs.add_argument('-db', '--original_qPCR_Ct_dir', default='', help='the dataframe that contains your library concentrations in nM')
    GenArgs.add_argument('-pn', '--project_name', default='', help='project name (required by IGI core sequencer)')
    GenArgs.add_argument('-ip', '--index_plate', default='NEBNext E6441', help='the plate where i5/i7 were taken from')
    GenArgs.add_argument('-l', '--library_type', default='Metagenomic', help='the type of libraries (required by IGI core sequencer)')
    GenArgs.add_argument('-s', '--sample_status', default='Pooled Library',help='sample status (required by IGI core sequencer)')
    GenArgs.add_argument('-idb', '--index_db', default='sheets/sample2indices.csv',help='2-column sample and indices table (column names: User Sample Name and Well Location)')
    GenArgs.add_argument('-n', '--nebindex_indexes_sheet', default='sheets/NEBNEXT_indexes_Sheet1.csv',help='NEBNext indexes sheet')
    GenArgs.add_argument('-o','--outfile', default='NGS_Library_Import_Form.csv', help='Output file')
    args = parser.parse_args()

    IGI_Seq_Table_generator(args.molarity,
                            args.original_qPCR_Ct_dir,
                            args.project_name,
                            args.index_plate,
                            args.library_type,
                            args.sample_status,
                            args.index_db,
                            args.nebindex_indexes_sheet,
                            args.outfile
                            )



