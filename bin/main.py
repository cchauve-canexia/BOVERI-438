"""
Collecting data to compare the indel results of strelka on a list of NextSeq
runs with the results on the matching MiSeq runs.
Identifying the dump variants file and, if none, the output and raw data files
for each pair of matching runs.
Computing missing dump files.
Generating files with all reportable indel calls, one CSV file and one TSV file
where variants are matched
"""

# Standard imports
from collections import defaultdict
import os
import subprocess
import sys

# Third-party imports
import pandas as pd
import numpy as np

# Reading parameters
PARAMETERS_FILE = open(sys.argv[1], 'r').readlines()
PARAMETERS = {}
for l in PARAMETERS_FILE:
    if len(l) > 2:
        l1 = l.rstrip().split('\t')
        PARAMETERS[l1[0]] = l1[1]
# Log file
LOG_FILE = sys.argv[2]
# Commands to run
CMDS = sys.argv[3:]

# Parameters description

'''
Log file
'''
LOG = open(LOG_FILE, 'w')
'''
File linking matching NextSeq and MiSeq runs, generated from RUN_SUMMARIES, FORMAT (tsv):
NextSeq_run: NextSeq run ID
MiSeq_run: matching MiSeq run ID
'''
NEXTSEQ_2_MISEQ_TSV = PARAMETERS['NEXTSEQ_2_MISEQ_TSV']
'''
Files providing data information for MiSeq/NextSeq runs without latest run, FORMAT (tsv):
Run#: run ID
Raw Data: raw data ID
Pipeline Output in S3 bucket rd-output: latest run from RUN_SUMMARIES
'''
NEXTSEQ_DATA_TSV = PARAMETERS['NEXTSEQ_DATA_TSV']
MISEQ_DATA_TSV = PARAMETERS['MISEQ_DATA_TSV']
'''
Command to generate latest runs information
'''
LATEST_RUNS_CMD = PARAMETERS['LATEST_RUNS_CMD']
'''
Files providing data information for MiSeq/NextSeq runs with latest, FORMAT (tsv):
Run#: run ID
Raw Data: raw data ID
Pipeline Output in S3 bucket rd-output: latest run from RUN_SUMMARIES
latest_run: latest run from script
'''
NEXTSEQ_LATEST_RUNS_TSV = PARAMETERS['NEXTSEQ_LATEST_RUNS_TSV']
MISEQ_LATEST_RUNS_TSV = PARAMETERS['MISEQ_LATEST_RUNS_TSV']
'''
File aggregating all information about matched runs, FORMAT (tsv):
NextSeq_run: NextSeq run ID
NextSeq_raw: raw data for NextSeq run
NextSeq_latest_RUN: latest run for NextSeq from RUN_SUMMARIES
NextSeq_latest_FD: latest run for NextSeq from script
MiSeq_run: matching MiSeq run ID
MiSeq_raw: raw data for matching MiSeq run
MiSeq_latest_RUN: latest run for matching MiSeq from RUN_SUMMARIES
MiSeq_latest_FD: latest run for matching MiSeq from script
'''
MAIN_DF_FILE = PARAMETERS['MAIN_DF']
'''
DUMP_CMD: generic dumping command
DUMP_IN: input file to dumping command
DUMP_PARAMETERS: parameters to dumping command
DUMP_OUT: output file of dumping command
DUMP_LIST_TSV: list of dump files, FORMAT (tsv)
run_pair_id, run, run_type, dump_file
'''
DUMP_CMD = PARAMETERS['DUMP_CMD'].split(' ')
DUMP_IN = PARAMETERS['DUMP_IN']
DUMP_OUT = PARAMETERS['DUMP_OUT']
DUMP_PARAMETERS = PARAMETERS['DUMP_PARAMETERS'].split(' ')
DUMP_FILES_TSV = PARAMETERS['DUMP_FILES_TSV']
'''
Strelka reportable indel thresholds
'''
STRELKA_MIN_VAF = float(PARAMETERS['STRELKA_MIN_VAF'])
STRELKA_MIN_SCORE = int(PARAMETERS['STRELKA_MIN_SCORE'])
STRELKA_MIN_COVERAGE = int(PARAMETERS['STRELKA_MIN_COVERAGE'])
'''
CSV file where all reportable variants from Strelka are reported in dumpvariants format
'''
VARIANTS_CSV = PARAMETERS['VARIANTS_CSV']
VARIANTS_CSV_HEADER = ['sample', 'chromosome', 'pos', 'mutation_type', 'ref', 'alt']
VARIANTS_CSV_HEADER += ['coverage', 'score', 'tumour_ref_allele_count', 'tumour_alt_allele_count', 'vaf']
VARIANTS_CSV_HEADER += ['gene', 'exon', 'cdna_change', 'codon', 'transcript', 'annotation_type', 'coding', 'impact', 'filter']
'''
TSV  file where all reportable variants from Strelka are reported, FORMAT (tsv)
v_id: variant ID
NextSeq, MiSeq: NextSeq run, MiSeq run
vaf_n, vaf_m: VAF in NextSeq/MiSeq (empty if not called)
score_n, score_m: score in NextSeq/MiSeq (empty if not called)
coverage_n, coverage_m: coverage in NextSeq/MiSeq (empty if not called)
'''
VARIANTS_TSV = PARAMETERS['VARIANTS_TSV']


# Generating the latest runs files
if 'latest_runs' in CMDS:
    latest_runs_cmd = LATEST_RUNS_CMD.replace('RUNS_IN', NEXTSEQ_DATA_TSV).replace('RUNS_OUT', NEXTSEQ_LATEST_RUNS_TSV).split(' ')
    LOG.write('LATEST_RUNS\t' + ' '.join(latest_runs_cmd) + '\n')
    subprocess.call(latest_runs_cmd)
    latest_runs_cmd = LATEST_RUNS_CMD.replace('RUNS_IN', MISEQ_DATA_TSV).replace('RUNS_OUT', MISEQ_LATEST_RUNS_TSV).split(' ')
    LOG.write('LATEST_RUNS\t' + ' '.join(latest_runs_cmd) + '\n')
    subprocess.call(latest_runs_cmd)

# Reading data files and creating a consolidated data frame of all input information
if 'aggregate' in CMDS:
    NEXTSEQ_2_MISEQ_DF = pd.read_csv(NEXTSEQ_2_MISEQ_TSV, sep='\t')
    NEXTSEQ_RUNS_DF = pd.read_csv(NEXTSEQ_LATEST_RUNS_TSV, sep='\t')
    NEXTSEQ_RUNS_DF.rename(
        columns={'Run#': 'NextSeq_run',
        'Raw Data': 'NextSeq_raw',
        'Pipeline Output in S3 bucket rd-output': 'NextSeq_latest_RUN',
        'latest_run': 'NextSeq_latest_FD'},
        inplace=True
    )
    MISEQ_RUNS_DF = pd.read_csv(MISEQ_LATEST_RUNS_TSV, sep='\t')
    MISEQ_RUNS_DF.rename(
        columns={'Run#': 'MiSeq_run',
        'Raw Data': 'MiSeq_raw',
        'Pipeline Output in S3 bucket rd-output': 'MiSeq_latest_RUN',
        'latest_run': 'MiSeq_latest_FD'},
        inplace=True
    )
    MAIN_DF = pd.merge(pd.merge(NEXTSEQ_2_MISEQ_DF, NEXTSEQ_RUNS_DF, on='NextSeq_run'), MISEQ_RUNS_DF, on='MiSeq_run')
    MAIN_DF.to_csv(MAIN_DF_FILE, sep='\t')
else:
    MAIN_DF = pd.read_csv(MAIN_DF_FILE, sep='\t')

# Creating dumpfiles
if 'dump_files' in CMDS:
    def create_dump_file(run_pair_id, latest_run, dump_files, run, run_type):
        if pd.isnull(latest_run):
            dump_files.write('\n' + str(run_pair_id) + '\t' + run + '\t' + run_type + '\t' + str(np.nan))
            LOG.write('DUMP_FILES\t' + run + '\tno_latest_run\t' + str(np.nan) + '\n')
        else:
            dump_out = DUMP_OUT.replace('RUN_NAME', run)
            if not os.path.isfile(dump_out):
                dump_cmd = DUMP_CMD + ['-d', DUMP_IN.replace('LATEST_RUN', latest_run), '-o', dump_out] + DUMP_PARAMETERS
                LOG.write('DUMP_FILES\t' + run + '\tcomputes\t' +  ' '.join(dump_cmd) + '\n')
                subprocess.call(dump_cmd)
            else:
                LOG.write('DUMP_FILES\t' + run + '\texists\t' +  dump_out + '\n')
            dump_files.write('\n' + str(run_pair_id) + '\t' + run + '\t' + run_type + '\t' + dump_out)

    def get_latest_run(data, run_type, run):
        latest_run = data[run_type + '_latest_FD']
        if run_type == 'NextSeq' and 'v51' in run: latest_run = data['NextSeq_latest_RUN']
        if pd.isnull(latest_run):
            LOG.write('DUMP_FILES\tlatest_run_FD\t' + run + '\t' + str(np.nan) + '\n')
            latest_run =  data[run_type + '_latest_RUN']
            if pd.isnull(latest_run):
                LOG.write('DUMP_FILES\tlatest_run_RUN\t' + '\t' + str(np.nan) + '\n')
        return latest_run

    DUMP_FILES = open(DUMP_FILES_TSV, 'w')
    DUMP_FILES.write('run_pair_id\trun\trun_type\tdump_file')
    for run_pair_id, data in MAIN_DF.iterrows():
        NextSeq_run = data['NextSeq_run']
        latest_run = get_latest_run(data, 'NextSeq', NextSeq_run)
        create_dump_file(run_pair_id, latest_run, DUMP_FILES, NextSeq_run, 'NextSeq')
        MiSeq_run = data['MiSeq_run']
        latest_run = get_latest_run(data, 'MiSeq', MiSeq_run)
        create_dump_file(run_pair_id, latest_run, DUMP_FILES, MiSeq_run, 'MiSeq')
    DUMP_FILES.close()

# Extracting variants to compare
if 'variants' in CMDS or 'variants_MiSeq' in CMDS or 'variants_NextSeq' in CMDS or 'variants_both' in CMDS:
    if 'variants_MiSeq' in CMDS: 
        VARIANTS_CSV = VARIANTS_CSV.replace('.csv', '_MiSeq.csv')
        VARIANTS_TSV = VARIANTS_TSV.replace('.tsv', '_MiSeq.csv')
    elif 'variants_NextSeq' in CMDS: 
        VARIANTS_CSV = VARIANTS_CSV.replace('.csv', '_NextSeq.csv')
        VARIANTS_TSV = VARIANTS_TSV.replace('.tsv', '_NextSeq.csv')
    elif 'variants_both' in CMDS: 
        VARIANTS_CSV = VARIANTS_CSV.replace('.csv', '_MiSeq_NextSeq.csv')
        VARIANTS_TSV = VARIANTS_TSV.replace('.tsv', '_MiSeq_NextSeq.csv')

    VARIANTS_CSV_OUT = open(VARIANTS_CSV, 'w')
    VARIANTS_CSV_OUT.write(','.join(VARIANTS_CSV_HEADER))
    DUMP_FILES_DF = pd.read_csv(DUMP_FILES_TSV, sep='\t')

    # ALl indel calls per run pair id, indexed per run type for a run pair
    INDELS_DF, RUN_PAIRS = defaultdict(dict), defaultdict(dict)
    for _, data in DUMP_FILES_DF.iterrows():
        RUN_PAIRS[data['run_pair_id']][data['run_type']] = data['run']
        if pd.notnull(data['dump_file']):
            LOG.write('VARIANTS\t' + str(data['run_pair_id']) + '\t' + data['run'] + '\t' + data['dump_file'] + '\n')
            MUTATIONS_DF_ALL = pd.read_csv(data['dump_file'], sep='\t')
            INDELS_DF[data['run_pair_id']][data['run_type']] = MUTATIONS_DF_ALL.loc[MUTATIONS_DF_ALL['mutation_type'] == 'indel']
        else:
            LOG.write('VARIANTS\t' + str(data['run_pair_id']) + '\t' + data['run'] + '\tMissing dump file\n')
            INDELS_DF[data['run_pair_id']][data['run_type']] = None

    def add_variant(variants_dict, variants_csv, variant, run_pair, run_type):
        run_type_ext = {'MiSeq': '_m', 'NextSeq': '_n'}
        ext = run_type_ext[run_type]
        sample_1 = variant['sample'].split('-CG001')[0]
        sample_2 = variant['sample'].split('_S')[1].split('_somatic')[0]
        sample = sample_1 + '_' + sample_2 + '_S' + sample_2
        v_id = '_'.join([run_pair['NextSeq'], run_pair['MiSeq'], sample, variant['chromosome'], str(variant['pos']), variant['ref'], variant['alt']])
        variants_dict[v_id]['sample'] = sample
        variants_dict[v_id]['MiSeq'] = run_pair['MiSeq']
        variants_dict[v_id]['NextSeq'] = run_pair['NextSeq']
        variants_dict[v_id]['vaf' + ext] = variant['vaf']
        variants_dict[v_id]['coverage' + ext] = variant['coverage']
        variants_dict[v_id]['score' + ext] = variant['score']
        variants_csv.write('\n' + ','.join([str(x) for x in variant]))

    def extract_MiSeq_variants(variants_dict):
        variants_to_remove = []
        for v_id, variant in variants_dict.items():
            if 'vaf_m' not in variant.keys(): variants_to_remove.append(v_id)
        for v_id in variants_to_remove:
            LOG.write('VARIANTS_MISEQ delete ' +' '.join([str(x) for x in variants_dict[v_id]]))
            del variants_dict[v_id]
            
    def extract_NextSeq_variants(variants_dict):
        variants_to_remove = []
        for v_id, variant in variants_dict.items():
            if 'vaf_n' not in variant.keys(): variants_to_remove.append(v_id)
        for v_id in variants_to_remove:
            LOG.write('VARIANTS_MISEQ delete ' +' '.join([str(x) for x in variants_dict[v_id]]))
            del variants_dict[v_id]
            

    def filter_variants_thresholds(variants_dict):
        variants_to_remove = []
        for v_id, variant in variants_dict.items():
            test_m = 'vaf_m' not in variant.keys() or variant['vaf_m'] < STRELKA_MIN_VAF or variant['score_m'] < STRELKA_MIN_SCORE or variant['coverage_m'] < STRELKA_MIN_COVERAGE
            test_n = 'vaf_n' not in variant.keys() or variant['vaf_n'] < STRELKA_MIN_VAF or variant['score_n'] < STRELKA_MIN_SCORE or variant['coverage_n'] < STRELKA_MIN_COVERAGE
            if test_m and test_n:
                variants_to_remove.append(v_id)
        for v_id in variants_to_remove:
            LOG.write('VARIANTS delete ' +' '.join([str(x) for x in variants_dict[v_id]]))
            del variants_dict[v_id]

    # Iterating over matched NextSeq/MiSeq run pairs
    VARIANTS_DICT = defaultdict(dict)
    for run_pair_id, run_pair in RUN_PAIRS.items():
        NEXTSEQ_INDELS_DF = INDELS_DF[run_pair_id]['NextSeq']
        MISEQ_INDELS_DF = INDELS_DF[run_pair_id]['MiSeq']
        if (NEXTSEQ_INDELS_DF is not None) and (MISEQ_INDELS_DF is not None):
            for _, variant in NEXTSEQ_INDELS_DF.iterrows():
                add_variant(VARIANTS_DICT, VARIANTS_CSV_OUT, variant, run_pair, 'NextSeq')
            for _, variant in MISEQ_INDELS_DF.iterrows():
                add_variant(VARIANTS_DICT, VARIANTS_CSV_OUT, variant, run_pair, 'MiSeq')
        else:
            LOG.write('VARIANTS\t' + str(run_pair_id) + '\t' + ','.join([run_pair['NextSeq'], run_pair['MiSeq']]) + '\t' + 'missing dump file\n')
    
    if 'variants_MiSeq' in CMDS:
        extract_MiSeq_variants(VARIANTS_DICT)
    elif 'variants_NextSeq' in CMDS:
        extract_NextSeq_variants(VARIANTS_DICT)
    elif 'variants_both' in CMDS:
        pass
    else:
        filter_variants_thresholds(VARIANTS_DICT)

    VARIANTS_DF = pd.DataFrame.from_dict(VARIANTS_DICT, orient='index').reset_index()
    if 'variants_MiSeq' in CMDS or 'variants_NextSeq' in CMDS or 'variants_both' in CMDS:
        VARIANTS_DF.to_csv(VARIANTS_TSV, sep=',')
    else:
        VARIANTS_DF.to_csv(VARIANTS_TSV, sep='\t')
    VARIANTS_CSV_OUT.close()

LOG.close()
