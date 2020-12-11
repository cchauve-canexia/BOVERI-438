'''
Collecting data to compare the indel results of strelka on a list of NextSeq
runs with the results on the matching MiSeq runs.
Identifying the dump variants file and, if none, the output and raw data files
for each pair of matching runs.
Computing missing dump files.
Generating a CSV file of matched and unmatched  indel calls
'''

# Standard imports
from collections import defaultdict
import os
import subprocess
import sys

# Third-party imports
import pandas as pd
import numpy as np

# Local imports
import utils
from data_format import *

# Computing the latest runs files for all runs of matching pairs
def compute_latest_runs(parameters, log):
    MANIFESTS = list(utils.get_manifests(parameters).keys())
    # Data files, indexed by manifest
    NEXTSEQ_RUNS = utils.get_NextSeq_runs_file(parameters)
    MISEQ_RUNS = utils.get_MiSeq_runs_file(parameters)
    # Results files, indexed by manifest
    NEXTSEQ_LATEST_RUNS = utils.get_NextSeq_latest_runs_file(parameters)
    MISEQ_LATEST_RUNS = utils.get_MiSeq_latest_runs_file(parameters)
    # Command to generate result files
    LATEST_RUNS_CMD = utils.get_latest_run_cmd(parameters)
    for manifest in MANIFESTS:
        # NextSeq run
        nextseq_in, nextseq_out = NEXTSEQ_RUNS[manifest], NEXTSEQ_LATEST_RUNS[manifest]
        latest_runs_cmd = utils.get_latest_run_cmd(parameters, nextseq_in, nextseq_out)
        log.write('LATEST_RUNS\t' + ' '.join(latest_runs_cmd) + '\n')
        subprocess.call(latest_runs_cmd)
        # MiSeq run
        miseq_in, miseq_out = MISEQ_RUNS[manifest], MISEQ_LATEST_RUNS[manifest]
        latest_runs_cmd = utils.get_latest_run_cmd(parameters, miseq_in, miseq_out)
        log.write('LATEST_RUNS\t' + ' '.join(latest_runs_cmd) + '\n')
        subprocess.call(latest_runs_cmd)

# Computing the dump files
def compute_dump_files(parameters, log):
    MANIFESTS = list(utils.get_manifests(parameters).keys())
    # Input files indicating directory of latest run
    NEXTSEQ_LATEST_RUNS = utils.get_NextSeq_latest_runs_file(parameters)
    MISEQ_LATEST_RUNS = utils.get_MiSeq_latest_runs_file(parameters)
    # Output files
    DUMP_FILES = open(utils.get_dump_tsv_file(parameters), 'w')
    DUMP_FILES.write(DUMP_HEADER)
    # List of runs for which the dump file will be generated
    RUNS_TO_DO = []
    for manifest in MANIFESTS:
        # NextSeq runs
        for _, run_data in pd.read_csv(NEXTSEQ_LATEST_RUNS[manifest], sep=',').iterrows():
            RUNS_TO_DO.append(run_data)
        # MiSeq runs
        for _, run_data in pd.read_csv(MISEQ_LATEST_RUNS[manifest], sep=',').iterrows():
            RUNS_TO_DO.append(run_data)
        # Generating dump files
        for run_data in RUNS_TO_DO:
            run, latest_run = run_data[RUN_ID], utils.get_latest_run_file(run_data)
            if pd.isnull(latest_run): dump_cmd, dump_out = 'no_latest_run', str(np.nan)
            else:
                dump_out = utils.get_dump_out(parameters, run)
                if not os.path.isfile(dump_out):
                    dump_cmd = utils.get_dump_cmd(parameters, run, latest_run)
                    subprocess.call(dump_cmd)
                else:
                    dump_cmd = 'file exists'
            DUMP_FILES.write('\n' + run + '\t' + dump_out)
            log.write('DUMP_FILES\t' + run + '\t' + dump_cmd + '\t' + dump_out + '\n')
    DUMP_FILES.close()

# Extracting matched variants into a csv file
def compute_variants(parameters, log):

    def get_run_pairs(runs_df):
        '''Extract all run pairs to consider from dataframe of runs'''
        RUN_PAIRS =  {}
        for _, run_data in runs_df.iterrows():
            # For each pair, both orders are recorded
            RUN_PAIRS[run_data[RUN_ID]] = run_data[MISEQ_RUN]
            RUN_PAIRS[run_data[MISEQ_RUN]] = run_data[RUN_ID]
        return RUN_PAIRS

    def get_excluded_samples(excluded_samples_df, run_pairs, manifest):
        '''Get excluded samples for a manifest and a list of run pairs, indexed
        by run pair'''
        EXCLUDED_SAMPLES = {}
        RUNS_LIST = list(run_pairs.keys())
        for run_1, run_2 in run_pairs.items():
            # By default, no excluded sample per run pair
            EXCLUDED_SAMPLES[(run_1, run_2)] = []
        for _, sample_data in excluded_samples_df.iterrows():
            # Sample = integer DNA-19157-CG001Qv40Run179-4 -> 4
            run, sample = sample_data[EXC_RUN], sample_data[EXC_SAMPLE].split('-')[-1]
            if manifest in run and run in RUNS_LIST:
                EXCLUDED_SAMPLES[(RUN_PAIRS[run], run)].append(sample)
                EXCLUDED_SAMPLES[(run, RUN_PAIRS[run])].append(sample)
        return EXCLUDED_SAMPLES

    def get_MSI_coordinates(manifest_file):
        '''Get a dictionary indexed by chromosome of a list (start,end) of amplicons
        coordinates for amplicons labeled as MSI amplicons'''
        amplicons_df = pd.read_csv(manifest_file, sep='\t')
        MSI_df = amplicons_df.loc[amplicons_df[AMP_MUT]==AMP_MSI]
        MSI_coords = defaultdict(list)
        for _, amplicon in MSI_df.iterrows():
            MSI_coords[amplicon[AMP_CHR]].append((int(amplicon[AMP_START]), int(amplicon[AMP_END])))
        return MSI_coords

    def test_out_of_MSI(variant, MSI_coords):
        '''True if variant not in an MSI amplicon'''
        chromosome, position = variant[DUMP_CHR], variant[DUMP_POS]
        if chromosome in MSI_coords.keys():
            for (a_start, a_end) in MSI_coords[chromosome]:
                if a_start <= position <= a_end: return False
        return True

    def add_variant(variants_dict, variant, run_pair, run_type, excluded_samples, MSI_coords):
        '''Add variant to variants_dict if it is in the target mutations (indels)
        not in an excluded sample and not in an MSI amplicon'''
        if variant[DUMP_MUT] in DUMP_MUT_TO_SELECT:
            run_type_ext = {'MiSeq': '_m', 'NextSeq': '_n'}
            ext = run_type_ext[run_type]
            nextseq_run, miseq_run = run_pair['NextSeq'], run_pair['MiSeq']
            sample_run = variant[DUMP_SAMPLE].replace(DUMP_INDEL_SUFFIX, '')
            sample_1 = sample_run.split('-CG001')[0]
            sample_2 = sample_run.split('_S')[1]
            test_not_excluded_sample = sample_2 not in excluded_samples[(nextseq_run, miseq_run)]
            test_not_in_MSI = test_out_of_MSI(variant, MSI_coords)
            if test_not_excluded_sample and test_not_in_MSI:
                sample = sample_1 + '_' + sample_2
                index = '_'.join([variant[DUMP_CHR], str(variant[DUMP_POS]), variant[DUMP_REF], variant[DUMP_ALT]])
                v_id = '_'.join([nextseq_run, miseq_run, sample, index])
                variants_dict[v_id]['index'] = index
                variants_dict[v_id][run_type + '_sample'] = sample_run
                variants_dict[v_id]['gene'] = variant[DUMP_GENE]
                variants_dict[v_id]['codon'] = variant[DUMP_CODON]
                variants_dict[v_id]['cDNA_change'] = variant[DUMP_CDNA]
                variants_dict[v_id]['MiSeq_run'] = miseq_run
                variants_dict[v_id]['NextSeq_run'] = nextseq_run
                variants_dict[v_id]['sample'] = sample
                variants_dict[v_id]['vaf' + ext] = variant[DUMP_VAF]
                variants_dict[v_id]['coverage' + ext] = variant[DUMP_COV]
                variants_dict[v_id]['score' + ext] = variant[DUMP_SCORE]
                return 0
            else: return 1
        else: return 0

    def add_variants_run_pair(run_data, dump_files_dict, excluded_samples, MSI_coords, variants_dict, log):
        '''Add to variants_dict all variants for run pair run_pair'''
        nextseq_run, miseq_run = run_data[RUN_ID], run_data[MISEQ_RUN]
        nextseq_dump, miseq_dump = dump_files_dict[nextseq_run], dump_files_dict[miseq_run]
        if pd.isnull(nextseq_dump) or pd.isnull(miseq_dump):
            log.write('VARIANTS\t' + ','.join([nextseq_run, miseq_run ]) + '\t' + 'missing dump file\n')
        else:
            nb_excluded_variants = 0
            run_pair = {'NextSeq': nextseq_run, 'MiSeq': miseq_run}
            for _, variant in pd.read_csv(nextseq_dump, sep='\t').iterrows():
                nb_excluded_variants += add_variant(variants_dict, variant, run_pair, 'NextSeq', excluded_samples,  MSI_coords)
            for _, variant in pd.read_csv(miseq_dump, sep='\t').iterrows():
                nb_excluded_variants += add_variant(variants_dict, variant, run_pair, 'MiSeq', excluded_samples,  MSI_coords)
            log.write('VARIANTS\t' + ','.join([nextseq_run, miseq_run ]) + '\t' + str(nb_excluded_variants) + 'excluded variants\n')

    MANIFESTS_INFO = utils.get_manifests(parameters)
    DUMP_FILES_DF = pd.read_csv(utils.get_dump_tsv_file(parameters), sep='\t')
    DUMP_FILES_DICT = {x[DUMP_RUN]: x[DUMP_FILE] for _, x in DUMP_FILES_DF.iterrows()}
    NEXTSEQ_LATEST_RUNS_FILES = utils.get_NextSeq_latest_runs_file(parameters)
    EXCLUDED_SAMPLES_FILES = utils.get_excluded_samples_file(parameters)
    for manifest, manifest_file in MANIFESTS_INFO.items():
        MSI_COORDS = get_MSI_coordinates(manifest_file)
        RUNS_DF = pd.read_csv(NEXTSEQ_LATEST_RUNS_FILES[manifest], sep=',')
        RUN_PAIRS = get_run_pairs(RUNS_DF)
        EXCLUDED_SAMPLES_DF = pd.read_csv(EXCLUDED_SAMPLES_FILES[manifest], sep=',')
        EXCLUDED_SAMPLES = get_excluded_samples(EXCLUDED_SAMPLES_DF, RUN_PAIRS, manifest)
        VARIANTS_DICT = defaultdict(dict)
        for _, run_data in RUNS_DF.iterrows():
            add_variants_run_pair(run_data, DUMP_FILES_DICT, EXCLUDED_SAMPLES, MSI_COORDS, VARIANTS_DICT, log)
        VARIANTS_DF = pd.DataFrame.from_dict(VARIANTS_DICT, orient='index')
        VARIANTS_DF.to_csv(utils.get_variants_csv_file(parameters, manifest), sep=',')

# Reading parameters
PARAMETERS = utils.read_parameters(sys.argv[1])
# Log file
LOG = open(sys.argv[2], 'w')
# Commands to run
CMDS = sys.argv[3:]

if 'latest_runs' in CMDS:
    compute_latest_runs(PARAMETERS, LOG)
if 'dump_files' in CMDS:
    compute_dump_files(PARAMETERS, LOG)
if 'variants' in CMDS:
    compute_variants(PARAMETERS, LOG)

LOG.close()
