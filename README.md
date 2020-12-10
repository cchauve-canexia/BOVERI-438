# BOVERI-438: Investigate performance of Strelka in NextSeq data
Comparison of the strelka results between matching NextSeq and MiSeq runs.

A list of matching NextSeq/MiSeq runs for amplicons manifests v4.0 and v5.1 was obtained from <a href="https://docs.google.com/spreadsheets/d/1vIIKf5DTvHQy_mG7rCq0oCy-z82kj7bmv88wYcsWWCo/edit#gid=779463918">RUN_SUMMARIES</a>, on December 9, 2020, using <a href="https://github.com/contextual-genomics/Bioinformatics/blob/dev/Operations/data_to_exclude/development_runs/v5_miseq/dev_RD_excluded_runs.tsv">v5_miseq_RD_excluded_runs</a> and <a href="https://github.com/contextual-genomics/Bioinformatics/blob/dev/longitudinal_analysis/cohort_details/master_RD_excluded_RUNS.csv">v4_RD_excluded_RUNS</a> to exclude runs. This resulted in four files:
- ```data/RUN_SUMMARIES_NextSeq_v40.csv and data/RUN_SUMMARIES_NextSeq_v51.csv```: information on NextSeq runs, including matching MiSeq run, per manifest;
- ```data/RUN_SUMMARIES_MiSeq_v40.csv and data/RUN_SUMMARIES_MiSeq_v51.csv```: information on MiSeq runs.

Lists of excluded samples were obtained from <a href="https://github.com/contextual-genomics/Bioinformatics/blob/dev/Operations/data_to_exclude/development_runs/v5_miseq/dev_RD_failed_samples.tsv">v5_miseq_RD_failed_samples</a> and <a href="https://github.com/contextual-genomics/Bioinformatics/blob/dev/longitudinal_analysis/cohort_details/master_RD_excluded_samples.csv">v4_RD_excluded_samples</a>. This resulted in two files, ```data/EXCLUDED_SAMPLES_v40.csv``` and ```EXCLUDED_SAMPLES_v51.csv```.

For each pair of NextSeq/MiSe matching runs, the latest_run information was obtained using the script <a href="https://github.com/contextual-genomics/Bioinformatics/blob/master/Operations/get_pipeline_latest_output.py">get_pipeline_latest_output.py</a>. This was done through the command  
```python bin/main.py parameters.tsv log/latest_runs_v40_v51.log latest_runs > latest_runs_v40_v51.out 2>&1```

The variants dump files were obtained using a fixed version of the script <a href="https://github.com/contextual-genomics/biosys/blob/rd/rd_analysis/dump_variants.py">dump_variants.py</a>. 
This was done through the command  
```python bin/main.py parameters.tsv log/dump_files_v40_v51.log dump_files > dump_files_v40_v51.out 2>&1```

