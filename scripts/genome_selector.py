# -*- coding: utf-8 -*-
#!/usr/bin/python

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Release date: 2020-04-30
# Last update: 2021-07-01

import pandas as pd
import numpy as np
from epiweeks import Week
import time
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Select list of genomes following pre-defined sampling scheme",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", required=True, help="TSV metadata file")
    parser.add_argument("--keep", required=False, help="List of samples to keep, in all instances")
    parser.add_argument("--remove", required=False, help="List of samples to remove, in all instances")
    parser.add_argument("--scheme", required=True, help="Subsampling scheme")
    parser.add_argument("--report", required=False, help="Report listing samples per category")
    args = parser.parse_args()

    metadata = args.metadata
    keep = args.keep
    remove = args.remove
    scheme = args.scheme
    report = args.report



    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_variants/nextstrain/runX_20210630_fromGISAID/'
    # import os
    # os.chdir(path)
    # metadata = path + 'pre-analyses/metadata_gisaid.tsv'
    # keep = 'config/keep.txt'
    # remove = None
    # scheme = path + 'config/subsampling_scheme.tsv'
    # report = path + 'report.tsv'
    # output = path + 'selected_strains.txt'



    today = time.strftime('%Y-%m-%d', time.gmtime())

    # force genomes to be kept in final dataset
    if keep == None:
        to_keep = []
    else:
        to_keep = [strain.strip() for strain in open(keep, 'r').readlines() if strain[0] not in ['#', '\n']]

    # subsampling scheme
    dfS = pd.read_csv(scheme, encoding='utf-8', sep='\t', dtype=str)

    results = {'strain': {}, 'gisaid_epi_isl': {}}

    ### IGNORING SAMPLES
    # list of rows to be ignored
    ignore = {}
    for idx, val in dfS.loc[dfS['purpose'] == 'ignore', 'value'].to_dict().items():
        key = dfS.iloc[idx]['filter']
        if key not in ignore.keys():
            ignore[key] = [val]
        else:
            ignore[key].append(val)
    # print(ignore)

    print('\n* Loading metadata...')
    # nextstrain metadata
    dfN = pd.read_csv(metadata, encoding='utf-8', sep='\t', dtype=str)
    dfN.fillna('', inplace=True) # replace empty values by blank
    dfN['strain'] = dfN['strain'].str.replace('hCoV-19/', '', regex=False)


    # print(dfN['strain'].tolist())
    name_accno = {}
    for genome in to_keep:
        if genome in dfN['strain'].tolist():
            accno = dfN.loc[dfN['strain'] == genome, 'gisaid_epi_isl'].values[0]
            name_accno[genome] = accno

    # drop lines if samples are set to be ignored
    for column, names in ignore.items():
        dfN = dfN[~dfN[column].isin(names)]

    print('* Removing unwanted samples, if any is listed...')
    # drop rows listed in remove.txt
    if remove == None:
        to_remove = []
    else:
        to_remove = [strain.strip() for strain in open(remove, 'r').readlines()]
    dfN = dfN[~dfN['strain'].isin(to_remove)]

    # prevent samples already selected in keep.txt from being resampled
    dfN = dfN[~dfN['strain'].isin(to_keep)]

    print('* Dropping sequences with incomplete date...')
    # drop rows with incomplete dates
    dfN = dfN[dfN['date'].apply(lambda x: len(x.split('-')) == 3)] # accept only full dates
    dfN = dfN[dfN['date'].apply(lambda x: 'X' not in x)] # exclude -XX-XX missing dates

    # convert string dates into date format
    dfN['date'] = pd.to_datetime(dfN['date']) # coverting to datetime format
    dfN = dfN.sort_values(by='date')  # sorting lines by date
    start, end = dfN['date'].min(), today

    print('* Assigning epiweek column...')
    # get epiweek end date, create column
    dfN['date'] = pd.to_datetime(dfN['date'], errors='coerce')
    dfN['epiweek'] = dfN['date'].apply(lambda x: Week.fromdate(x, system="cdc").enddate())


    ## SAMPLE FOCAL AND CONTEXTUAL SEQUENCES
    print('* Loading sampling scheme...')
    purposes = ['focus', 'context']
    subsamplers = [] # list of focal and contextual categories
    for category in purposes:
        query = {} # create a dict for each 'purpose'
        for idx, val in dfS.loc[dfS['purpose'] == category, 'value'].to_dict().items():
            key = dfS.iloc[idx]['filter']
            if key not in query.keys():
                query[key] = [val]
            else:
                query[key].append(val)
            if query not in subsamplers:
                subsamplers.append(query)

    print('* Performing genome selection...')
    # perform subsampling
    for scheme_dict in subsamplers:
        # print(scheme_dict)
        # group by level
        for level, names in scheme_dict.items():
            # create dictionary for level
            for id in results.keys():
                if level not in results[id]:
                    results[id][level] = {}

            for id in results.keys():
                for name in names:
                    # add place name as key in its corresponding level in dict
                    if name not in results[id][level].keys():
                        results[id][level][name] = []

            glevel = dfN.groupby(level)
            for name, dfLevel in glevel:
                if name in names: # check if name is among focal places
                    # apply secondary filtering, if any exists
                    if dfS.loc[dfS['value'] == name, 'filter2'].values[0] not in [None, np.nan]:
                        filter2 = dfS.loc[dfS['value'] == name, 'filter2'].values[0]
                        value2 = dfS.loc[dfS['filter2'] == filter2, 'value2'].values[0]
                        dfLevel = dfLevel[dfLevel[filter2].isin([value2])]

                    # define new temporal boundaries, if provided
                    min_date, max_date = start, end
                    new_start = dfS.loc[dfS['value'] == name, 'start'].values[0]
                    new_end = dfS.loc[dfS['value'] == name, 'end'].values[0]
                    if not pd.isna(new_start):
                        min_date = new_start
                    if not pd.isna(new_end):
                        max_date = new_end

                    # drop any row with dates outside the start/end dates
                    mask = (dfLevel['date'] > min_date) & (dfLevel['date'] <= max_date)
                    dfLevel = dfLevel.loc[mask]  # apply mask

                    gEpiweek = dfLevel.groupby('epiweek')
                    # print(dfS.loc[dfS['value'] == name, 'size'])
                    sample_size = int(dfS.loc[dfS['value'] == name, 'sample_size'])
                    total_genomes = dfLevel[level].count()
                    for epiweek, dfEpiweek in gEpiweek:
                        bin_pool = dfEpiweek['epiweek'].count() # genomes in bin
                        sampled = int(np.ceil((bin_pool/total_genomes) * sample_size)) # proportion sampled from bin
                        # print(epiweek, '-', sampled, '/', bin_pool)
                        if sampled > bin_pool: # if # requested samples higher than available genomes, get all
                            sampled = bin_pool

                        # selector
                        random_subset = dfEpiweek.sample(n=sampled)
                        for id in results.keys():
                            selected = random_subset[id].to_list()
                            results[id][level][name] = results[id][level][name] + selected

                    # drop pre-selected samples to prevent duplicates
                    dfN = dfN[~dfN[level].isin([name])]

    ### EXPORT RESULTS
    print('\n\n# Genomes sampled per category in subsampling scheme\n')
    exported = []

    outfile_names = open('selected_names.txt', 'w')
    outfile_names.write('# Genomes selected on ' + today + '\n')
    outfile_accno = open('selected_accession.txt', 'w')
    outfile_accno.write('# Genomes selected on ' + today + '\n')


    if report != None:
        outfile2 = open(report, 'w')
        outfile2.write('sample_size' + '\t' + 'category' + '\n')

    # export list selected genomes
    reported = []
    genome_count = ''
    for id in results.keys():
        genome_count = 0
        for level, name in results[id].items():
            for place, entries in name.items():
                if len(entries) > 1:
                    genome_count += len(entries)
                    entry = str(len(entries)) + '\t' + place + ' (' + level + ')'
                    if entry not in reported:
                        print('\t' + entry)
                    if report != None and entry not in reported:
                        outfile2.write(entry + '\n')
                    reported.append(entry)

                    for genome in entries:
                        if genome not in [exported + to_keep]:
                            if id == 'strain':
                                outfile_names.write(genome + '\n')
                            else:
                                outfile_accno.write(genome + '\n')
                            exported.append(genome)

    # report selected samples listed in keep.txt
    not_found = []
    if len(to_keep) > 0:
        print('\t- ' + str(len(to_keep)) + ' genome(s) added from pre-selected list\n')
        outfile_names.write('\n# Pre-existing genomes listed in keep.txt\n')
        outfile_accno.write('\n# Pre-existing genomes listed in keep.txt\n')
        for genome in to_keep:
            if genome not in exported:
                outfile_names.write(genome + '\n')
                if genome in name_accno:
                    genome = name_accno[genome]
                    outfile_accno.write(genome + '\n')
                else:
                    not_found.append(genome)
                exported.append(genome)


        warning = 0
        for level, name in results['strain'].items():
            for place, entries in name.items():
                if len(entries) == 0:
                    if warning < 1:
                        if report != None:
                            outfile2.write('\n# Categories not found in metadata\n')
                        print('\n# Genomes matching the criteria below were not found')
                        warning += 1

                    entry = str(len(entries)) + '\t' + place + ' (' + level + ')'
                    print('\t' + entry)
                    if report != None:
                        outfile2.write(entry + '\n')

    if len(not_found) > 0:
        print('\n# ' + str(len(not_found)) + ' out of ' + str(len(to_keep)) + ' pre-existing genomes(s) were not found:')
        for name in not_found:
            print('\t- ' + name)

    print('\n' + str(genome_count) + ' genome(s) exported according to subsampling scheme\n')

