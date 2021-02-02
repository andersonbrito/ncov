# -*- coding: utf-8 -*-

import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata1", required=True, help="Metadata file from nextstrain")
    parser.add_argument("--metadata2", required=True, help="Metadata file from GISAID")
    parser.add_argument("--output", required=True, help="Merged metadata file")
    args = parser.parse_args()

    metadata1 = args.metadata1
    metadata2 = args.metadata2
    output = args.output

    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_variants/nextstrain/run4_20210116_variants/pre-analyses/'
    # metadata1 = path + 'metadata_nextstrain.tsv'
    # metadata2 = path + 'gisaid_hcov-19_2021_01_18_07.tsv'
    # output = path + 'metadata_merged.tsv'

    separator1 = ''
    if str(metadata1).split('.')[-1] == 'tsv':
        separator1 = '\t'
    elif str(metadata1).split('.')[-1] == 'csv':
        separator1 = ','

    separator2 = ''
    if str(metadata2).split('.')[-1] == 'tsv':
        separator2 = '\t'
    elif str(metadata2).split('.')[-1] == 'csv':
        separator2 = ','


    # nextstrain metadata
    dfN = pd.read_csv(metadata1, encoding='utf-8', sep=separator1, dtype='str')
    dfN.fillna('', inplace=True)

    # GISAID metadata
    dfG = pd.read_csv(metadata2, encoding='utf-8', sep=separator2, dtype='str')
    dfG = dfG.rename(columns={'Virus name': 'strain', 'Accession ID': 'gisaid_epi_isl', 'Collection date': 'date',
                              'Host': 'host', 'Gender': 'sex', 'Patient age': 'age', 'Lineage': 'pangolin_lineage'})

    list_columns = [col for col in dfG.columns.to_list() if col in dfN.columns.to_list() + ['Location']]  # list of columns in common
    dfG = dfG[list_columns]

    # fix columns
    dfG['strain'] = dfG['strain'].str.replace('hCoV-19/', '')
    dfG['region'] = dfG['Location'].str.split(' / ').str[0]
    dfG['country'] = dfG['Location'].str.split(' / ').str[1]
    try:
        dfG['division'] = dfG['Location'].str.split(' / ').str[2]
    except:
        dfG['division'] = ''

    try:
        dfG['location'] = dfG['Location'].str.split(' / ').str[3]
    except:
        dfG['location'] = ''

    dfG = dfG.drop(columns=['Location'])


    found = []
    for entry in dfG['strain'].to_list():
        if entry in dfN['strain'].to_list():
            found.append(entry)

    dfG = dfG[~dfG['strain'].isin(found)] # avoid duplicates

    if len(found) > 0:
        c = 1
        print('\n The following entries were already present in the original metadata file:\n')
        for i in found:
            print('\t' + str(c) + '. ' + i)
            c += 1

    # merge frames
    frames = [dfN, dfG]
    result = pd.concat(frames)
    result.fillna('', inplace=True)
    result.to_csv(output, sep='\t', index=False)

    print('\nTSV metadata files successfully merged.\n')
