# -*- coding: utf-8 -*-

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Release date: 2020-03-24
# Last update: 2021-01-30


import pycountry_convert as pyCountry
import pycountry
from Bio import SeqIO
import pandas as pd
from epiweeks import Week
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genomes", required=True, help="FASTA file genomes to be used")
    parser.add_argument("--metadata1", required=True, help="Metadata file from NextStrain")
    parser.add_argument("--metadata2", required=False, help="Custom lab metadata file")
    parser.add_argument("--output1", required=True, help="Filtered metadata file")
    parser.add_argument("--output2", required=True, help="Reformatted, final FASTA file")
    args = parser.parse_args()

    genomes = args.genomes
    metadata1 = args.metadata1
    metadata2 = args.metadata2
    output1 = args.output1
    output2 = args.output2

    #     genomes = path + 'temp_sequences.fasta'
    #     metadata1 = path + 'metadata_nextstrain.tsv'
    #     metadata2 = path + 'GLab_SC2_sequencing_data.xlsx'
    #     output1 = path + 'metadata_filtered.tsv'
    #     output2 = path + 'sequences.fasta'

    # create a dict of existing sequences
    sequences = {}
    for fasta in SeqIO.parse(open(genomes), 'fasta'):  # as fasta:
        id, seq = fasta.description, fasta.seq
        if id not in sequences.keys():
            sequences[id] = str(seq)

    # get ISO alpha3 country codes
    isos = {}


    def get_iso(country):
        global isos
        if country not in isos.keys():
            try:
                isoCode = pyCountry.country_name_to_country_alpha3(country, cn_name_format="default")
                isos[country] = isoCode
            except:
                try:
                    isoCode = pycountry.countries.search_fuzzy(country)[0].alpha_3
                    isos[country] = isoCode
                except:
                    isos[country] = ''
        return isos[country]


    # create epiweek column
    def get_epiweeks(date):
        date = pd.to_datetime(date)
        epiweek = str(Week.fromdate(date, system="cdc"))  # get epiweeks
        epiweek = epiweek[:4] + '_' + 'EW' + epiweek[-2:]
        return epiweek


    # add state code

    us_state_abbrev = {
        'Alabama': 'AL',
        'Alaska': 'AK',
        'American Samoa': 'AS',
        'Arizona': 'AZ',
        'Arkansas': 'AR',
        'California': 'CA',
        'Colorado': 'CO',
        'Connecticut': 'CT',
        'Delaware': 'DE',
        'District of Columbia': 'DC',
        'Washington DC': 'DC',
        'Florida': 'FL',
        'Georgia': 'GA',
        'Guam': 'GU',
        'Hawaii': 'HI',
        'Idaho': 'ID',
        'Illinois': 'IL',
        'Indiana': 'IN',
        'Iowa': 'IA',
        'Kansas': 'KS',
        'Kentucky': 'KY',
        'Louisiana': 'LA',
        'Maine': 'ME',
        'Maryland': 'MD',
        'Massachusetts': 'MA',
        'Michigan': 'MI',
        'Minnesota': 'MN',
        'Mississippi': 'MS',
        'Missouri': 'MO',
        'Montana': 'MT',
        'Nebraska': 'NE',
        'Nevada': 'NV',
        'New Hampshire': 'NH',
        'New Jersey': 'NJ',
        'New Mexico': 'NM',
        'New York': 'NY',
        'North Carolina': 'NC',
        'North Dakota': 'ND',
        'Northern Mariana Islands': 'MP',
        'Ohio': 'OH',
        'Oklahoma': 'OK',
        'Oregon': 'OR',
        'Pennsylvania': 'PA',
        'Puerto Rico': 'PR',
        'Rhode Island': 'RI',
        'South Carolina': 'SC',
        'South Dakota': 'SD',
        'Tennessee': 'TN',
        'Texas': 'TX',
        'Utah': 'UT',
        'Vermont': 'VT',
        'Virgin Islands': 'VI',
        'Virginia': 'VA',
        'Washington': 'WA',
        'West Virginia': 'WV',
        'Wisconsin': 'WI',
        'Wyoming': 'WY'
    }

    # nextstrain metadata
    dfN = pd.read_csv(metadata1, encoding='utf-8', sep='\t', dtype='str')
    try:
        dfN.insert(4, 'iso', '')
    except:
        pass
    dfN['update'] = ''
    dfN['variants'] = ''

    dfN.fillna('', inplace=True)
    list_columns = dfN.columns.values  # list of column in the original metadata file
    # print(dfN)

    # Lab genomes metadata
    dfL = pd.read_excel(metadata2, index_col=None, header=0, sheet_name=0,
                        converters={'Sample-ID': str, 'Collection-date': str, 'Update': str})
    dfL.fillna('', inplace=True)

    dfL = dfL.rename(
        columns={'Sample-ID': 'id', 'Collection-date': 'date', 'Country': 'country', 'Division (state)': 'division',
                 'Location (county)': 'location', 'Country of exposure': 'country_exposure',
                 'State of exposure': 'division_exposure', 'Lineage': 'pango_lineage', 'Source': 'originating_lab',
                 'Update': 'update', 'Authors':'authors','Submitting lab':'submitting_lab'})
    if 'id' in dfL.columns.to_list():
        dfL = dfL[~dfL['id'].isin([''])]

    # add inexistent columns
    for col in list_columns:
        if col not in dfL.columns:
            dfL[col] = ''

    # output dataframe
    outputDF = pd.DataFrame(columns=list_columns)
    found = []
    lab_label = {}

    # process metadata from excel sheet
    for idx, row in dfL.iterrows():
        id = dfL.loc[idx, 'id']
        if id in sequences:
            dict_row = {}
            for col in list_columns:
                dict_row[col] = ''
                if col in row:
                    dict_row[col] = dfL.loc[idx, col]  # add values to dictionary

            if dict_row['location'] in ['', None]:
                dict_row['location'] = dfL.loc[idx, 'location']

            collection_date = ''
            if len(str(dict_row['date'])) > 1:
                collection_date = dict_row['date'].split(' ')[0].replace('.', '-').replace('/', '-')
                dict_row['date'] = collection_date

            # fix exposure
            columns_exposure = ['country_exposure', 'division_exposure']
            for level_exposure in columns_exposure:
                level = level_exposure.split('_')[0]
                dict_row[level_exposure] = dfL.loc[idx, level_exposure]
                if dict_row[level_exposure] in ['', None]:
                    if level_exposure == 'country_exposure':
                        dict_row[level_exposure] = dict_row[level]
                    else:
                        if dict_row['country_exposure'] != dfL.loc[idx, 'country']:
                            dict_row[level_exposure] = dict_row['country_exposure']
                        else:
                            dict_row[level_exposure] = dict_row[level]

            code = ''
            if dict_row['division'] in us_state_abbrev:
                code = us_state_abbrev[dict_row['division']] + '-'

            strain = dfL.loc[idx, 'country'] + '/' + code + dfL.loc[idx, 'id'] + '/' + collection_date.split('-')[
                0]  # set the strain name
            dict_row['strain'] = strain
            dict_row['iso'] = get_iso(dict_row['country'])
            dict_row['submitting_lab'] = dict_row['submitting_lab']
            dict_row['authors'] = dict_row['authors']

            # add lineage
            lineage = ''
            if dfL.loc[idx, 'pango_lineage'] != '':
                lineage = dfL.loc[idx, 'pango_lineage']
            dict_row['pango_lineage'] = lineage

            if dfL.loc[idx, 'update'] == 'variant':
                dict_row['variants'] = get_epiweeks(collection_date)
                dict_row['update'] = ''
            else:
                dict_row['update'] = 'Update' + str('0' * (2 - len(dfL.loc[idx, 'update']))) + dfL.loc[idx, 'update']
                dict_row['variants'] = ''
            found.append(strain)
            lab_label[id] = strain

            outputDF = outputDF.append(dict_row, ignore_index=True)

    # process metadata from TSV
    dfN = dfN[dfN['strain'].isin(sequences.keys())]
    for idx, row in dfN.iterrows():
        strain = dfN.loc[idx, 'strain']
        if strain in sequences:
            if strain in outputDF['strain'].to_list():
                continue
            dict_row = {}
            for col in list_columns:
                dict_row[col] = ''
                if col in row:
                    dict_row[col] = dfN.loc[idx, col]

            # fix exposure
            columns_exposure = ['country_exposure', 'division_exposure']
            for level_exposure in columns_exposure:
                level = level_exposure.split('_')[0]
                if dict_row[level_exposure] in ['', None]:
                    dict_row[level_exposure] = dict_row[level]

            dict_row['iso'] = get_iso(dict_row['country'])
            found.append(strain)

            outputDF = outputDF.append(dict_row, ignore_index=True)

    # write new metadata files
    outputDF = outputDF.drop(columns=['region'])
    outputDF.to_csv(output1, sep='\t', index=False)

    # write sequence file
    exported = []
    with open(output2, 'w') as outfile2:
        # export new metadata lines
        for id, sequence in sequences.items():
            if id in lab_label:
                if lab_label[id] not in exported:
                    entry = '>' + lab_label[id] + '\n' + sequence + '\n'
                    outfile2.write(entry)
                    print('* Exporting newly sequenced genome and metadata for ' + id)
                    exported.append(lab_label[id])
            else:
                if id not in exported:
                    entry = '>' + id + '\n' + sequence + '\n'
                    outfile2.write(entry)
                    exported.append(id)

print('\nMetadata file successfully reformatted and exported!\n')
