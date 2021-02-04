# -*- coding: utf-8 -*-

import pycountry_convert as pyCountry
import pycountry
import pandas as pd
import argparse
from uszipcode import SearchEngine


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Reformat metadata file by adding column with subcontinental regions based on the UN geo-scheme",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", required=True, help="Nextstrain metadata file")
    parser.add_argument("--geoscheme", required=True, help="XML file with geographic classifications")
    parser.add_argument("--output", required=True, help="Updated metadata file")
    args = parser.parse_args()

    metadata = args.metadata
    geoscheme = args.geoscheme
    output = args.output

    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_variants/nextstrain/run6_20210202_b117/ncov/'
    # metadata = path + 'pre-analyses/metadata_filtered.tsv'
    # geoscheme = path + 'config/geoscheme.tsv'
    # output = path + 'pre-analyses/metadata_geo.tsv'
    #
    focus = ['USA', 'United Kingdom', 'Maine', 'New Hampshire',
             'Massachusetts', 'Connecticut', 'Vermont', 'New York']

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

    # parse subcontinental regions in geoscheme
    scheme_list = open(geoscheme, "r").readlines()[1:]
    geoLevels = {}
    c = 0
    for line in scheme_list:
        if not line.startswith('\n'):
            id = line.split('\t')[2]
            type = line.split('\t')[0]
            if type == 'region':
                members = line.split('\t')[5].split(',') # elements inside the subarea
                for country in members:
                    iso = get_iso(country.strip())
                    geoLevels[iso] = id

            # parse subnational regions for countries in geoscheme
            if type == 'country':
                members = line.split('\t')[5].split(',') # elements inside the subarea
                for state in members:
                    if state.strip() not in geoLevels.keys():
                        geoLevels[state.strip()] = id

            # parse subareas for states in geoscheme
            if type == 'location':
                members = line.split('\t')[5].split(',')  # elements inside the subarea
                for zipcode in members:
                    if zipcode.strip() not in geoLevels.keys():
                        geoLevels[zipcode.strip()] = id


    # open metadata file as dataframe
    dfN = pd.read_csv(metadata, encoding='utf-8', sep='\t')
    try:
        dfN.insert(4, 'region', '')
    except:
        pass
    dfN['region'] = dfN['iso'].map(geoLevels) # add 'column' region in metadata
    dfN['us_region'] = ''

    # convert sets of locations into sub-locations
    print('\nApplying geo-schemes...')
    dfN.fillna('', inplace=True)
    for idx, row in dfN.iterrows():

        # flatten divison names as country names, for countries that are not a focus of study
        country = dfN.loc[idx, 'country']
        if country not in focus:
            dfN.loc[idx, 'division'] = country

        # assign US region
        if country not in ['USA']:
            if 'Europe' in dfN.loc[idx, 'region']:
                dfN.loc[idx, 'us_region'] = 'Europe'
            else:
                dfN.loc[idx, 'us_region'] = 'Global'
        if country == 'USA' and dfN.loc[idx, 'us_region'] == '':
            dfN.loc[idx, 'us_region'] = geoLevels[dfN.loc[idx, 'division']]

        # convert sets of states into subnational regions
        division = dfN.loc[idx, 'division']
        if division not in ['', 'unknown']:
            if division in geoLevels.keys():
                dfN.loc[idx, 'country'] = geoLevels[dfN.loc[idx, 'division']]

        # flatten location names as division names for divisions that are not a focus of study
        if division not in focus:
            dfN.loc[idx, 'location'] = division
        print('Processing metadata for... ' + row['strain'])


    dfN = dfN.drop_duplicates(subset=['strain'])
    dfN.to_csv(output, sep='\t', index=False)

print('\nMetadata file successfully reformatted applying geo-scheme!\n')
