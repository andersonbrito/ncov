# -*- coding: utf-8 -*-

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Release date: 2020-03-24
# Last update: 2020-12-31


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

    # metadata = path + 'metadata_filtered.tsv'
    # geoscheme = path + "geoscheme.tsv"
    # output = path + 'metadata_geo.tsv'

    focus = ['USA', 'Canada', 'United Kingdom', 'Maine', 'New Hampshire',
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
            if type == 'region_exposure':
                members = line.split('\t')[5].split(',') # elements inside the subarea
                for country in members:
                    iso = get_iso(country.strip())
                    geoLevels[iso] = id

            # parse subnational regions for countries in geoscheme
            if type == 'country_exposure':
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
        dfN.insert(4, 'region_exposure', '')
    except:
        pass
    dfN['region_exposure'] = dfN['iso'].map(geoLevels) # add 'column' region in metadata

    notfound = []
    # convert sets of locations into sub-locations
    print('\nApplying geo-schemes...')
    dfN.fillna('', inplace=True)
    search = SearchEngine(simple_zipcode=True)
    for idx, row in dfN.iterrows():

        # flatten divison names as country names, for countries that are not a focus of study
        country = dfN.loc[idx, 'country_exposure']
        if country not in focus:
            dfN.loc[idx, 'division_exposure'] = country

        # convert sets of states into subnational regions
        division = dfN.loc[idx, 'division_exposure']
        if division not in ['', 'unknown']:
            if division in geoLevels.keys():
                dfN.loc[idx, 'country_exposure'] = geoLevels[dfN.loc[idx, 'division_exposure']]

        # convert sets of cities into sub-state regions
        location = dfN.loc[idx, 'location']
        # print(location)
        if location not in ['', 'unknown'] and division == 'Connecticut':
            try:
                res = search.by_city_and_state(location, "CT")
                area_zip = res[0].zipcode
                if area_zip in geoLevels.keys():
                    dfN.loc[idx, 'location'] = geoLevels[area_zip]
                else:
                    print(row['location'] + ' has a zip code (' + area_zip + ') not found in the geo-scheme.)')
                    notfound.append(location)
            except:
                notfound.append(location)
                dfN.loc[idx, 'location'] = ''

        # flatten location names as division names for divisions that are not a focus of study
        if division not in focus:
            dfN.loc[idx, 'location'] = division
        print('Processing metadata for... ' + row['strain'])

    # report errors
    if len(notfound) > 0:
        print('\nSome locations were not assigned to sub-locations, and were not exported. Check for typos.')
        for entry in notfound:
            print('- ' + entry)

    dfN.to_csv(output, sep='\t', index=False)
print('\nMetadata file successfully reformatted applying geo-scheme!\n')
