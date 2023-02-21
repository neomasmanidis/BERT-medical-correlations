import csv
from tqdm import tqdm
import requests
import pandas as pd
import xmltodict

print("Loading pathway maps list..")
r = requests.get('http://rest.kegg.jp/list/pathway')
# Returned as plain text, separated by new lines
s = r.text.split("\n")[:-1]  # (Last line is empty)
maps = []
for entry in s:
    maps.append(entry[8:13])  # keep only map number, organism will be human only (hsa)

not_hsa = []        # Will store maps that do not provide kgml file
only_entries = []   # Will store maps that provide an entries table, but don't provide a relation or reactions table
useful_maps = []    # Will store maps that provided a kgml file, entries and relations or reactions table

print("Loading data from each map..")
status = tqdm(maps)
for m in status:
    status.set_description("Map %s" % m)
    # Request kgml file from kegg api
    r = requests.get('http://rest.kegg.jp/get/hsa' + m + '/kgml', allow_redirects=True)
    if r.status_code != 404:
        data = xmltodict.parse(r.text)

        # ----- Handle Entities Data -----
        df_entries = pd.DataFrame(data['pathway']['entry'])
        df_entries = df_entries.rename(columns={"@id": "id", "@name": "name", "@type": "type", "@link": "link",
                                                "@reaction": "reaction"})
        if 'reaction' not in df_entries.keys():  # Not all entities have reaction data
            df_entries['reaction'] = ''
        # Extract data from nested "graphics" entry
        df_entries = pd.concat([df_entries.drop(['graphics'], axis=1),
                                df_entries['graphics'].apply(pd.Series)], axis=1)
        # Drop unnecessary data
        df_entries = df_entries.drop(['@fgcolor', '@bgcolor', '@type', '@x', '@y', '@width', '@height'], axis=1)
        # Organize
        df_entries = df_entries.rename({"@name": "gene_names"}, axis=1)
        df_entries = df_entries[['id', 'name', 'type', 'link', 'gene_names', 'reaction']]
        # Store to csv in "data dump" folder
        df_entries.to_csv('data dump/hsa' + m + ' entries.csv')

        # ----- Handle Relations data -----
        if 'relation' in data['pathway'].keys():
            # Check if single entry or list of dictionaries
            if isinstance(data['pathway']['relation'], list):
                df_relations = pd.DataFrame.from_dict(data['pathway']['relation'])
            else:
                df_relations = pd.DataFrame.from_dict(data['pathway']['relation'], orient='index').T

            if 'subtype' in df_relations.keys():  # If there are subtype data (Relation type name/value)
                # Check for multiple subtypes
                df_supp = pd.DataFrame()  # DataFrame to hold extracted rows (supplementary rows)
                drop_list = []  # Which rows to drop later
                for i in range(0, len(df_relations.index)):
                    # If there are multiple subtypes, there will be a list
                    if isinstance(df_relations['subtype'][i], list):
                        for j in range(0, len(df_relations['subtype'][i])):
                            tempdf = pd.DataFrame({"@entry1": [df_relations['@entry1'][i]],
                                                   "@entry2": [df_relations['@entry2'][i]],
                                                   "@type": [df_relations['@type'][i]],
                                                   "@name": [df_relations['subtype'][i][j]["@name"]],
                                                   "@value": [df_relations['subtype'][i][j]["@value"]]})
                            df_supp = pd.concat([df_supp, tempdf],
                                                ignore_index=True)
                        drop_list.append(i)
                # If there were multiple subtypes in the same row, drop them. they are stored in df_supp and added later
                df_relations.drop(drop_list, axis=0, inplace=True)
                # Extract subtype data to separate columns
                df_relations = pd.concat([df_relations.drop(['subtype'], axis=1),
                                          df_relations['subtype'].apply(pd.Series)], axis=1)
                # Append supplementary rows
                df_relations = pd.concat([df_relations, df_supp], ignore_index=True)

            else:  # If no subtype data, then fill in blanks
                df_relations['@name'] = ''
                df_relations['@value'] = ''

            # Organize data
            df_relations = df_relations.rename(columns={"@entry1": "entry1", "@entry2": "entry2", "@type": "link",
                                               "@name": "name", "@value": "value"})
            df_relations = df_relations[['entry1', 'entry2', 'link', 'value', 'name']]
            # Add relations pathway source
            df_relations['pathway'] = 'hsa' + m
            # Store to csv in "data dump" folder
            df_relations.to_csv('data dump/hsa' + m + ' relations.csv')

        # ----- Handle Reactions data -----
        if 'reaction' in data['pathway'].keys():
            # Check if single entry or list of dictionaries
            if isinstance(data['pathway']['reaction'], list):
                df_reactions = pd.DataFrame.from_dict(data['pathway']['reaction'])
            else:
                df_reactions = pd.DataFrame.from_dict(data['pathway']['reaction'], orient='index').T

            df_reactions_relations = pd.DataFrame()
            for i, row in df_reactions.iterrows():
                if isinstance(df_reactions['substrate'][i], list):
                    substrate_list = df_reactions['substrate'][i]
                else:
                    substrate_list = [df_reactions['substrate'][i]]

                if isinstance(df_reactions['product'][i], list):
                    product_list = df_reactions['product'][i]
                else:
                    product_list = [df_reactions['product'][i]]

                for j in range(0, len(substrate_list)):
                    for k in range(0, len(product_list)):
                        # If multiple similar ids separated by space, keep the 1st
                        ls = substrate_list[j]['@name'].split(" ")
                        head_id = ls[0]
                        ls = product_list[k]['@name'].split(" ")
                        tail_id = ls[0]
                        temp = pd.DataFrame({'head id': head_id,
                                             'head name': substrate_list[j]['@name'],
                                             'tail id': tail_id,
                                             'tail name': product_list[k]['@name'], 'link type': 'reaction',
                                             'relation name': df_reactions['@name'][i],
                                             'relation value': df_reactions['@type'][i],
                                             'entry1': substrate_list[j]['@id'],
                                             'entry2': product_list[k]['@id'],
                                             'pathway': 'hsa'+m},
                                            index=[0])

                        df_reactions_relations = pd.concat([df_reactions_relations, temp], ignore_index=True)

            df_reactions_relations.to_csv('data dump/hsa' + m + ' reactions.csv')

        if 'relation' not in data['pathway'].keys() and 'reaction' not in data['pathway'].keys():
            only_entries.append(m)
        else:
            useful_maps.append('hsa' + m)

    else:
        # print("No kgml file found")
        not_hsa.append(m)


with open('Useful maps.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    for m in useful_maps:
        writer.writerow([m])

with open('Only entries maps.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    for m in only_entries:
        writer.writerow([m])

with open('Not human maps.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    for m in not_hsa:
        writer.writerow([m])
