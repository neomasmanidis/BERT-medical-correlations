import csv
import os
import pandas as pd
from tqdm import tqdm

# List of pathway maps to process
maps = []
with open('Useful maps.csv', 'r') as file:
    csvreader = csv.reader(file)
    for row in csvreader:
        maps.append(''.join(row))

df_curated_data = pd.DataFrame(columns=['head id', 'head name', 'tail id', 'tail name', 'pathway'])

print("Processing pathways data..")
status = tqdm(maps)
for m in status:
    status.set_description("Map %s" % m)
    # Load map data
    entries_table = pd.read_csv('data dump/' + m + ' entries.csv')
    if os.path.exists('data dump/'+m+' relations.csv'):
        relations_table = pd.read_csv('data dump/'+m+' relations.csv')

        # Iterate, for every relation in the relations table:
        for index, relation in relations_table.iterrows():
            link = relation['link']
            rel_name = relation['name']
            rel_value = relation['value']

            # Find the entry info from the entry table
            entry1 = int(relation['entry1'])  # This is a "local" id, only for data from the same kgml file
            entry2 = int(relation['entry2'])

            entry1_id = entries_table.loc[entries_table['id'] == entry1]['name'].to_list()[0]  # List only has 1 element
            entry2_id = entries_table.loc[entries_table['id'] == entry2]['name'].to_list()[0]

            if entry1_id == 'undefined':  # Some entries don't link to a database id
                continue
            if entry2_id == 'undefined':
                continue

            # In case of multiple ids (usually mutations), keep first id
            ls = entry1_id.split(" ")
            entry1_id = ls[0]
            ls = entry2_id.split(" ")
            entry2_id = ls[0]

            entry1_name = entries_table.loc[entries_table['id'] == entry1]['gene_names'].to_list()[0]
            ls = entry1_name.split(",")
            entry1_name = ls[0]  # Keep only 1st Alias as the entry name
            ls = entry1_name.split(" ")
            entry1_name = ls[0]

            entry2_name = entries_table.loc[entries_table['id'] == entry2]['gene_names'].to_list()[0]
            ls = entry2_name.split(",")
            entry2_name = ls[0]
            ls = entry2_name.split(" ")
            entry2_name = ls[0]

            if relation['name'] == 'compound':
                # In compound relations, there are 2 relations (Direction only assumed):
                # Relation 1: entry1 -> compound
                # Relation 2: compound -> entry2

                # -- Get compound name and id from entries table --
                comp = int(relation['value'])
                comp_id = entries_table.loc[entries_table['id'] == comp]['name'].to_list()[0]
                if comp_id == 'undefined':
                    continue
                ls = comp_id.split(" ")
                comp_id = ls[0]

                comp_name = entries_table.loc[entries_table['id'] == comp]['gene_names'].to_list()[0]
                ls = comp_name.split(",")
                comp_name = ls[0]

                # Relation 1 = entry1 -> compound
                temp = pd.DataFrame({'head id': entry1_id, 'head name': entry1_name, 'tail id': comp_id,
                                     'tail name': comp_name, 'link type': link, 'relation name': rel_name,
                                     'relation value': rel_value, 'entry1': entry1, 'entry2': entry2, 'pathway': m},
                                    index=[0])
                df_curated_data = pd.concat([df_curated_data, temp], ignore_index=True)

                # Relation 2 = compound - > entry2
                temp = pd.DataFrame({'head id': comp_id, 'head name': comp_name, 'tail id': entry2_id,
                                     'tail name': entry2_name, 'link type': link, 'relation name': rel_name,
                                     'relation value': rel_value, 'entry1': entry1, 'entry2': entry2, 'pathway': m},
                                    index=[0])
                df_curated_data = pd.concat([df_curated_data, temp], ignore_index=True)
            else:
                # In any other case, Relation = entry1 -> entry 2
                temp = pd.DataFrame({'head id': entry1_id, 'head name': entry1_name, 'tail id': entry2_id,
                                     'tail name': entry2_name, 'link type': link, 'relation name': rel_name,
                                     'relation value': rel_value, 'entry1': entry1, 'entry2': entry2, 'pathway': m},
                                    index=[0])
                df_curated_data = pd.concat([df_curated_data, temp], ignore_index=True)

    # Lastly, add the reaction data. Since they do not refer to the entries table, they don't need to be processed
    if os.path.exists('data dump/' + m + ' reactions.csv'):
        reactions_table = pd.read_csv('data dump/' + m + ' reactions.csv')
        df_curated_data = pd.concat([df_curated_data, reactions_table], ignore_index=True)


df_curated_data = df_curated_data.drop_duplicates(subset=['head id', 'tail id', 'relation name', 'relation value'])
df_curated_data.to_csv("All_relations-Curated.csv")
