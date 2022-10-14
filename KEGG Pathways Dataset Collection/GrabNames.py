import requests
from tqdm import tqdm
import pandas as pd


def extract_gene_name(text):
    ls = text.split('\n')
    symbol = ls[1]
    name = ls[2]
    # -- Gene symbol --
    symbol = symbol.strip()
    symbol = symbol[12:]
    ls = symbol.split(",")
    symbol = ls[0]
    # -- Gene name --
    name = name.strip()
    name = name[21:]

    return symbol+', '+name


def extract_compound_name(text):
    ls = text.split('\n')
    name = ls[1]
    name = name.strip()
    name = name[12:]
    name = name.strip(";")
    return name


def extract_glycan(text):
    ls = text.split('\n')
    name = ls[3]
    name = name.strip()
    # Some glycans are the same as compounds, therefore store them as compounds
    if name[0:21] == "REMARK      Same as: ":
        return "cpd:"+name[21:27]
    else:
        name = ls[1]
        name = name[12:]
        name = name.strip(";")
        return name


df_final_data = pd.read_csv("All_relations-Curated.csv")
df_final_data['head full name'] = ''
df_final_data['tail full name'] = ''
found_names_dict = {}

print("Grabbing full names from Kegg..")
status = tqdm(df_final_data.iterrows(), total=df_final_data.shape[0])
for index, entry in status:

    # search dictionary first for head name
    if entry['head id'] in found_names_dict.keys():
        df_final_data.at[index, 'head full name'] = found_names_dict[entry['head id']]
    else:
        # Request info from Kegg api
        r1 = requests.get('http://rest.kegg.jp/get/' + entry['head id'])

        if entry['head id'][0:3] == 'hsa':  # If human gene
            head_name = extract_gene_name(r1.text)
            df_final_data.at[index, 'head full name'] = head_name
        elif entry['head id'][0:2] == 'gl':  # If glycan
            head_name = extract_glycan(r1.text)
            if head_name[0:3] == 'cpd':  # If this glycan is also a compound
                r1 = requests.get('http://rest.kegg.jp/get/' + head_name)
                # change head id to cpd
                entry['head id'] = head_name
                # find cpd name
                head_name = extract_compound_name(r1.text)
            df_final_data.at[index, 'head full name'] = head_name
        else:  # In every other case try as compound
            try:
                head_name = extract_compound_name(r1.text)
                df_final_data.at[index, 'head full name'] = head_name
            except:  # If unable to extract info, fill blanks
                df_final_data.at[index, 'head full name'] = ""

        found_names_dict[entry['head id']] = df_final_data.at[index, 'head full name']

    status.set_description("%s -> %s " % (entry['head id'], df_final_data.at[index, 'head full name']))

    # search dictionary first for tail name
    if entry['tail id'] in found_names_dict.keys():
        df_final_data.at[index, 'tail full name'] = found_names_dict[entry['tail id']]
    else:
        # Request info from Kegg api
        r2 = requests.get('http://rest.kegg.jp/get/' + entry['tail id'])

        if entry['tail id'][0:3] == 'hsa':
            tail_name = extract_gene_name(r2.text)
            df_final_data.at[index, 'tail full name'] = tail_name
        elif entry['tail id'][0:2] == 'gl':
            head_name = extract_glycan(r2.text)
            if head_name[0:3] == 'cpd':
                r2 = requests.get('http://rest.kegg.jp/get/' + head_name)
                # change head id to cpd
                entry['tail id'] = head_name
                # find cpd name
                head_name = extract_compound_name(r2.text)
            df_final_data.at[index, 'tail full name'] = head_name
        else:
            try:
                tail_name = extract_compound_name(r2.text)
                df_final_data.at[index, 'tail full name'] = tail_name
            except:
                df_final_data.at[index, 'tail full name'] = ""

        found_names_dict[entry['tail id']] = df_final_data.at[index, 'tail full name']

    status.set_description("%s -> %s " % (entry['tail id'], df_final_data.at[index, 'tail full name']))

df_final_data.to_csv("All_relations-Curated-full_names.csv")
