# BERT-medical-correlations
Using BERT models to find correlations between human genes and compounds

How to use:

1. KEGG Pathways Dataset Collection <br>
Run 'Dataset Collection Handler.py' to collect and process the data from KEGG.
The final data will be stored in 'All_relations-Curated-full_names.csv'

2. KEGG Undirected Graph Dataset<br>
Run 'KEGG Undirected Graph Split Louvain.ipynb' notebook to split the main graph into 2 groups from which seperate datasets are created.

3. Feature extraction<br>
Run 'Save embeddings.ipynb' notebook to create an embedding dictionary based on a pretrained BERT model.
Run 'feature-extraction-model.ipynb' notebook to create and train a neural network with the afformentioned embeddings as inputs according to the datasets created in step 2.

3. Finetuning<br>
Run 'bert-finetuning.ipynb' notebook to create and train a model based on a pretrained BERT model with the datasets created in step 2.

