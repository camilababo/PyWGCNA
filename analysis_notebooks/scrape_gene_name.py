import pandas as pd
import requests
from bs4 import BeautifulSoup
import re

from bs4 import BeautifulSoup
import re

def scrapeGeneNameFromEnsemblPlant(data_expression_file, index_col=0, header=0, sep=',', batch_size=50):
    """
    Reads a gene expression csv file, retrieving the index column as the locus tags to search for in Ensembl Plant.
    :param sep: the separator used in the csv file
    :param header: the header row for the csv file
    :param index_col: the index column for the csv file
    :param data_expression_file: the gene expression csv file
    :param batch_size: the number of locus tags to search for in each batch
    :return: a txt file with the gene names and the locus tags
    """

    print(f'Function scrapeGeneNameFromEnsemblPlant() is running...')

    # Read the csv file
    df = pd.read_csv(data_expression_file, index_col=index_col, header=header, sep=sep)
    locus_tags = df.index.tolist()

    # Scrape the gene names from Ensembl Plant in batches
    gene_names = []
    gene_locus_tags = []  # Keep track of locus tags separately
    for i in range(0, len(locus_tags), batch_size):
        batch_locus_tags = locus_tags[i:i + batch_size]

        # Create a list to store the HTML responses
        responses = []
        for locus_tag in batch_locus_tags:
            url = f'https://plants.ensembl.org/Arabidopsis_thaliana/Gene/Summary?db=core;g={locus_tag};r=5:18586870-18588792'
            response = requests.get(url)
            responses.append(response)

        # Process the HTML responses outside the inner loop
        for response in responses:
            if response.status_code == 200:
                soup = BeautifulSoup(response.content, 'html.parser')
                gene_name_div = soup.find('h1', class_='summary-heading')
                if gene_name_div:
                    # Use regular expression to extract the gene name from the H1 element
                    gene_name_match = re.search(r'Gene:\s(.+)', gene_name_div.text)
                    if gene_name_match:
                        gene_name = gene_name_match.group(1)
                        gene_names.append(gene_name)
                    else:
                        print(f"Gene name not found for locus tag: {locus_tag}")
                        gene_names.append(locus_tag)  # Append locus tag when gene name is not found
                else:
                    print(f"Gene information div not found for locus tag: {locus_tag}")
                    gene_names.append(locus_tag)  # Append locus tag when gene information div is not found
                gene_locus_tags.append(locus_tag)  # Always append locus tag to the gene_locus_tags list

    # Save the gene names and locus tags to a txt file, comma separated and with headers
    output_path = data_expression_file.replace('.csv', '_geneList.txt')
    with open(output_path, 'w') as f:
        for name in gene_names:
            name_parts = name.split(' ', 1)
            if len(name_parts) == 2:
                gene_symbol = name_parts[0]
                locus_tag = name_parts[1]
                f.write(f'{gene_symbol} {locus_tag}\n')
            else:
                f.write(f'{name} {name}\n')


    print(f' Function scrapeGeneNameFromEnsemblPlant() has finished running.')
    print(f'Gene names and locus tags saved to {output_path}')

    return output_path


if __name__ == '__main__':
    input_path = '../data/E-MEXP-1304/2000_top_dev_genes_filtered_df.csv'
    scrapeGeneNameFromEnsemblPlant(data_expression_file=input_path, batch_size=100)
