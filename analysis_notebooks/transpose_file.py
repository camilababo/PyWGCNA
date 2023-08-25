import pandas as pd

def transpose_csv(input_file, output_file):
    # Read the CSV file
    df = pd.read_csv(input_file, index_col=0, header=0, sep=',')

    # Transpose the data
    transposed_df = df.transpose()

    # Save the transposed data to a new CSV file
    transposed_df.to_csv(output_file, index=True, header=True, sep=',')

    print("Transposed CSV saved successfully!")

def write_expression_dataset_index_to_csv(expression_dataset):
    """
    Write the index of the expression dataset to a csv file.
    :param expression_dataset: The expression dataset
    :return: A csv file with the index of the expression dataset
    """

    df = pd.read_csv(expression_dataset, index_col=0, header=0, sep=',')
    df_index = df.index.tolist()

    with open('expression_dataset_index.csv', 'w') as f:
        for item in df_index:
            f.write(f'{item}\n')

    print("Expression dataset index saved successfully!")

def process_dataset_to_keep_only_firt_locus_tag(expression_dataset, output_file):
    """
    Reads the expression dataset and removes any locus tag after the first one.
    :param expression_dataset: A csv file with the expression dataset
    :return: A csv file with the expression dataset with only the first locus tag
    """

    df = pd.read_csv(expression_dataset, index_col=0, header=0, sep=',')
    df_index = df.index.tolist()

    filtered_index = [locus_tag.split(';')[0] for locus_tag in df_index]

    df.index = filtered_index

    df.to_csv(output_file, index=True, header=True, sep=',')

    print("Expression dataset index saved successfully!")

input_file_path = "../data/E-MEXP-728/eMEXP728_filtered.csv"
output_file_path = "../data/E-MEXP-728/eMEXP728_filtered_transposed.csv"
transpose_csv(input_file_path, output_file_path)
# rite_expression_dataset_index_to_csv(output_file_path)
# process_dataset_to_keep_only_firt_locus_tag(input_file_path, output_file_path)