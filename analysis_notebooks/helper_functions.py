import numpy as np
import pandas as pd


def log2_transform_transpose(expression_dataset, output_file, add_constant=False, constant=1):
    """
    Performs a log2 transformation and transpose the expression dataset.
    :param expression_dataset: The expression dataset
    :param output_file: The output file
    :param add_constant: Whether to add a constant to the Log2 transformation
    :param constant: The constant to add
    :return: The log2 transformed and transposed expression dataset
    """

    df = pd.read_csv(expression_dataset, index_col=0, header=0, sep=';')

    if add_constant:
        df = df + constant

    log2_df = np.log2(df)
    transposed_df = log2_df.transpose()
    transposed_df.to_csv(output_file, index=True, header=True, sep=',')

    print("Log2 transformed and transposed CSV saved successfully!")

def transpose_csv(input_file, output_file):
    """
    Transpose a csv file.
    :param input_file: The csv file to transpose
    :param output_file: The output file
    :return: The transposed csv file
    """

    df = pd.read_csv(input_file, index_col=0, header=0, sep=',')
    transposed_df = df.transpose()
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

def process_dataset_to_keep_only_first_locus_tag(expression_dataset, output_file):
    """
    Reads the expression dataset and removes any locus tag after the first one.
    :param expression_dataset: A csv file with the expression dataset
    :param output_file: The output file
    :return: A csv file with the expression dataset with only the first locus tag
    """

    df = pd.read_csv(expression_dataset, index_col=0, header=0, sep=',')
    df_index = df.index.tolist()

    filtered_index = [locus_tag.split(';')[0] for locus_tag in df_index]

    df.index = filtered_index

    df.to_csv(output_file, index=True, header=True, sep=',')

    print("Expression dataset index saved successfully!")

if __name__ == '__main__':
    input_file_path = "../data/E-MEXP-728/eMEXP728_processed_locus_tag.csv"
    output_file_path = "../data/E-MEXP-728/eMEXP728_non-filtered.csv"
    # transpose_csv(input_file_path, output_file_path)
    # rite_expression_dataset_index_to_csv(output_file_path)
    # process_dataset_to_keep_only_first_locus_tag(input_file_path, output_file_path)
    log2_transform_transpose(input_file_path, output_file_path, add_constant=True, constant=1)