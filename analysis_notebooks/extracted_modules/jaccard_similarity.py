import fileinput
import os
from typing import List

import math
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches


class CompareModulesBetweenApproaches:
    """
    This class is used to compare the modules between the different approaches.
    
    :param module_files: A list of paths to the files containing the modules
    :type module_files: List
    """""

    def __init__(self, module_files: List):
        self.module_files = module_files
        self.geneModule = {}
        self.moduleColors = list(self.geneModule.keys())

        self.jaccard_similarity = None
        self.fraction = None

    def create_combined_module_file(self, remove_non_overlapping_genes: bool = False):
        """
        Creates a txt file containing a dictionary with the gene modules.
        Does so by reading the original module txt files that contain the genes and the name of the module they belong to.
        :param remove_non_overlapping_genes: If True, only genes that are present in all files will be included.
        :return: A txt file containing a dictionary with the geneby  modules.
        """

        all_genes = set()  # To store all gene names
        module_files_data = {}  # To store data from each file

        # First, read all genes from all files
        for file_path in self.module_files:
            network_name = os.path.splitext(os.path.basename(file_path))[0]  # Extract network/module name

            genes = []
            colors = []
            with fileinput.input(files=file_path) as f:
                for line in f:
                    try:
                        gene, color = line.strip().split()  # Assuming each line has "gene color"
                        genes.append(gene)
                        colors.append(color)
                    except ValueError:
                        print(f"Line {f.filelineno()} in file {f.filename()} is not in the expected format.")

            module_files_data[network_name] = {"genes": genes, "colors": colors}
            all_genes.update(genes)

        # If remove_non_overlapping_genes is True, find common genes
        if remove_non_overlapping_genes:
            common_genes = all_genes.copy()
            for data in module_files_data.values():
                common_genes.intersection_update(data["genes"])
        else:
            common_genes = all_genes

        # Now, create DataFrames for only common genes
        for network_name, data in module_files_data.items():
            common_genes_indices = [i for i, gene in enumerate(data["genes"]) if gene in common_genes]
            common_genes = [data["genes"][i] for i in common_genes_indices]
            common_colors = [data["colors"][i] for i in common_genes_indices]

            module_df = pd.DataFrame({
                "moduleColors": common_colors,
                "index": common_genes
            })

            self.geneModule[network_name] = module_df

        return self.geneModule

    def calculate_jaccard_similarity(self):
        """
        This function calculates the jaccard similarity between the modules.
        :return: A pandas dataframe with the jaccard similarity between the modules
        """

        num = 0
        names = []

        for module in self.geneModule.keys():
            num = num + len(self.geneModule[module].moduleColors.unique())
            tmp = [f'{module}:' + s for s in self.geneModule[module].moduleColors.unique().tolist()]
            names = names + tmp

        jaccard_similarity = pd.DataFrame(0.0, columns=names, index=names)

        for network1 in self.geneModule.keys():
            for network2 in self.geneModule.keys():
                if network1 != network2:
                    modules1 = self.geneModule[network1].moduleColors.unique().tolist()
                    modules2 = self.geneModule[network2].moduleColors.unique().tolist()
                    for module1 in modules1:
                        for module2 in modules2:
                            list1 = self.geneModule[network1]["index"][
                                self.geneModule[network1].moduleColors == module1].tolist()
                            list2 = self.geneModule[network2]["index"][
                                self.geneModule[network2].moduleColors == module2].tolist()
                            jaccard_similarity.loc[f"{network1}:{module1}", f"{network2}:{module2}"] = \
                                self.jaccard_similarity_score(list1, list2)
                else:
                    modules = self.geneModule[network1].moduleColors.unique().tolist()
                    for module in modules:
                        jaccard_similarity.loc[f"{network1}:{module}", f"{network1}:{module}"] = 1.0

        self.jaccard_similarity = jaccard_similarity

        return jaccard_similarity

    def calculate_fraction_common_genes(self):
        """
        This function calculates the fraction of common genes between the modules.
        :return: A pandas dataframe with the fraction of common genes between the modules
        """

        num = 0
        names = []
        for network in self.geneModule.keys():
            num = num + len(self.geneModule[network].moduleColors.unique())
            tmp = [f"{network}:" + s for s in self.geneModule[network].moduleColors.unique().tolist()]
            names = names + tmp
        fraction = pd.DataFrame(0, columns=names, index=names)

        for network1 in self.geneModule.keys():
            for network2 in self.geneModule.keys():
                if network1 != network2:
                    modules1 = self.geneModule[network1].moduleColors.unique().tolist()
                    modules2 = self.geneModule[network2].moduleColors.unique().tolist()
                    for module1 in modules1:
                        for module2 in modules2:
                            list1 = self.geneModule[network1].index[
                                self.geneModule[network1].moduleColors == module1].tolist()
                            list2 = self.geneModule[network2].index[
                                self.geneModule[network2].moduleColors == module2].tolist()
                            num = np.intersect1d(list1, list2)
                            fraction.loc[f"{network1}:{module1}", f"{network2}:{module2}"] = len(num) / len(list2) * 100
                else:
                    modules = self.geneModule[network1].moduleColors.unique().tolist()
                    for module in modules:
                        fraction.loc[f"{network1}:{module}", f"{network1}:{module}"] = 1.0

        self.fraction = fraction

        return fraction


    def plot_jaccard_similarity(self,
                                color=None,
                                cutoff=0.1,
                                fig_size=None,
                                save=True,
                                plot_show=True,
                                plot_format="png",
                                file_name="jaccard_similarity",
                                network_colors=None
                                ):

        """
        This function plots the jaccard similarity between the modules.
        :param color: Possibility to color for each node of the network separately.
        :type color: Dict
        :param cutoff: The cutoff for the jaccard similarity
        :type cutoff: float
        :param fig_size: The size of the figure
        :type fig_size: tuple
        :param save: If True, the plot is saved
        :type save: bool
        :param plot_show: If True, the plot is shown
        :type plot_show: bool
        :param plot_format: The format of the plot
        :type plot_format: str
        :param file_name: The name of the file
        :type file_name: str
        :param network_colors: The colors of the networks
        :type network_colors: dict
        :return: A plot of the jaccard similarity between the modules
        """
        df = self.jaccard_similarity
        np.fill_diagonal(df.values, 0)
        df = pd.DataFrame(df.stack())
        df.reset_index(inplace=True)
        df = df[df[0] >= cutoff]
        df.columns = ['source', 'dest', 'weight']
        if df.shape[0] == 0:
            print(f"None of the connections pass the cutoff")
            return None

        G = nx.from_pandas_edgelist(df, 'source', 'dest', 'weight')
        node_labels = {}
        nodes = list(G.nodes())
        for i in range(len(nodes)):
            node_labels[nodes[i]] = nodes[i].split(":")[1]
        edges = G.edges()
        weights = [G[u][v]['weight'] * 10 for u, v in edges]
        edge_labels = {}
        for u, v in edges:
            edge_labels[u, v] = str(round(G[u][v]['weight'], 2))

        color_map = []
        if color is None:
            color_map = None
        else:
            for node in G:
                color_map.append(color[node.split(":")[0]])

        if fig_size is None:
            fig_size = (len(G.nodes()) / 2, len(G.nodes()) / 2)
        fig, ax = plt.subplots(figsize=fig_size, facecolor='white')
        pos = nx.spring_layout(G, k=1 / math.sqrt(len(G.nodes()) / 2))

        # Create a list of unique colors, one for each network
        networks = {node.split(":")[0] for node in G.nodes()}
        colors = list(mcolors.CSS4_COLORS.keys())  # Get a list of CSS4 color names
        color_map = {network: colors[i % len(colors)] for i, network in enumerate(networks)}

        node_colors = []  # New list

        for node in G:
            network_name = node.split(":")[0]
            if network_colors is not None and network_name in network_colors:
                color = network_colors[network_name]
            else:
                color = "default_color"  # Assign a default color if network_colors is not provided
            node_colors.append(color)

        nx.draw_networkx(G,
                         pos=pos,
                         node_color=node_colors,  # Use node_colors here
                         width=weights,
                         labels=node_labels,
                         font_size=8,
                         node_size=500,
                         with_labels=True,
                         ax=ax)

        legend_patches = [mpatches.Patch(color=color, label=network) for network, color in color_map.items()]

        nx.draw_networkx_edge_labels(G,
                                     pos,
                                     edge_labels=edge_labels,
                                     font_size=7)

        if color is not None:
            for label in color:
                color = {str(key): value for key, value in network_colors.items()}
                ax.plot([0], [0], color=color[label], label=label)

        plt.legend(handles=legend_patches)
        plt.tight_layout()
        if save:
            plt.savefig(f"{file_name}.{plot_format}")
        if plot_show:
            plt.show()
        else:
            plt.close()



    @staticmethod

    def jaccard_similarity_score(list1, list2):
        """
        This function calculates the jaccard similarity between two lists.
        :param list1: The first list of genes
        :type list1: list
        :param list2: The second list of genes
        :type list2: list
        :return: The jaccard similarity between the two lists
        """
        intersection = len(list(set(list1).intersection(list2)))
        union = (len(list1) + len(list2)) - intersection
        return float(intersection) / union


