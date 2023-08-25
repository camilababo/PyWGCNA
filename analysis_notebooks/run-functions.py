from pathlib import Path

import PyWGCNA

def run_jupyter_functions(dataset_path):
    PyWGCNA_1304 = PyWGCNA.WGCNA(name='E-MEXP-1304',
                                 species='arabidopsis thaliana',
                                 geneExpPath=dataset_path,
                                 outputPath='',
                                 save=True)
    PyWGCNA_1304.geneExpr.to_df()
    PyWGCNA_1304.preprocess()
    PyWGCNA_1304.findModules()



if __name__ == "__main__":
    dataset_path = Path("../data/emexp1304/emexp-1304-filtered.csv")
    run_jupyter_functions(dataset_path)