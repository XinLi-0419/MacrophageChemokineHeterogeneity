#https://stackoverflow.com/questions/33960051/unable-to-import-a-module-from-python-notebook-in-jupyter
#https://stackoverflow.com/questions/58550576/importing-module-in-jupyter-notebook
#https://medium.com/@nrk25693/how-to-add-your-conda-environment-to-your-jupyter-notebook-in-just-4-steps-abeab8b8d084
# Applicable for terminal IDE direct run or Jupyter Notebook
# Run under pyscenic Environment with "conda activate pyscenic" in terminal or choosing pyscenic in Jupyter Notebook
import os
import glob
import pickle
import pandas as pd
import numpy as np
import csv

from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

# Directory
MM_TFS_FNAME = os.path.join("/Volumes/mLu/Analysis/scRNA_seq/Datasets/SCENIC/Transcription_Factors", "mm_mgi_tfs.txt")
DATABASES_GLOB = os.path.join("/Volumes/mLu/Analysis/scRNA_seq/Datasets/SCENIC/Regulon_Databases", "mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join("/Volumes/mLu/Analysis/scRNA_seq/Datasets/SCENIC/Motif_Annotation", "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")
SC_EXP_FNAME = os.path.join("/Volumes/mLu/Analysis/scRNA_seq/Datasets/SCENIC", "mLu.combined.IMs.exprMat.txt")
REGULONS_FNAME = os.path.join("/Volumes/mLu/Analysis/scRNA_seq/Datasets/SCENIC", "mLu.combined.IMs.500.100.regulons.p")
MOTIFS_FNAME = os.path.join("/Volumes/mLu/Analysis/scRNA_seq/Datasets/SCENIC", "mLu.combined.IMs.500.100.motifs.csv")
AUC_FNAME=os.path.join("/Volumes/mLu/Analysis/scRNA_seq/Datasets/SCENIC", "mLu.combined.IMs.500.100.auc.tsv")

# Expression matrix with rownames=names(cells) and colnames=names(genes)
ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0).T
ex_matrix.shape

# Import TFS
tf_names = load_tf_names(MM_TFS_FNAME)

# Import database
# db_fnames = glob.glob(DATABASES_GLOB)
db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
dbs

# Two lines would complete all
adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True) #Time-consuming
modules = list(modules_from_adjacencies(adjacencies, ex_matrix))

# Calculate a list of enriched motifs and the corresponding target genes for all modules.
with ProgressBar():
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

# Create regulons from this table of enriched motifs.
regulons = df2regulons(df)

# Save the enriched motifs and the discovered regulons to disk.
df.to_csv(MOTIFS_FNAME)
with open(REGULONS_FNAME, "wb") as f:
    pickle.dump(regulons, f)
    
# Instead of saving as loom. file, save AUC as csv. file
regulons = [r.rename(r.name.replace('(+)',' ('+str(len(r))+'g)')) for r in regulons]
auc_mtx = aucell(ex_matrix, regulons)
auc_mtx.to_csv(AUC_FNAME)

REGULONS_FNAME = os.path.join("/Volumes/mLu/Analysis/scRNA_seq/Datasets/SCENIC", "mLu.combined.IMs.500.100.regulons.csv")

# Specify the file name and mode ('w' for writing) of the CSV file
with open(REGULONS_FNAME, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Write each item in the list as a row in the CSV file
    for item in regulons:
        writer.writerow([item])



    
