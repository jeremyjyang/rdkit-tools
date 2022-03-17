#!/usr/bin/env python3
###
"""
Scipy/Scikit-learn hierarchical, agglomerative clustering (e.g. Ward) designed for
chemical fingerprints.

https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html
https://scikit-learn.org/stable/auto_examples/cluster/plot_agglomerative_dendrogram.html
https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html
"""
import sys,os,time,argparse,logging
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot
import scipy
from scipy.cluster.hierarchy import dendrogram,linkage

import sklearn
from sklearn.cluster import AgglomerativeClustering

#############################################################################
def PlotDendrogram(model, orientation, **kwargs):
  """
scipy.cluster.hierarchy.dendrogram(Z, p=30, truncate_mode=None, color_threshold=None, get_leaves=True, orientation='top', labels=None, count_sort=False, distance_sort=False, show_leaf_counts=True, no_plot=False, no_labels=False, leaf_font_size=None, leaf_rotation=None, leaf_label_func=None, show_contracted=False, link_color_func=None, ax=None, above_threshold_color='C0')
"""
  # Create linkage matrix (Z) and then plot the dendrogram.
  # Create the counts of samples under each node.
  counts = np.zeros(model.children_.shape[0])
  n_samples = len(model.labels_)
  for i, merge in enumerate(model.children_):
    current_count = 0
    for child_idx in merge:
      if child_idx < n_samples:
        current_count += 1  # leaf node
      else:
        current_count += counts[child_idx - n_samples]
    counts[i] = current_count

  # Stack 1-D arrays as columns into a 2-D array.
  Z = np.column_stack([model.children_, model.distances_, counts]).astype(float)

  # Generate the corresponding dendrogram
  ddg = dendrogram(Z, orientation=orientation, labels=model.labels_, **kwargs)
  return ddg

#############################################################################
def ClusterFeaturevectors(X, affinity, linkage):
  logging.info(f"X: {X.shape[0]}x{X.shape[1]}")
  model = AgglomerativeClustering(
	linkage=linkage, affinity=affinity,
	compute_distances=True,
	distance_threshold=0, n_clusters=None) # distance_threshold=0 ensures full tree.
  model = model.fit(X)
  model.labels_ = X.index #IDs
  return model

#############################################################################
def DescribeModel(model):
  logging.info(f"clusters: {model.n_clusters_}")
  logging.info(f"leaves: {model.n_leaves_}")
  labels = [str(label) for label in model.labels_]
  logging.debug(f"""labels ({len(labels)}): "{'","'.join(labels[:10])}"...""")
  logging.debug(f"feature_names_in ({len(model.feature_names_in_)}): {model.feature_names_in_}") #ndarray of shape (n_features_in_,)
  for i,row in enumerate(model.children_): #array-like of shape (n_samples-1,2)
    logging.debug(f"{i+1:5d}: merge {row[0]:5d}, {row[1]:5d} (d:{model.distances_[i]:6.3f})") #distances_: array-like of shape (n_nodes-1,)

#############################################################################
def WriteDistanceFile(model, ofile):
  fout = open(ofile, "w+")
  for i,nodePair in enumerate(model.children_):
    nodeA,nodeB = nodePair
    fout.write(f"{nodeA}\t{model.labels_[nodeA] if nodeA<len(model.labels_) else ''}\t{nodeB}\t{model.labels_[nodeB] if nodeB<len(model.labels_) else ''}\t{model.distances_[i]:.3f}\n")
  fout.close()
  logging.info(f"Output distance matrix written: {ofile}")

#############################################################################
def DescribeDendrogram(ddg):
  logging.debug(f"Dendrogram.color_list: {ddg['color_list']}")
  logging.debug(f"Dendrogram.ivl: {ddg['ivl'][:10]}...")

#############################################################################
if __name__=="__main__":
  parser = argparse.ArgumentParser(description="Hierarchical, agglomerative clustering by Scikit-learn", epilog="")
  AFFINITIES = ["euclidean", "l1", "l2", "manhattan", "cosine", "precomputed"]
  LINKAGES = ["ward", "complete", "average", "single"]
  ORIENTATIONS = ["left", "top", "right", "bottom"]
  OPS = ["cluster", "demo", ]
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("--i", dest="ifile", help="input file, TSV")
  parser.add_argument("--o", dest="ofile", help="output file, TSV")
  parser.add_argument("--o_dist", dest="ofile_dist", help="output distance file, TSV")
  parser.add_argument("--o_vis", dest="ofile_vis", default="/tmp/cluster_dendrogram.html", help="output file, PNG or HTML")
  parser.add_argument("--scratchdir", default="/tmp")
  parser.add_argument("--idelim", default="\t", help="delim for input TSV")
  parser.add_argument("--odelim", default="\t", help="delim for output TSV")
  parser.add_argument("--affinity", choices=AFFINITIES, default="euclidean", help="")
  parser.add_argument("--linkage", choices=LINKAGES, default="ward", help="")
  parser.add_argument("--truncate_level", type=int, default=4, help="Level from root of hierarchy for clusters and dendrogram.")
  parser.add_argument("--iheader", action="store_true", help="input TSV has header")
  parser.add_argument("--oheader", action="store_true", help="output TSV has header")
  parser.add_argument("--dendrogram_orientation", choices=ORIENTATIONS, default="left")
  parser.add_argument("--display", action="store_true", help="display dendrogram")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.debug(f"Scikit-learn version: {sklearn.__version__}")
  logging.debug(f"Scipy version: {scipy.__version__}")
  logging.debug(f"Numpy version: {np.__version__}")
  logging.debug(f"Matplotlib version: {mpl.__version__}")
  mpl.rcParams['figure.figsize'] = [10, 8]

  t0=time.time()

  if args.op=="demo":
    from sklearn.datasets import load_iris
    iris = load_iris()
    X = pd.DataFrame(iris.data, columns=iris.feature_names)
    model = ClusterFeaturevectors(X, args.affinity, args.linkage)
    DescribeModel(model)
    if args.display:
      pyplot.title(f"Hierarchical Clustering Dendrogram ({args.linkage}/{args.affinity})")
      ddg = PlotDendrogram(model, orientation=args.dendrogram_orientation, truncate_mode="level", p=args.truncate_level)
      DescribeDendrogram(ddg)
      pyplot.xlabel("Number of points in node (or index of point if no parenthesis).")
      pyplot.show()

  elif args.op=="cluster":
    if args.ifile is None:
      logging.error(f"--ifile required for operation: {args.op}")
    X = pd.read_csv(args.ifile, sep="\t")
    logging.info(f"Input file {args.ifile}: {X.shape[0]} x {X.shape[1]}")
    # Convert 1st col to index (must be unique IDs).
    X.columns = ["Name"]+list(X.columns[1:])
    X.set_index("Name", drop=True, verify_integrity=True, inplace=True)
    model = ClusterFeaturevectors(X, args.affinity, args.linkage)
    DescribeModel(model)
    if args.ofile_dist: WriteDistanceFile(model, args.ofile_dist)
    if args.display:
      pyplot.title(f"Hierarchical Clustering Dendrogram ({args.linkage}/{args.affinity})")
      ddg = PlotDendrogram(model, orientation=args.dendrogram_orientation, truncate_mode="level", p=args.truncate_level)
      DescribeDendrogram(ddg)
      if args.dendrogram_orientation in ("top", "bottom"):
        pyplot.yticks(ticks=[])
        pyplot.margins(y=.3)
      else:
        pyplot.xticks(ticks=[])
        pyplot.margins(x=.3)
      pyplot.xlabel("(Labels: N_in_cluster, or ID if singleton).")
      pyplot.show()

  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info(f"""Elapsed: {time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))}""")


