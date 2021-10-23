import matplotlib

from table_info import BiomTable,CSVTable
import sklearn.cluster
from sklearn.metrics import silhouette_score
from matplotlib import pyplot as plt
import numpy as np

def find_nd_lines(df, genus_name, MIN_GENUS_COUNT):
    print("Starting Samples", df.shape[0])
    df_sum = df.sum(axis=1)

    # Filter points close to origin, they're just noise.
    filtered_df = df[df_sum >= MIN_GENUS_COUNT]
    print("Samples with", genus_name, " >", MIN_GENUS_COUNT, ": ", filtered_df.shape[0])

    # Project all remaining points to simplex
    simplex_df = filtered_df.divide(filtered_df.sum(axis=1), axis=0)

    best_clustering = None
    best_silhouette = -1
    # It's unlikely our samples have more species of this genus than exist in the reference db
    # Let's see what the optimal number of clusters is
    for k in range(2, simplex_df.shape[1] + 1):
        clusterer = sklearn.cluster.KMeans(n_clusters = k)
        clusterer.fit(simplex_df)
        labels = clusterer.labels_
        sil_score = silhouette_score(simplex_df, labels, metric='euclidean')
        if sil_score > best_silhouette:
            best_silhouette = sil_score
            best_clustering = clusterer

    print("Num", genus_name, "clusters in samples: ", len(best_clustering.cluster_centers_))
    print("Cluster Centers")
    print(best_clustering.cluster_centers_)
    print("Rough Cluster Names")
    for cluster in best_clustering.cluster_centers_:
        max_val = np.max(cluster)
        max_index = np.argmax(cluster)
        name = df.columns[max_index] + "-" + str(int(max_val * 100)) + "%"
        print(name)

    for col1_index in range(len(filtered_df.columns)):
        for col2_index in range(len(filtered_df.columns)):
            if col1_index == col2_index:
                continue
            c1 = filtered_df.columns[col1_index]
            c2 = filtered_df.columns[col2_index]
            maxx = max(filtered_df[c1])
            maxy = max(filtered_df[c2])
            if maxx < MIN_GENUS_COUNT or maxy < MIN_GENUS_COUNT:
                continue

            plt.subplot(1,2,1)
            plt.scatter(filtered_df[c1], filtered_df[c2], c=best_clustering.labels_, cmap="Set1")

            for line_index in range(len(best_clustering.cluster_centers_)):
                line = best_clustering.cluster_centers_[line_index]
                x = maxx
                y = line[col2_index] * maxx / line[col1_index]
                if y > maxy:
                    x = line[col1_index] * maxy / line[col2_index]
                    y = maxy

                cmap = matplotlib.cm.get_cmap('Set1')
                plt.plot([0, x], [0, y], c=cmap(line_index / (len(best_clustering.cluster_centers_)-1)))

            plt.axis('equal')
            plt.title('Clustered Reads:' + genus_name)
            plt.xlabel(c1)
            plt.ylabel(c2)

            plt.subplot(1,2,2)
            centers_x = []
            centers_y = []
            plt.scatter(simplex_df[c1], simplex_df[c2], c=best_clustering.labels_, cmap="Set1")
            for cluster_index in range(len(best_clustering.cluster_centers_)):
                centers_x.append(best_clustering.cluster_centers_[cluster_index][col1_index])
                centers_y.append(best_clustering.cluster_centers_[cluster_index][col2_index])

            plt.scatter(
                centers_x, centers_y,
                s=80, c=range(len(centers_x)),
                cmap="Set1", marker='X',
                edgecolors='black')
            plt.xlabel(c1)
            plt.show()


woltka_table = BiomTable("none")
metadata_table = CSVTable("./woltka_metadata.tsv", delimiter="\t")

woltka_table = woltka_table.load_dataframe()
metadata_table = metadata_table.load_dataframe()

for genus in ["Bacteroides"]:
    genus_ids = metadata_table[metadata_table['genus'] == genus].sort_values("species")[['#genome', 'genus', 'species']]
    genus_ids = [g for g in genus_ids['#genome'] if g in woltka_table.columns]
    df = woltka_table[genus_ids]
    find_nd_lines(df, genus, 500)