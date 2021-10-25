import matplotlib

from table_info import BiomTable,CSVTable
import sklearn.cluster
from sklearn.metrics import silhouette_score
from matplotlib import pyplot as plt
import numpy as np


def plot_scatter(filtered_df, simplex_df, best_clustering, title, col1_index, col2_index):
    if col1_index == col2_index:
        print(title, "Plot Skipped (same columns)")
        return
    c1 = filtered_df.columns[col1_index]
    c2 = filtered_df.columns[col2_index]
    maxx = max(filtered_df[c1])
    maxy = max(filtered_df[c2])

    plt.subplot(1,2,1)
    if best_clustering is None:
        plt.scatter(filtered_df[c1], filtered_df[c2])
    else:
        plt.scatter(filtered_df[c1], filtered_df[c2], c=best_clustering.labels_, cmap="Set1")

    if best_clustering is not None:
        for line_index in range(len(best_clustering.cluster_centers_)):
            line = best_clustering.cluster_centers_[line_index]
            lx = line[col1_index]
            ly = line[col2_index]

            if lx == 0 and ly == 0:
                x = 0
                y = 0
            elif lx == 0:
                x = 0
                y = maxy
            elif ly == 0:
                x = maxx
                y = 0
            else:
                x = maxx
                y = ly * maxx / lx
                if y > maxy:
                    x = lx * maxy / ly
                    y = maxy

            cmap = matplotlib.cm.get_cmap('Set1')
            plt.plot([0, x], [0, y], c=cmap(line_index / (len(best_clustering.cluster_centers_)-1)))

    plt.axis('equal')
    plt.title('Clustered Reads:' + title)
    plt.xlabel(c1)
    plt.ylabel(c2)

    plt.subplot(1,2,2)
    centers_x = []
    centers_y = []
    if best_clustering is None:
        plt.scatter(simplex_df[c1], simplex_df[c2])
    else:
        plt.scatter(simplex_df[c1], simplex_df[c2], c=best_clustering.labels_, cmap="Set1")

    if best_clustering is not None:
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


def find_nd_lines(df, genus_name, MIN_GENUS_COUNT):
    print("Starting Samples", df.shape[0])
    df_sum = df.sum(axis=1)

    # Determine how many samples are relevant
    any_df = df[df_sum >= 0]
    print("Samples with", genus_name, " >", 0, ": ", any_df.shape[0])

    # Can potentially use a dynamic filter here, rather than a fixed threshold.
    # Filter points close to origin, they're just noise.
    filtered_df = df[df_sum >= MIN_GENUS_COUNT]
    print("Samples with", genus_name, " >", MIN_GENUS_COUNT, ": ", filtered_df.shape[0])

    # Project all remaining points to simplex
    simplex_df = filtered_df.divide(filtered_df.sum(axis=1), axis=0)

    best_clustering = None
    best_silhouette = -1
    # It's unlikely our samples have more species of this genus than exist in the reference db
    # Let's see what the optimal number of clusters is
    # At least 10 samples per line
    MIN_POINTS_PER_LINE = 10
    if simplex_df.shape[1] < 2:
        print("Only one reference genome, skip")
        return -1
    elif int(simplex_df.shape[0] / MIN_POINTS_PER_LINE) < 2:
        print("Not enough samples containing", genus_name, "to identify ratios, skip")
        return -1

    max_clusters = min(int(simplex_df.shape[0] / MIN_POINTS_PER_LINE), simplex_df.shape[1])

    MAX_KMEANS_ITER = 300
    for k in range(2, max_clusters + 1):
        clusterer = sklearn.cluster.KMeans(n_clusters = k, max_iter=MAX_KMEANS_ITER)
        clusterer.fit(simplex_df)
        if clusterer.n_iter_ == MAX_KMEANS_ITER:
            print("k=",k,": Convergence Failed")
            continue
        labels = clusterer.labels_
        sil_score = silhouette_score(simplex_df, labels, metric='euclidean')
        if sil_score > best_silhouette:
            best_silhouette = sil_score
            best_clustering = clusterer

    if best_clustering == None:
        print("No clustering could be identified, may indicate only 1 axis.  TODO: Implement fallback")
    else:
        print("Num", genus_name, "clusters in samples: ", len(best_clustering.cluster_centers_))
        print("Cluster Centers")
        print(best_clustering.cluster_centers_)
        print("Rough Cluster Names")
        for cluster in best_clustering.cluster_centers_:
            max_val = np.max(cluster)
            max_index = np.argmax(cluster)
            name = df.columns[max_index] + "-" + str(int(max_val * 100)) + "%"
            print(name)

    for col1_index in range(len(filtered_df.columns)-1):
        pass
        # plot_scatter(
        #     filtered_df, simplex_df,
        #     best_clustering,
        #     genus_name + " Score:" + str(best_silhouette),
        #     col1_index, col1_index+1
        # )

    return best_silhouette


if __name__ == '__main__':
    woltka_table = BiomTable("none")
    metadata_table = CSVTable("./woltka_metadata.tsv", delimiter="\t")

    woltka_table = woltka_table.load_dataframe()
    metadata_table = metadata_table.load_dataframe()

    sil_map = {}
    all_genera = metadata_table['genus'].unique()
    all_genera = [str(g) for g in all_genera]
    all_genera.sort()
    for genus in all_genera:
        genus_ids = metadata_table[metadata_table['genus'] == genus].sort_values("species")[['#genome', 'genus', 'species']]
        genus_ids = [g for g in genus_ids['#genome'] if g in woltka_table.columns]
        if len(genus_ids) < 2:
            continue
        df = woltka_table[genus_ids]
        sil_map[genus] = find_nd_lines(df, genus, 500)

    plt.hist(list(sil_map.values()))
    plt.show()

    failures = 0
    good_ones = 0
    for g in sil_map:
        if sil_map[g] > 0.85:
            good_ones += 1
        if sil_map[g] == -1:
            failures += 1

    print("Failures: ", failures)
    print("Good Ones: ", good_ones)
    print("Total Ones: ", all_genera)
