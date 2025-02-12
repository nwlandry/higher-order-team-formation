import pandas as pd
import xgi

from src import *

size = 3

conferences = [
    ("AES", 2017),
    ("CMC", 2018),
    ("MCL", 2015),
    ("TDA", 2015),
]
dataset_folder = "higher-order-datasets"

interaction_networks = {}
for c in conferences:
    series, year = c
    interaction_networks[c] = xgi.read_hif(
        f"{dataset_folder}/interaction_network_{series}_{year}.json"
    )

collaboration_networks = {}
for c in conferences:
    series, year = c
    collaboration_networks[c] = xgi.read_hif(
        f"{dataset_folder}/collaboration_network_{series}_{year}.json"
    )

t_c = {}
t_nc = {}
a_c = {}
a_nc = {}
for c in conferences:
    series, year = c
    H1 = interaction_networks[(series, year)]
    H2 = collaboration_networks[(series, year)]
    t_c[c], t_nc[c], a_c[c], a_nc[c] = utilities.compute_subgroup_interaction_times(
        H1, H2, size, return_aggregates=True
    )
    A, idx = create_data_matrix(t_nc[c], t_c[c], size)
    print(idx)
    columns = ["sync1", "async1", "async2", "async3", "collaborated"]
    df = pd.DataFrame(data=A, columns=columns)
    df.to_csv(f"data/{c[0]}_{c[1]}_matrix.csv", index=False)
