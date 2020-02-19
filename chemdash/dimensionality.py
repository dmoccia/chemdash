import math
from typing import List
import lap
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
import time
import umap
from sklearn import manifold


def add_dummy_rows(df:pd.DataFrame, id_col:str, bit_columns:List) -> (np.ndarray, List):

    num_dummy_rows = int(math.pow(math.ceil(math.sqrt(df.shape[0])), 2) - df.shape[0])
    df_dummy = pd.DataFrame(0, index=np.arange((num_dummy_rows)), columns=bit_columns)
    df_dummy[id_col] = 'Placeholder'
    return pd.concat([df, df_dummy]).reset_index(drop=True)


def separate_metadata(df:pd.DataFrame, prefix_bit_cols:str):

    # just want to get the fp columns list...
    features = df.filter(regex=prefix_bit_cols).columns
    # get remaining columns
    metadata = (list(set(df.columns) - set(features)))

    # Separating out the features from metadata
    x = df.loc[:, features].values

    # Standardizing the features  (do not need this for  fingerprints...)
    #x = StandardScaler().fit_transform(x)

    return x, metadata



def generate_tsne(df:pd.DataFrame, x:np.ndarray, metadata:List):
    # do tSNE
    tsne = manifold.TSNE(n_components=2, init='pca', random_state=0)
    trans_data = tsne.fit_transform(x)

    trans_data = pd.DataFrame(data=trans_data
                              , columns=['tsne_x', 'tsne_y'])

    print('tsne is complete')

    return pd.concat([trans_data, df[metadata]], axis=1)

def generate_umap(df:pd.DataFrame, x:np.ndarray, metadata:List):
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(x)
    #print (embedding.shape)
    embedding = pd.DataFrame(data=embedding
                              , columns=['umap_x', 'umap_y'])
    print ('umap is complete')

    return pd.concat([embedding, df[metadata]], axis=1)



def generate_lap_grid(df, x_col, y_col):

    side = int(math.sqrt(df.shape[0]))
    data2d = np.asarray(list(map(list,zip(df[x_col],df[y_col]))))

    #normalized output...
    data2d -= data2d.min(axis=0)
    data2d /= data2d.max(axis=0)

    #generate grid for mapping...
    xv, yv = np.meshgrid(np.linspace(0, 1, side), np.linspace(0, 1, side))
    grid = np.dstack((xv, yv)).reshape(-1, 2)


    #generate cost matrix
    print('cdist')
    print(time.time())
    cost = cdist(grid, data2d, 'sqeuclidean')

    #scale this cost matrix
    cost *= 100000 / cost.max()

    #useful tuturial on this...https://github.com/kylemcdonald/CloudToGrid/blob/master/CloudToGrid.ipynb

    #lapify...help('lap')
    print('lapjv')
    print(time.time())
    min_cost, row_assigns, col_assigns = lap.lapjv(np.copy(cost))

    grid_jv = grid[col_assigns]
    df_coords = pd.DataFrame(data=grid_jv, columns=['lap_x', 'lap_y'])

    return df.join(df_coords)

