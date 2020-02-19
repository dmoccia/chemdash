import pandas as pd
import numpy as np
import multiprocessing
from multiprocessing import Pool

num_cores = multiprocessing.cpu_count()
print ('Number of cores...' + str(num_cores))
num_partitions = num_cores

def parallelize_dataframe_func(df: pd.DataFrame, func, *args) -> object:
    """
    #http://www.racketracer.com/2016/07/06/pandas-in-parallel/
    Args:
        df: input dataframe which will be parallelize
        func: the function you would like to call on the dataframe
        num_partitions: partitions to split df

    Returns:

    """
    df_split = np.array_split(df, num_partitions)
    pool = Pool(num_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df