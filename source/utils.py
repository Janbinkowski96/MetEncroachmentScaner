import os
import multiprocessing
from typing import Tuple

import pandas as pd

def rename_axis(df: pd.DataFrame, mapper: dict) -> pd.DataFrame:
    """Function to rename dataframe axis"""

    if 'index' in mapper.keys():
        df.index.set_names("".join(mapper.values()), inplace=True)
    
    else:
        df = df.rename(mapper, axis=1)
    
    return df

def interval_size(bidirectional: bool, search_range: int, base_cpg_loc: int) -> Tuple[int, int]:
    if bidirectional:
        start_loc, stop_loc = base_cpg_loc - search_range, base_cpg_loc + search_range
    
    else:
        start_loc, stop_loc = base_cpg_loc, base_cpg_loc + search_range
        
    return start_loc, stop_loc

def check_cores_number(process_number: int) -> None:
    available_cpu = multiprocessing.cpu_count()
    if process_number > available_cpu:
        raise Exception("Selected number of processes exceeds the number of available cores.")

def create_path(path: str) -> str:
    return os.path.join(path, "report.csv")