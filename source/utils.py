import os
import multiprocessing
from contextlib import suppress

import pandas as pd


def rename_axis(df: pd.DataFrame, mapper: dict) -> pd.DataFrame:
    """Function to rename dataframe axis"""

    if 'index' in mapper.keys():
        df.index.set_names("".join(mapper.values()), inplace=True)
    
    else:
        df = df.rename(mapper, axis=1)
    
    return df


def check_cores_number(process_number: int) -> None:
    available_cpu = multiprocessing.cpu_count()
    if process_number > available_cpu:
        raise Exception("Selected number of processes exceeds the number of available cores.")


def create_path(path: str, name="report.csv") -> str:
    return os.path.join(path, name)


def create_dir(path:str, filename: str) -> None:
    with suppress(FileExistsError):
        os.mkdir(create_path(path, filename))
