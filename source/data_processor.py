import pandas as pd
from tqdm import tqdm

from .utils import create_path


class DataProcessor:
    def __init__(self, manifest_type):
        self.manifest_type = manifest_type
        self.island_annotations = None
        self.raw_mynorm = None
        self.cg_per_chr = None
        self.mynorm_std = None
        self.manifest = None
        self.mynorm = None

    def set_manifest(self) -> None:
        if self.manifest_type == "E":
            self.manifest = pd.read_csv("resources/EPIC/EPIC.csv", index_col=0, low_memory=False)
        else:
            self.manifest = pd.read_csv("resources/450K/450K.csv", index_col=0, low_memory=False)

    @staticmethod
    def first_load(mynorm_path: str) -> dict:
        df = pd.read_csv(mynorm_path, nrows=5)

        dtypes = df.dtypes
        col_names = dtypes.index
        types = [i.name for i in dtypes.values]
        column_types = dict(zip(col_names, types))

        return column_types

    def load_mynorm(self, mynorm_path: str, column_types: dict) -> None:
        self.raw_mynorm = pd.read_csv(mynorm_path, index_col=0, encoding="latin1", dtype=column_types)
        self.mynorm = self.raw_mynorm.mean(axis=1).to_frame(name="beta-values")
        self.mynorm_std = self.raw_mynorm.std(axis=1).to_frame(name="beta-values std")

    def select_cpg(self) -> None:
        mynorm_cpg = set(self.mynorm.index)
        overlapped = set.intersection(mynorm_cpg, set(self.manifest.index))
        self.manifest = self.manifest.loc[overlapped, :]

    def split_per_chromosome(self) -> None:
        """Function return list of dfs with CpGs and MAPINFO per chromosome"""
        cpg_per_chr = []

        for chr_ in tqdm(self.manifest["CHR"].unique()):
            cgs = self.manifest[self.manifest["CHR"] == chr_]["MAPINFO"].astype(int)
            cpg_per_chr.append(cgs)

        self.cpg_per_chr = cpg_per_chr

    @staticmethod
    def export_to_csv(df: pd.DataFrame, path) -> None:
        path = create_path(path, name="CpG_in_island.csv")
        df.to_csv(path)

    @staticmethod
    def annotate(df: pd.DataFrame, annotations: pd.DataFrame) -> pd.DataFrame:
        """Function to annotate base CpG and nearby CpG"""
        conseq_cpg_annotations = annotations.loc[df.index, :]
        conseq_cpg_annotations = conseq_cpg_annotations.add_prefix('Nearby CpG: ')
        df = pd.concat((df, conseq_cpg_annotations), axis=1, sort=False)

        df = df.reset_index()
        df = df.set_index("Base CpG")

        base_cpg_annotations = annotations.loc[df.index, :]
        base_cpg_annotations = base_cpg_annotations.add_prefix('Base CpG: ', )
        df = pd.concat((df, base_cpg_annotations), axis=1, sort=False)

        df = df.sort_values(by=["Base CpG: CHR", "Base CpG: MAPINFO"])
        return df
