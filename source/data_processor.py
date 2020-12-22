import pandas as pd
from tqdm import tqdm

class DataProcessor:
    def __init__(self, manifest_type):
        self.manifest_type = manifest_type
        self.cg_per_chr = None
        self.mynorm_std = None
        self.manifest = None
        self.mynorm = None
                
    def set_manifest(self) -> None:
        if self.manifest_type == "E":
            self.manifest = pd.read_csv("resources/EPIC/EPIC.csv", index_col=0, low_memory=False)
        else:
            self.manifest = pd.read_csv("resources/450K/450K.csv", index_col=0, low_memory=False)
    
    def load_mynorm(self, mynorm_path: str) -> None:
        mynorm = pd.read_csv(mynorm_path, index_col=0, encoding="latin1")
        self.mynorm = mynorm.mean(axis=1).to_frame(name="beta-values")
        self.mynorm_std = mynorm.std(axis=1).to_frame(name="beta-values std")
            
    def select_cpg(self) -> None:
        mynorm_cpg = set(self.mynorm.index)
        overlapped = set.intersection(mynorm_cpg, set(self.manifest.index))
        self.manifest = self.manifest.loc[overlapped, :]
    
    def load_islands(self):
        return pd.read_csv("resources/hg19/Annotated_islands.csv", index_col=0)
     
    def split_per_chromosome(self) -> None:
        """Function return list of dfs with CpGs and MAPINFO per chromosome"""
        cpg_per_chr = []
    
        for chr_ in tqdm(self.manifest["CHR"].unique()):
            cgs = self.manifest[self.manifest["CHR"] == chr_]["MAPINFO"].astype(int)
            cpg_per_chr.append(cgs)
    
        self.cpg_per_chr = cpg_per_chr
        
    @staticmethod
    def annotate(df: pd.DataFrame, annotations: pd.DataFrame) -> pd.DataFrame:
        """Function to annotate base CpG and nearby CpG"""
        conseq_cpg_annotations = annotations.loc[df.index, :]
        conseq_cpg_annotations = conseq_cpg_annotations.add_prefix('Nearby CpG: ')
        df = pd.concat((df, conseq_cpg_annotations), axis=1, sort=False)
        
        df = df.reset_index()
        df = df.set_index("Base CpG")
        
        base_cpg_annotations = annotations.loc[df.index, :]
        base_cpg_annotations = base_cpg_annotations.add_prefix('Base CpG: ',)
        df = pd.concat((df, base_cpg_annotations), axis=1, sort=False)

        df = df.sort_values(by=["Base CpG: CHR", "Base CpG: MAPINFO"])
        return df
