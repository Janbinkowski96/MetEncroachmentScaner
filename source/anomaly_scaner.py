import multiprocessing

import pandas as pd
from tqdm import tqdm

from source.utils import rename_axis, interval_size

class AnomalyScaner:
    def __init__(self, mynorm: pd.DataFrame, mynorm_std: pd.DataFrame, search_range: int, annotations: pd.DataFrame):
        self.search_range = search_range
        self.annotations = annotations
        self.mynorm_std = mynorm_std
        self.mynorm = mynorm
    
    def search_for_anomalies(self, chromosome_cpg_collection: pd.Series) -> pd.DataFrame:
        batch = []
        
        for cpg in tqdm(chromosome_cpg_collection.index, desc=multiprocessing.current_process().name):
            
            # Set first CpG as base
            base_cpg_loc = chromosome_cpg_collection.loc[cpg]
            
            # Calculate interval size
            start_loc, stop_loc = interval_size(self.search_range, base_cpg_loc)

            # Find CpG in nearby
            cpg_in_range = chromosome_cpg_collection[(chromosome_cpg_collection < stop_loc) & 
                                                     (chromosome_cpg_collection > start_loc)].to_frame()
            
            # If there are CpGs in specific interval around base-CpG, find their annotations
            if not cpg_in_range.empty:

                cpg_in_range = rename_axis(cpg_in_range, {"index": "Nearby CpG"})
                cpg_in_range = rename_axis(cpg_in_range, {"MAPINFO": "Nearby CpG position"})

                # Get methylation levels and std for nearby CpG
                cpg_in_range["Nearby CpG b-values"] = self.mynorm.loc[cpg_in_range.index, :]
                cpg_in_range["Nearby CpG b-values std"] = self.mynorm_std.loc[cpg_in_range.index, :]
                
                # Set index to base cg and add base_loc
                cpg_in_range["Base CpG"] = cpg
                cpg_in_range["Base CpG position"] = base_cpg_loc

                # Check methylation level for base
                cpg_in_range["Base CpG b-values"] = float(self.mynorm.loc[cpg, :])
                cpg_in_range["Base CpG b-values std"] = float(self.mynorm_std.loc[cpg, :])

                # Count real diff met level
                cpg_in_range["Delta"] = cpg_in_range["Nearby CpG b-values"] - cpg_in_range["Base CpG b-values"] 
                cpg_in_range["Delta"] = cpg_in_range["Delta"]
                
                # Count distance
                cpg_in_range["Distance"] = cpg_in_range["Nearby CpG position"] - cpg_in_range["Base CpG position"]
                cpg_in_range = cpg_in_range[(cpg_in_range["Delta"].abs().round(2) > 0.3)]
                if not cpg_in_range.empty:
                    if cpg_in_range.shape[0] > 1:
                        cpg_in_range = cpg_in_range[cpg_in_range["Distance"] == cpg_in_range["Distance"].min()]
                    batch.append(cpg_in_range)
                    
                        
        batch = pd.concat(batch, axis=0, sort=False)
        return batch
