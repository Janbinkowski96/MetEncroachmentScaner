import pandas as pd
from tqdm import tqdm
import plotly.graph_objects as go

from .island_collector import IslandCollector

class IslandAnalyzer:
    def __init__(self, cpgs: list, mynorm: pd.DataFrame, manifest: pd.DataFrame, islands_annotations: pd.DataFrame):
        self.cpgs = cpgs
        self.mynorm = mynorm
        self.manifest = manifest
        self.base_cpgs_index = cpgs.index
        self.islands_annotations = islands_annotations
        self.related_islands = []
        
    def find_related_island(self) -> None:
        for base_cpg in tqdm(self.base_cpgs_index):
            
            nearby_cpg = self.cpgs.loc[base_cpg, "Nearby CpG"]
            cpg_loc = self.manifest.loc[base_cpg, ["MAPINFO", "CHR"]]
            mapinfo = int(cpg_loc["MAPINFO"])
            chr_ = f"chr{cpg_loc['CHR']}"
            
            related_islands = self.islands_annotations[(self.islands_annotations["Island and shore start"] <= mapinfo) & 
                                                 (self.islands_annotations["Island and shore end"] >= mapinfo) &
                                                 (self.islands_annotations.index == chr_)]
            
            island = IslandCollector(chr_, base_cpg, nearby_cpg, related_islands)
            self.related_islands.append(island)
                    
    def plot(self, raw_mynorm: pd.DataFrame) -> None:
        for island in tqdm(self.related_islands):
            
            base_cpg, nearby_cpg = island.get_anomaly_pair

            if not island.empty:
                
                island_start, island_end = island.coordinates
                cpgs_in_island = self.manifest[(self.manifest["MAPINFO"] >= island_start) & 
                                               (self.manifest["MAPINFO"] <= island_end)]
                
                cpgs_in_island = cpgs_in_island["MAPINFO"].to_frame()
                cpgs_in_island = pd.concat((cpgs_in_island, raw_mynorm), axis=1, sort=False).dropna()

                fig = go.Figure()

                fig.add_vrect(x0 = island_end - 2000,
                          x1 = island_end,
                          fillcolor='rgb(8,81,156)', opacity=0.2, 
                          name="N shore", layer="below", 
                          line_width=0)
                
                fig.add_vrect(x0 = island_start - 2000,
                          x1 = island_start,
                          fillcolor='rgb(8,81,156)', opacity=0.2, 
                          name="S Shore", layer="below", 
                          line_width=0)
                
                for cpg in cpgs_in_island.index:
                    
                    if cpg == base_cpg or cpg == nearby_cpg:
                        color = "rgb(255,0,0)"
                    
                    else:
                        color = "grey"
                        
                    y = cpgs_in_island.drop("MAPINFO", axis=1).loc[cpg, :].values
                    x = [int(cpgs_in_island.loc[cpg, "MAPINFO"])] * len(y)
                    
                    fig.add_trace(go.Box(x=x, boxmean=True, boxpoints='all',
                                         y=y, fillcolor=color, marker_color=color, 
                                         line_color=color, name=f"{cpg}"))
                
                fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='grey', tickformat=".0f", title_text="MAPINFO")
                fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='grey', title_text='Methylation level [b-values]', range=[0, 1])
                fig.update_layout(showlegend=False)
                fig.write_image(f"TEST/Plots/{base_cpg}.png")