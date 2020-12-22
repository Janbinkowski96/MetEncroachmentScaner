import pandas as pd
from tqdm import tqdm
import plotly.graph_objects as go

class Visualize():
    def __init__(self, cpgs: list, mynorm: pd.DataFrame, manifest: pd.DataFrame, islands_annotations: pd.DataFrame):
        self.base_cpgs_index = cpgs.index
        self.cpgs = cpgs
        self.mynorm = mynorm
        self.manifest = manifest
        self.islands_annotations = islands_annotations
        
        
    def plot(self):
        for base_cpg in tqdm(self.base_cpgs_index):
            
            nearby_cpg = self.cpgs.loc[base_cpg, "Nearby CpG"]
            
            cpg_loc = self.manifest.loc[base_cpg, ["MAPINFO", "CHR"]]
            mapinfo = int(cpg_loc["MAPINFO"])
            chr_ = f"chr{cpg_loc['CHR']}"
            
            related_island = self.islands_annotations[(self.islands_annotations["CpG island start"] <= mapinfo) & 
                                                 (self.islands_annotations["CpG island end"] >= mapinfo) &
                                                 (self.islands_annotations.index == chr_)]
            if not related_island.empty:
                island_start, island_end = int(related_island["CpG island start"].min()), int(related_island["CpG island end"].max())

                cpgs_in_island = self.manifest[(self.manifest["MAPINFO"] > island_start - 2000) &
                                           (self.manifest["MAPINFO"] < island_end + 2000)]["MAPINFO"].to_frame()
                
                cpgs_in_island = pd.concat((cpgs_in_island, self.mynorm), axis=1, sort=False).dropna()
                
                fig = go.Figure()
                
                fig.add_vrect(x0 = island_start - 2000,
                          x1 = island_start  - 1,
                          fillcolor='rgb(8,81,156)', opacity=0.2, 
                          name="N shore", layer="below", 
                          line_width=0)
                
                fig.add_vrect(x0 = island_end + 1,
                          x1 = island_end + 2000,
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
                    
                    fig.add_trace(go.Box(x=x, boxpoints = "all", boxmean=True,
                             y=y, fillcolor=color, marker_color=color, line_color=color,
                             name=f"{cpg}"))
                
                fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='grey', tickformat=".0f", title_text="MAPINFO")
                fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='grey', title_text='Methylation level [b-values]', range=[0, 1])
                fig.update_layout(showlegend=False)
                fig.write_image(f"Plots/{base_cpg}.png")