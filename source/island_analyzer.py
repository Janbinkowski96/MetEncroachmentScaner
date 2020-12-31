import pandas as pd
from tqdm import tqdm
import plotly.graph_objects as go

from .island_collector import IslandCollector


class IslandAnalyzer:
    def __init__(self, cpgs: pd.DataFrame, mynorm: pd.DataFrame, manifest: pd.DataFrame):
        self.cpgs = cpgs
        self.mynorm = mynorm
        self.manifest = manifest
        self.base_cpgs_index = cpgs.index
        self.related_islands = []
        self.islands_annotations = None

    def select_islands(self):
        manifest = self.manifest[["Relation_to_UCSC_CpG_Island", "UCSC_CpG_Islands_Name", "CHR"]].dropna()
        manifest["UCSC_CpG_Islands_Name"] = manifest["UCSC_CpG_Islands_Name"].map(lambda x: x.split(":")[1]).to_frame()

        manifest["Island start"] = manifest["UCSC_CpG_Islands_Name"].map(lambda x: int(x.split("-")[0]))
        manifest["Island end"] = manifest["UCSC_CpG_Islands_Name"].map(lambda x: int(x.split("-")[1]))

        manifest["Island and shore start"] = manifest["Island start"] - 2000
        manifest["Island and shore end"] = manifest["Island end"] + 2000

        manifest = manifest[manifest["Relation_to_UCSC_CpG_Island"] == "Island"]
        self.islands_annotations = manifest[["CHR", "Island start", "Island end",
                                             "Island and shore start", "Island and shore end"]]

    def find_related_island(self) -> None:
        for base_cpg in tqdm(self.base_cpgs_index):
            nearby_cpg = self.cpgs.loc[base_cpg, "Nearby CpG"]
            cpg_loc = self.manifest.loc[base_cpg, ["MAPINFO", "CHR"]]
            mapinfo = int(cpg_loc["MAPINFO"])
            chr_ = str(cpg_loc['CHR'])

            related_islands = self.islands_annotations[(self.islands_annotations["Island and shore start"] <= mapinfo) &
                                                       (self.islands_annotations["Island and shore end"] >= mapinfo) &
                                                       (self.islands_annotations["CHR"] == chr_)]

            island = IslandCollector(chr_, base_cpg, nearby_cpg, related_islands)
            self.related_islands.append(island)

    def plot(self, raw_mynorm: pd.DataFrame) -> None:
        for island in tqdm(self.related_islands):

            if not island.empty:
                base_cpg, nearby_cpg = island.get_anomaly_pair
                island_start, island_end, island_chr_ = island.coordinates

                cpgs_in_island = self.manifest[(self.manifest["MAPINFO"] >= island_start) &
                                               (self.manifest["MAPINFO"] <= island_end) &
                                               (self.manifest["CHR"] == island_chr_)]

                cpgs_in_island = cpgs_in_island["MAPINFO"].to_frame()
                cpgs_in_island = pd.concat((cpgs_in_island, raw_mynorm), axis=1, sort=False).dropna()

                fig = go.Figure()

                fig.add_vrect(x0=island_start + 2000,
                              x1=island_end - 2000,
                              fillcolor='rgb(8,81,156)', opacity=0.2,
                              name="Island", layer="below",
                              line_width=0)

                for cpg in cpgs_in_island.index:

                    if cpg == base_cpg or cpg == nearby_cpg:
                        color = "rgb(255,0,0)"

                    else:
                        color = "grey"

                    y = cpgs_in_island.drop("MAPINFO", axis=1).loc[cpg, :].values
                    x = [int(cpgs_in_island.loc[cpg, "MAPINFO"])] * len(y)

                    fig.add_trace(go.Scatter(x=x, mode='markers',
                                             y=y, fillcolor=color, marker_color=color,
                                             line_color=color, name=f"{cpg}"))

                fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='grey', tickformat=".0f", title_text="MAPINFO",
                                 range=[island_start, island_end])
                fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='grey',
                                 title_text='Methylation level [b-values]', range=[0, 1])

                fig.update_layout(showlegend=False)
                fig.write_image(f"TEST/Plots/{base_cpg}.png")
