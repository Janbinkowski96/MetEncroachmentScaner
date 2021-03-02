from dataclasses import dataclass
from typing import Tuple

import pandas as pd


@dataclass
class IslandCollector:
    chr_: str
    base_cpg: str
    nearby_cpg: str
    related_islands: pd.DataFrame

    @property
    def get_anomaly_pair(self) -> Tuple[int, int]:
        return self.base_cpg, self.nearby_cpg

    @property
    def coordinates(self) -> Tuple[int, int, str]:
        island_start = int(self.related_islands["Island and shore start"].min())
        island_end = int(self.related_islands["Island and shore end"].max())

        return island_start, island_end, self.chr_

    @property
    def empty(self) -> bool:

        if self.related_islands.empty:
            return True

        else:
            return False
