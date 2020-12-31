from dataclasses import dataclass
import pandas as pd


@dataclass
class CpGCollector:
    base_cpg_id: str
    search_range: int
    chromosome_cpg_collection: pd.Series
    start: int = 0
    end: int = 0
    cpg_in_range: pd.DataFrame = pd.DataFrame()

    def interval_size(self) -> None:
        base_cpg_loc = self.chromosome_cpg_collection[self.base_cpg_id]
        self.start, self.end = base_cpg_loc, base_cpg_loc + self.search_range

    def find_nearby_cpg(self) -> None:
        cpg_in_range = self.chromosome_cpg_collection[(self.chromosome_cpg_collection > self.start)
                                                      & (self.chromosome_cpg_collection < self.end)]
        self.cpg_in_range = cpg_in_range.to_frame()

    @property
    def empty(self) -> bool:
        if self.cpg_in_range.empty:
            return True
        else:
            return False
