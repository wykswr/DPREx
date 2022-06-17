from pathlib import Path
from typing import Union
import pandas as pd


def get_rna_property(bed: Union[str, Path], work_dir: Union[str, Path]) -> pd.DataFrame:
    """
    :param bed: the path of the input repeats file
    :param work_dir: the directory to save the intermediate file
    :return: the dataframe, with column names: rloop, st_energy, ed_energy, DNA_repair_0, DNA_repair_1, DNA_repair_2,
    DNA_repair_3. The order of the rows must be consistent with the input repeats file.
    """
    pass
