import pandas as pd
from src.utils import get_dictonary_of_scores_maveDB, get_list_to_add_in_dataframe


class MAVE_GOLD_STANDARD:
    @staticmethod
    def get_dataframe_for_mave_gs_data(mave_gs_file_path, column_names):
        """
        Given the column names and the path of the mave db gold standard data
        it returns the pandas dataframe.
        """

        dictionary_of_data = get_dictonary_of_scores_maveDB(mave_gs_file_path)

        rows = []
        for protein_name, value in dictionary_of_data.items():
            row = [protein_name] + get_list_to_add_in_dataframe(value)
            rows.append(row)
        return pd.DataFrame(rows, columns=column_names)