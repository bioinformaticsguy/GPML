from pathlib import Path

from src.constants import TABLES_DIRECTORY_PATH, TABLE_FORMAT


class TableGenerator:
    @staticmethod
    def generate_table(df,
                       file_name,
                       table_path=TABLES_DIRECTORY_PATH,
                       table_format=TABLE_FORMAT,
                       columns=None):
        """
        Generate a table from a DataFrame.

        Parameters:
            df (DataFrame): The DataFrame containing the data.
            file_name (str): The name of the file to save the table as.
            table_path (str): The directory path to save the table.
            table_format (str): The format of the table file.
            columns (list, optional): List of column names to include in the table. If None, all columns are included.

        Returns:
            None
        """

        table_path = Path(table_path) / f"{file_name}.{table_format}"

        if columns:
            df = df[columns]

        df.to_csv(table_path, index=False)

        print(f"Table saved to: {table_path}")