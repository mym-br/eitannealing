import pathlib
import pandas as pd


#%%
def get_instance_executions(compilation_file_path: str) -> dict[str, int]:
    if not pathlib.Path(compilation_file_path).exists():
        return {}

    df = pd.read_csv(
        compilation_file_path,
        delimiter='\t',
        header=None,
    )
    df['Filename'] = df.iloc[:, 0].apply(lambda x: pathlib.Path(x).stem)
    return df.groupby(['Filename']).size().to_dict()


def count_total_exections(executions: dict[str, int], instances: set[str]):
    return sum([e for (i, e) in executions.items() if i in instances])


# %%
