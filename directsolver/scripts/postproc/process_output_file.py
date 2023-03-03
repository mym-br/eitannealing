#%%
import re
import pathlib
import pandas as pd


def get_log_files(folder: str):
    return [
        f for f in pathlib.Path(folder).iterdir() if f.suffix == ".txt"
        and "compilation" not in f.name and "malhaadpt" in f.name
    ]


def get_kernel_times(file: pathlib.Path) -> list[tuple[str]]:
    with open(file) as f:
        lines = [
            line for line in f.readlines() if re.match(".+breakdown", line)
        ]
    return [
        re.findall(
            "^Average (.+) iteration time breakdown: (\S+) \(.+\) (\S+) \(.+\) (\S+) \(remaining\) (\S+) \(total\)",
            l)[0] for l in lines
    ]


METHOD_MAPPING = {
    'serial': 'serial',
    'cuda': 'cuda',
    'Consolidated Cuda CG': 'consolidatedcuda',
    'Consolidated Cuda CG with cooperative groups': 'consolidatedcudaCG',
    'cublas/cusparse': 'cublas',
}


def process_kernel_by_method(method: str, kernel_times: list[tuple[str]]):
    COLUMNS_SUFFIXES = ["triangularsolver", "spmv", "remaining", "total"]
    times = [
        list(map(lambda x: float(x), t[1:])) for t in kernel_times
        if t[0] == method
    ]
    columns = [f"{METHOD_MAPPING[method]}{s}" for s in COLUMNS_SUFFIXES]
    return pd.DataFrame(times, columns=columns)


def process_kernel_times(kernel_times: list[tuple[str]], filename: str):
    # pd.concat([df3, df4], axis=1)
    df = pd.concat([
        process_kernel_by_method(m, kernel_times)
        for m in METHOD_MAPPING.keys()
    ],
                   axis=1)
    df["Filename"] = filename
    return df


def average_kernel_times(df_times):
    return df_times.groupby('Filename').mean()


def process_single_kernel(file: pathlib.Path):
    kernel_times_raw = get_kernel_times(file)
    return average_kernel_times(
        process_kernel_times(kernel_times_raw, file.stem))


def process_all_kernels(folder: str):
    log_files = get_log_files(folder)
    return pd.concat([process_single_kernel(f) for f in log_files])
