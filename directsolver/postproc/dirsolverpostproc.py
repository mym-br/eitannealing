#%%
# imports
import os
import sys
import scipy
import pathlib
import argparse
import pandas as pd
# from decimal import Decimal
from typing import Literal, TypedDict
#%%
# input arguments
if hasattr(sys, 'ps1'):
    # folder = "D:/aksato/ResearchResults/EITResults/EIT2022/MtxCudaTests4-2021"
    folder = "D:/aksato/ResearchResults/EITResults/EIT2022/MtxCudaTests5-2022"
else:
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", help="path to results folder")
    args = parser.parse_args()
    folder = args.folder

#%%
#  define types
CG_METHOD = Literal["serial", "cuda", "ccuda", "ccudacg", "cusparse"]


class InstanceData(TypedDict):
    folder: str
    instance: str
    methods: list[CG_METHOD]


ResultData = dict[str, dict[CG_METHOD, float]]

#%%
# define functions


# function to find instances data in results folder
def get_instances_data(results_folder: str) -> list[InstanceData]:
    mtx_files = [
        f.name for f in pathlib.Path(os.path.join(
            results_folder, "results")).iterdir() if f.suffix == ".mtx"
    ]
    mtx_files_basenames = list(set([f.split('_')[0] for f in mtx_files]))
    return [{
        "folder":
        results_folder,
        "instance":
        i,
        "methods":
        [s.split('_')[1][:-4] for s in mtx_files if s.split('_')[0] == i]
    } for i in mtx_files_basenames]


# function to calculate single residual
def compute_residual(x, A, rhs):
    return scipy.linalg.norm(A * x - rhs)


# function to calculate all residuals for instance
def get_instance_residuals(instance_data: InstanceData) -> ResultData:
    # read input matrix and right-hand side
    A = scipy.io.mmread(
        f"{instance_data['folder']}/mtx/eit/{instance_data['instance']}.mtx")
    rhs = scipy.io.mmread(
        f"{instance_data['folder']}/rhs/{instance_data['instance']}b.mtx")
    # read x solutions
    xs = [
        scipy.io.mmread(
            f"{instance_data['folder']}/results/{instance_data['instance']}_{suffix}.mtx"
        ) for suffix in instance_data["methods"]
    ]
    # compute residuals
    res = {}
    instance_res = [compute_residual(x, A, rhs) for x in xs]
    res[instance_data['instance']] = dict(
        zip(instance_data["methods"], instance_res))
    # [f"{Decimal(e):.6E}" for e in instance_res]))
    return res


# function to compute all instance residuals
def compute_residuals(folder: str) -> list[ResultData]:
    residuals = {}
    for r in get_instances_data(folder):
        residuals.update(get_instance_residuals(r))
    return residuals


# %%
# compute residuals for input folder
instances_res = compute_residuals(folder)
instances_res
# %%
# create dataframe with results
df_data = {key: [*val.values()] for (key, val) in instances_res.items()}
df_columns = list(max(instances_res.values(), key=lambda x: len(x)).keys())
df_res = pd.DataFrame.from_dict(df_data, orient='index', columns=df_columns)
df_res = df_res.reset_index().rename(columns={'index': 'instance'})
df_res
# %%
# save residual results
df_res.to_csv(f"{pathlib.Path(folder).name}_res.csv", index=False)

# %%
# define columns for compilation processing
columns = [
    "Filename", "N", "NNZ", "Serial Analyzer", "Serial Execution",
    "Serial Iterations", "Cuda Analyzer", "Cuda Execution", "Cuda Iterations",
    "Consolidated Cuda Analyzer", "Consolidated Cuda Execution",
    "Consolidated Cuda Iterations"
]
if "ccudacg" in df_columns:
    columns += [
        "Coop. Groups Cuda Analyzer", "Coop. Groups Cuda Execution",
        "Coop. Groups Cuda Iterations"
    ]
if "cusparse" in df_columns:
    columns += ["CUBLAS Analyzer", "CUBLAS Execution", "CUBLAS Iterations"]
columns

#%%
# Read compilation file to dataframe
df_comp = pd.read_csv(os.path.join(folder, "results", "compilation.txt"),
                      delimiter='\t',
                      names=columns)
df_comp
# %%
# Calculate mean of executions
df_comp = df_comp.groupby('Filename').mean().sort_values(by="N")
df_comp
# %%
# save compilation results
df_comp.to_csv(f"{pathlib.Path(folder).name}_comp.csv")
#%%