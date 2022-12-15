#%%
# imports
import os
import scipy
import pathlib
# from decimal import Decimal
from typing import Literal, TypedDict
#%%
# input arguments. TODO: receive via command line
folder = "D:/aksato/ResearchResults/EITResults/EIT2022/MtxCudaTests5-2022"

#%%
#  define types

CG_METHOD = Literal["serial", "cuda", "ccuda", "ccudacg", "cusparse"]


class InstanceData(TypedDict):
    folder: str
    instance: str
    methods: list[CG_METHOD]


ResultData = dict[str, dict[CG_METHOD, float]]


#%%
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


#%%
# function to calculate single residual
def compute_residual(x, A, rhs):
    return scipy.linalg.norm(A * x - rhs)


#%%
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


#%%
# function to compute all instance residuals
def compute_residuals(folder: str) -> list[ResultData]:
    return [get_instance_residuals(i) for i in get_instances_data(folder)]


# %%
# compute residuals for input folder
instances_res = compute_residuals(folder)
instances_res
# %%
