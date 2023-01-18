#%%
import pathlib
import pandas as pd
from functools import reduce


#%%
def get_compilation_columns(file_name: str):
    compilation_columns = {
        "compilation": [
            "Filename", "N", "NNZ", "Serial Analyzer", "Serial Execution",
            "Serial Iterations", "Cuda Analyzer", "Cuda Execution",
            "Cuda Iterations", "Consolidated Cuda Analyzer",
            "Consolidated Cuda Execution", "Consolidated Cuda Iterations",
            "Coop. Groups Cuda Analyzer", "Coop. Groups Cuda Execution",
            "Coop. Groups Cuda Iterations", "CUBLAS Analyzer",
            "CUBLAS Execution", "CUBLAS Iterations"
        ],
        "compilation_cusolver": [
            "Filename", "N", "NNZ", "Cusolver Analyzer", "Cusolver Execution",
            "Cusolver Iterations"
        ],
        "compilation_pardiso": [
            "Filename", "N", "NNZ", "Pardiso Analyzer", "Pardiso Execution",
            "Pardiso Iterations"
        ],
    }
    return compilation_columns[pathlib.Path(file_name).stem]


def get_sorted_averages(df):
    df = df.groupby('Filename').mean().sort_values(by="N")
    return df


def parse_single_compilation(compilation_file_path: str):
    df = pd.read_csv(compilation_file_path,
                     delimiter='\t',
                     header=None,
                     names=get_compilation_columns(compilation_file_path))
    df = get_sorted_averages(df)
    return df


def filter_eit_suitsparse(df):
    return df[df.index.str.contains("eit")], df[~df.index.str.contains("eit")]


def set_index_stem(df):
    df_stem = df.copy()
    df_stem.index = df_stem.index.map(lambda f: pathlib.Path(f).stem)
    return df_stem


def parse_compilation_files(compilation_file_paths: list[str]):
    parsed_dfs = map(lambda x: parse_single_compilation(x),
                     compilation_file_paths)
    df = reduce((lambda x, y: x.join(y.drop(columns=['N', 'NNZ']))),
                parsed_dfs)
    return [set_index_stem(d) for d in filter_eit_suitsparse(df)]


#%%
# import matplotlib.pyplot as plt

# df_eit, df_suitsparse = parse_compilation_files([
#     "D:/aksato/ResearchResults/EITResults/EIT2022/MtxCudaTests6-GoogleCloud-2022/results/compilation.txt",
#     "D:/aksato/ResearchResults/EITResults/EIT2022/MtxCudaTests6-GoogleCloud-2022/results/compilation_cusolver.txt",
#     "D:/aksato/ResearchResults/EITResults/EIT2022/MtxCudaTests6-GoogleCloud-2022/results/compilation_pardiso.txt",
# ])
# #%%
# # Figure 8 plot
# ax = df_eit.plot(x='N',
#                  y=[
#                      "Serial Execution", "Consolidated Cuda Execution",
#                      "Coop. Groups Cuda Execution"
#                  ],
#                  kind='line')
# ax.set_yscale('log')
# plt.legend(['Eigen3', 'ParallelPCG', 'CuFPPCGSolver'])
# plt.show()
# # %%
# # Figure 9 plot
# df_speedups = df_eit.copy()
# df_speedups["Consolidated Cuda Speedup"] = df_speedups[
#     "Serial Execution"] / df_speedups["Consolidated Cuda Execution"]
# df_speedups["Coop. Groups Cuda Speedup"] = df_speedups[
#     "Serial Execution"] / df_speedups["Coop. Groups Cuda Execution"]
# df_speedups.plot(x='N',
#                  y=["Consolidated Cuda Speedup", "Coop. Groups Cuda Speedup"],
#                  kind='line')
# plt.legend(['ParallelPCG', 'CuFPPCGSolver'])
# plt.show()
# # %%
# # Figure 11 plot
# ax = df_eit.plot(
#     x="N",
#     y=["Serial Execution", "Coop. Groups Cuda Execution", "CUBLAS Execution"],
#     kind="line")
# # ax.set_yscale('log')
# plt.legend(["Eigen3", "CuFPPCGSolver", "NVIDIA"])
# plt.show()
# # %%
# # Figure 12 plot
# ax = df_eit.plot(
#     x="N",
#     y=["Serial Analyzer", "Coop. Groups Cuda Analyzer", "CUBLAS Analyzer"],
#     kind="line")
# # ax.set_yscale('log')
# plt.legend(["Eigen3", "CuFPPCGSolver", "NVIDIA"])
# plt.show()
# #%%