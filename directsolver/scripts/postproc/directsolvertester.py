#%%
import os
import logging
import pathlib
import argparse
import subprocess
from typing import Optional
from alive_progress import alive_bar
from check_executions import get_instance_executions, count_total_exections


#%%
def get_mtx_files(folder: str):
    return [f for f in pathlib.Path(folder).iterdir() if f.suffix == ".mtx"]


def filter_mtx_files(files: list[pathlib.Path],
                     instances: set[str],
                     rhs=False):

    def get_instance_name(f: pathlib.Path) -> str:
        return f.stem if not rhs else f.stem[:-1]

    files_set = set([get_instance_name(f) for f in files])
    diff = list(instances - files_set)
    if len(diff) != 0:
        raise RuntimeError(
            f"could not find the following eit files: {', '.join(diff)}.")

    return [f for f in files if get_instance_name(f) in instances]


#%%

DEFAULT_REPETITIONS = 10

DEFAULT_COMPILATION_FILENAME = 'compilation.txt'

MTX_INSTANCES_SETS = {
    "eit":
    set([
        "malhaadpt-1-15", "malhaadpt-1-40", "malhaadpt-1-65", "malhaadpt-1-90",
        "malhaadpt-1-115", "malhaadpt-1-140", "malhaadpt-1-165",
        "malhaadpt-1-190", "malhaadpt-1-215", "malhaadpt-1-240",
        "malhaadpt-1-265", "malhaadpt-1-290", "malhaadpt-1-315",
        "malhaadpt-1-340", "malhaadpt-1-365", "malhaadpt-1-390",
        "malhaadpt-1-415", "malhaadpt-1-440", "malhaadpt-1-465",
        "malhaadpt-1-490", "malhaadpt-1-515", "malhaadpt-1-540",
        "malhaadpt-1-565", "malhaadpt-1-590", "malhaadpt-1-615",
        "malhaadpt-1-640", "malhaadpt-1-665", "malhaadpt-1-690"
    ]),
    "suitesparse":
    set([
        "Trefethen_20b", "Trefethen_20", "mesh1e1", "mesh1em1", "mesh1em6",
        "bcsstk01", "bcsstk02", "nos4", "bcsstk04", "bcsstk22", "lund_b",
        "lund_a", "bcsstk05", "mesh3e1", "mesh3em5", "mesh2e1", "mesh2em5",
        "mhdb416", "bcsstm07", "nos5", "494_bus", "662_bus", "nos6", "685_bus",
        "msc00726", "nos7", "nos3", "bcsstk08", "1138_bus", "bcsstk27",
        "bcsstm12", "nasa2146", "Chem97ZtZ", "mhd3200b", "mhd4800b",
        "crystm01", "bcsstk16", "s1rmq4m1", "Kuu", "Muu", "aft01", "fv1",
        "fv3", "fv2", "bundle1", "ted_B", "ted_B_unscaled", "t2dah_e",
        "crystm02", "Pres_Poisson", "Dubcova1", "gyro_m", "bodyy4", "bodyy5",
        "bodyy6", "crystm03", "wathen100", "wathen120", "jnlbrng1", "torsion1",
        "obstclae", "minsurfo", "gridgena"
    ])
}

DEFAULT_DIRECT_SOLVER = "cufppcgsolver.exe"
DIRECT_SOLVERS_EXECUTABLES = [
    DEFAULT_DIRECT_SOLVER, "cusolver.exe", "pardiso.exe"
]


def preprocess_mtx_files(folder: str,
                         expected_instances: set[str],
                         is_rhs=False,
                         set_name=""):
    try:
        mtx_files = get_mtx_files(folder)
        mtx_files = filter_mtx_files(mtx_files, expected_instances, rhs=is_rhs)
    except RuntimeError as err:
        raise RuntimeError(f"error checking {set_name} files: {err}")
    return mtx_files


def run_instance(executable: str,
                 mtx: pathlib.Path,
                 res: Optional[float],
                 maxits: float,
                 results_folder: str,
                 compilation_filename: str,
                 rhs: Optional[pathlib.Path] = None):
    args = [
        executable,
        mtx,
        f"--maxits={maxits}",
        f"--output={os.path.join(results_folder, mtx.stem)}",
        f"--compilation={os.path.join(results_folder, compilation_filename)}",
    ]
    if res:
        args.append(f"--res={res}")
    if rhs:
        args.append(f"-b{rhs}")

    result = subprocess.run(args, capture_output=True, text=True)
    result.check_returncode()
    with open(os.path.join(results_folder, f"{mtx.stem}.txt"), "a") as f:
        f.write(result.stdout)


def find_executables(folder: str):
    print(os.getcwd())
    print([f for f in pathlib.Path(folder).iterdir()])
    executables = [
        f for f in pathlib.Path(folder).iterdir()
        if f.suffix == ".exe" and f.name in DIRECT_SOLVERS_EXECUTABLES
    ]
    return executables


def batch_run_instances(instance_files: list[pathlib.Path],
                        executables: list[pathlib.Path],
                        repetitions: int,
                        total_executions: int,
                        executions: dict[str, int],
                        exec_compilation_filename: dict[str, str],
                        results_folder: str,
                        bar_text: str,
                        maxits: Optional[int] = 200,
                        res: Optional[float] = None):
    total_errors = 0
    with alive_bar(repetitions * len(instance_files) -
                   total_executions) as bar:
        bar.text = bar_text
        for mtx, rhs in instance_files:
            for _ in range(repetitions - executions.get(mtx.stem, 0)):
                try:
                    for executable in executables:
                        run_instance(
                            executable=executable,
                            mtx=mtx,
                            res=res,
                            maxits=maxits,
                            results_folder=results_folder,
                            compilation_filename=exec_compilation_filename[
                                executable.name],
                            rhs=rhs)
                except subprocess.CalledProcessError as err:
                    logging.error(
                        f'{err}\nstdout:{err.stdout}\nstderr:{err.stderr}\n\n')
                    total_errors = total_errors + 1
                finally:
                    bar()
    return total_errors


#%%
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("directsolver",
                        help="path with directsolver executables")
    parser.add_argument("eit", help="path to eit mtx files folder")
    parser.add_argument("suitesparse",
                        help="path to suitesparse mtx files folder")
    parser.add_argument("rhs", help="path to right hand side eit files folder")
    parser.add_argument("results", help="path to save results files")
    parser.add_argument("--repetitions",
                        help="number of executions per instance",
                        type=int,
                        default=10)
    parser.add_argument("--res",
                        help="residual for the CG solvers",
                        type=float)
    parser.add_argument("--maxits",
                        help="maximum number of iterations for CG solvers",
                        type=int)
    parser.add_argument("--compilation", help="name for the compilation file")
    parser.add_argument('--skip-cusolver',
                        help="skip execution of cusolver solver",
                        action=argparse.BooleanOptionalAction)
    parser.add_argument('--skip-pardiso',
                        help="skip execution of pardiso solver",
                        action=argparse.BooleanOptionalAction)
    args = parser.parse_args()

    repetitions = args.repetitions

    res = args.res

    eit_mtx_files = preprocess_mtx_files(args.eit,
                                         MTX_INSTANCES_SETS["eit"],
                                         set_name='eit mtx')
    eit_rhs_files = preprocess_mtx_files(args.rhs,
                                         MTX_INSTANCES_SETS["eit"],
                                         is_rhs=True,
                                         set_name='eit rhs')
    eit_files = list(zip(sorted(eit_mtx_files), sorted(eit_rhs_files)))

    suitesparse_mtx_files = preprocess_mtx_files(
        args.suitesparse,
        MTX_INSTANCES_SETS["suitesparse"],
        set_name='suitsparse')
    suitesparse_files = list(
        zip(suitesparse_mtx_files, [None] * len(suitesparse_mtx_files)))

    results_folder = args.results

    executables = find_executables(args.directsolver)
    if args.skip_cusolver:
        executables = [i for i in executables if i.name != "cusolver.exe"]
    if args.skip_pardiso:
        executables = [i for i in executables if i.name != "pardiso.exe"]

    logging.basicConfig(level=logging.INFO,
                        encoding='utf-8',
                        format="%(asctime)s [%(levelname)s] %(message)s",
                        handlers=[
                            logging.FileHandler(
                                os.path.join(results_folder, 'log.log')),
                            logging.StreamHandler()
                        ])

    compilation_filename = args.compilation if args.compilation else DEFAULT_COMPILATION_FILENAME
    exec_compilation_filename: dict[str, str] = {
        exec: f"{pathlib.Path(compilation_filename).stem}_{exec[:-4]}.txt"
        if exec != DEFAULT_DIRECT_SOLVER else compilation_filename
        for exec in DIRECT_SOLVERS_EXECUTABLES
    }

    # check for existing results
    executions = get_instance_executions(
        os.path.join(results_folder, compilation_filename))
    total_eit_executions = count_total_exections(executions,
                                                 MTX_INSTANCES_SETS["eit"])
    total_suitesparse_executions = count_total_exections(
        executions, MTX_INSTANCES_SETS["suitesparse"])

    # run eit instances
    eit_total_errors = batch_run_instances(
        eit_files,
        executables,
        repetitions,
        total_eit_executions,
        executions,
        exec_compilation_filename,
        results_folder,
        '-> Running EIT mtx instances (step 1/3), please wait...',
        res=res,
        maxits=args.maxits)
    logging.info(f'Step 1 concluded with {eit_total_errors} errors')

    # run suitesparse instances
    suitesparse_total_errors = batch_run_instances(
        suitesparse_files,
        executables,
        repetitions,
        total_suitesparse_executions,
        executions,
        exec_compilation_filename,
        results_folder,
        '-> Running suitesparse mtx instances (step 2/3), please wait...',
        res,
        maxits=args.maxits)
    logging.info(f'Step 2 concluded with {suitesparse_total_errors} errors')


if __name__ == "__main__":
    main()
