#%%
import os
import logging
import pathlib
import argparse
import subprocess
from typing import Optional
from alive_progress import alive_bar


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


#%%
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("directsolver", help="path to directsolver executable")
    parser.add_argument("eit", help="path to eit mtx files folder")
    parser.add_argument("suitesparse",
                        help="path to suitesparse mtx files folder")
    parser.add_argument("rhs", help="path to right hand side eit files folder")
    parser.add_argument("results", help="path to save results files")
    parser.add_argument("--repetitions",
                        help="number of executions per instance")
    parser.add_argument("--compilation", help="name for the compilation file")
    args = parser.parse_args()

    repetitions = args.repetitions if args.repetitions else DEFAULT_REPETITIONS

    eit_mtx_files = preprocess_mtx_files(args.eit,
                                         MTX_INSTANCES_SETS["eit"],
                                         set_name='eit mtx')
    eit_rhs_files = preprocess_mtx_files(args.rhs,
                                         MTX_INSTANCES_SETS["eit"],
                                         is_rhs=True,
                                         set_name='eit rhs')
    eit_files = list(zip(sorted(eit_mtx_files), sorted(eit_rhs_files)))
    suitesparse_files = preprocess_mtx_files(args.suitesparse,
                                             MTX_INSTANCES_SETS["suitesparse"],
                                             set_name='suitsparse')

    results_folder = args.results

    logging.basicConfig(filename=os.path.join(results_folder, 'log.log'),
                        encoding='utf-8',
                        level=logging.DEBUG)

    compilation_filename = args.compilation if args.compilation else DEFAULT_COMPILATION_FILENAME

    # run eit instances
    total_errors = 0
    with alive_bar(repetitions * len(eit_files)) as bar:
        bar.text = f'-> Running EIT mtx instances (step 1/3), please wait...'
        for mtx, rhs in eit_files:
            for _ in range(repetitions):
                try:
                    run_instance(executable=args.directsolver,
                                 mtx=mtx,
                                 res=None,
                                 maxits=200,
                                 results_folder=results_folder,
                                 compilation_filename=compilation_filename,
                                 rhs=rhs)
                except subprocess.CalledProcessError as err:
                    logging.error(
                        f'{err}\nstdout:{err.stdout}\nstderr:{err.stderr}\n\n')
                    total_errors = total_errors + 1
                finally:
                    bar()
    print(f'Step 1 concluded with {total_errors} errors')

    # run suitesparse instances
    total_errors = 0
    with alive_bar(repetitions * len(suitesparse_files)) as bar:
        bar.text = f'-> Running suitesparse mtx instances (step 2/3), please wait...'
        for mtx in suitesparse_files:
            for _ in range(repetitions):
                try:
                    run_instance(
                        executable=args.directsolver,
                        mtx=mtx,
                        res=0.001,
                        maxits=2000,
                        results_folder=results_folder,
                        compilation_filename=compilation_filename,
                    )
                except subprocess.CalledProcessError as err:
                    logging.error(
                        f'{err}\nstdout:{err.stdout}\nstderr:{err.stderr}\n\n')
                    total_errors = total_errors + 1
                finally:
                    bar()
    print(f'Step 2 concluded with {total_errors} errors')


if __name__ == "__main__":
    main()