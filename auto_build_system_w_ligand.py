#!/usr/bin/env python
from subprocess import Popen, PIPE
import sys
import os

def exec_popen(prog_name, args, cwd, env={}):
    popen_obj = Popen(
        args, cwd=cwd, stdout=PIPE, stderr=PIPE, shell=False, env=env, text=True
    )

    stdout, stderr = popen_obj.communicate()

    if stdout:
        print(stdout)

    if popen_obj.returncode != 0:
        raise Exception("Error while running '{}': {}".format(prog_name, stderr))

    return [stdout, stderr]


def exec_build_init_system(
    input_prot_pdb_file,
    input_lig_pdb_file,
    output_in_file,
    output_prmtop_file,
    output_crd_file,
    output_pdb_file,
    output_log_file,
    cwd,
    tleap_build_script="gen_model.sh",
    env={},
):
    args = [
        os.path.join(cwd,tleap_build_script),
        input_prot_pdb_file,
	input_lig_pdb_file,
        output_in_file,
        output_prmtop_file,
        output_crd_file,
        output_pdb_file,
    ]
    
    [stdout, stderr] = exec_popen("Build System Script", args, cwd, env=env)

    with open(output_log_file, "w") as output_log:
        output_log.write(stdout)


def exec_build_final_system(
    input_prot_pdb_file,
    input_lig_pdb_file,
    output_in_file,
    output_prmtop_file,
    output_crd_file,
    output_pdb_file,
    final_num_Na,
    final_num_Cl,
    cwd,
    tleap_build_script="gen_tleap_script.sh",
    env={},
):
    args = [
        tleap_build_script,
        input_prot_pdb_file,
	input_lig_pdb_file,
        output_in_file,
        output_prmtop_file,
        output_crd_file,
        output_pdb_file,
        str(final_num_Na),
        str(final_num_Cl),
    ]

    exec_popen("Build System Script", args, cwd, env=env)


def exec_get_waterions(
    log_file, cwd, get_num_waters_ions_script="get_leaplog_water_ion_numbers.py", env={}
):
    args = [get_num_waters_ions_script, log_file]

    [stdout, stderr] = exec_popen("Build System Script", args, cwd, env=env)

    [num_waters, num_Na, num_Cl] = [int(x.strip("'(),")) for x in stdout.split()]

    return [num_waters, num_Na, num_Cl]


def exec_calc_ions(
    ion_concentration,
    num_waters,
    num_Na,
    num_Cl,
    cwd,
    calc_ions_script="cal_required_ions_number_for_concentration.py",
    env={},
):
    args = [
        calc_ions_script,
        str(ion_concentration),
        str(num_waters),
        str(num_Na),
        str(num_Cl),
    ]

    [stdout, stderr] = exec_popen("Calc Ions Script", args, cwd, env=env)

    [final_num_Na, final_num_Cl] = [int(x.strip("',()")) for x in stdout.split()]

    return [final_num_Na, final_num_Cl]


def auto_build_system(input_prot_pdb_file, input_lig_pdb_file, ion_concentration=0.15, cwd='.', env={}):
    prefix = "complex"
    init_in_file = prefix + "_init.in"
    init_prmtop_file = prefix + "_init.prmtop"
    init_crd_file = prefix + "_init.crd"
    init_pdb_file = prefix + "_init.pdb"
    init_log_file = prefix + "_init.log"
    final_in_file = prefix + "_final.in"
    final_prmtop_file = prefix + "_final.prmtop"
    final_crd_file = prefix + "_final.crd"
    final_pdb_file = prefix + "_final.pdb"

    exec_build_init_system(
        input_prot_pdb_file,
	input_lig_pdb_file,
        init_in_file,
        init_prmtop_file,
        init_crd_file,
        init_pdb_file,
        init_log_file,
        cwd,
        env=env,
    )

    [num_waters, num_Na, num_Cl] = exec_get_waterions(init_log_file, cwd, env=env)

    [final_num_Na, final_num_Cl] = exec_calc_ions(
        ion_concentration, num_waters, num_Na, num_Cl, cwd, env=env
    )

    exec_build_final_system(
        input_prot_pdb_file,
	input_lig_pdb_file,
        final_in_file,
        final_prmtop_file,
        final_crd_file,
        final_pdb_file,
        final_num_Na,
        final_num_Cl,
        cwd,
        env=env,
    )


if __name__ == "__main__":
    input_prot_pdb_file = sys.argv[1]
    input_lig_pdb_file = sys.argv[2]
    ion_concentration = float(sys.argv[3])
    env = {
        'AMBERHOME':"/mnt/Tsunami_HHD/opt/amber20",
        'PATH':os.getcwd()+":/mnt/Tsunami_HHD/opt/amber20/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin",
        'LD_LIBRARY_PATH':"/mnt/Tsunami_HHD/opt/amber20/lib"
    }

    if len(sys.argv) > 4:

        for env_var in sys.argv[4:]:
            [var_name, value] = env_var.split("=")
            env[var_name] = value
    
    auto_build_system(input_prot_pdb_file, input_lig_pdb_file, ion_concentration=0.15, cwd='.', env=env)