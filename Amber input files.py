#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os, sys


# In[ ]:


#Energy Minimization Stage 1
with open("min0.in", "w") as file:
    file.write(
"""Type of Simulation Being Done: Energy Minimization, Stage1, RESTRAINING ALL HEAVY ATOMS EXCEPT WATER AND IONS,
 &cntrl
  ntxo=2, IOUTFM=1, !NetCDF Binary Format.
  imin=1, !Energy Minimization
  maxcyc=5000, !Total Minimization Cycles to be run. Steepest Decent First, then Conjugate Gradient Method if ncyc < maxcyc
  ncyc=2500, !Number of Steepest Decent Minimization Steps to run before switching to Conjugate Gradient
  cut=10.00000, !Cut Off Distance for Non-Bounded Interactions
  igb=0, !No Generalized Born
  ntp=0, !No pressure scaling (Default)
  ntf=1, !Complete Interactions are Calculated
  ntc=1, !SHAKE is NOT performed, DEFAULT
  ntpr=10, !Every 10 steps, energy information will be printed in human-readable form to files "mdout" and "mdinfo"
  ntwx=100, !Every 100 steps, the coordinates will be written to the mdcrd file
  ntwr=100, !Every 100 steps during dynamics, the restart file will be written, ensuring that recovery from a crash will not be so painful. #If ntwr < 0, a unique copy of the file, "restrt_<nstep>", is written every abs(ntwr) steps
  ntr=1, !Turn ON (Cartesian) Restraints
  restraintmask="!(@H=|:WAT|@Na+|@Cl-)", !Atoms to be Restrained are specified by a restraintmask
  restraint_wt=25.00, !Force Constant for Restraint, kcal/(mol * A^2)
/
  """)


# In[ ]:


#Energy Minimization Stage 2
with open("min1.in", "w") as file:
    file.write(
"""Type of Simulation Being Done: Energy Minimization, Stage2
 &cntrl
  ntxo=2, IOUTFM=1, !NetCDF Binary Format.
  imin=1, !Energy Minimization
  maxcyc=5000, !Total Minimization Cycles to be run. Steepest Decent First, then Conjugate Gradient Method if ncyc < maxcyc
  ncyc=2500, !Number of Steepest Decent Minimization Steps to run before switching to Conjugate Gradient
  cut=10.00000, !Cut Off Distance for Non-Bounded Interactions
  igb=0, !No Generalized Born
  ntp=0, !No pressure scaling (Default)
  ntf=1, !Complete Interactions are Calculated
  ntc=1, !SHAKE is NOT performed, DEFAULT
  ntpr=10, !Every 10 steps, energy information will be printed in human-readable form to files "mdout" and "mdinfo"
  ntwx=10, !Every 100 steps, the coordinates will be written to the mdcrd file
  ntwr=100, !Every 100 steps during dynamics, the restart file will be written, ensuring that recovery from a crash will not be so painful. #If ntwr < 0, a unique copy of the file, "restrt_<nstep>", is written every abs(ntwr) steps
  ntr=1, !Turn ON (Cartesian) Restraints
  restraintmask="@N,CA,C,O", !Atoms to be Restrained are specified by a restraintmask
  restraint_wt=10.00, !Force Constant for Restraint, kcal/(mol * A^2)
/
""")


# In[ ]:


#Energy Minimization Stage 3
with open("min2.in", "w") as file:
    file.write(
"""Type of Simulation Being Done: Energy Minimization, Stage3
 &cntrl
  ntxo=2, IOUTFM=1, !NetCDF Binary Format.
  imin=1, !Energy Minimization
  maxcyc=5000, !Total Minimization Cycles to be run. Steepest Decent First, then Conjugate Gradient Method if ncyc < maxcyc
  ncyc=2500, !Number of Steepest Decent Minimization Steps to run before switching to Conjugate Gradient
  cut=10.00000, !Cut Off Distance for Non-Bounded Interactions
  igb=0, !No Generalized Born
  ntp=0, !No pressure scaling (Default)
  ntf=1, !Complete Interactions are Calculated
  ntc=1, !SHAKE is NOT performed, DEFAULT
  ntpr=10, !Every 10 steps, energy information will be printed in human-readable form to files "mdout" and "mdinfo"
  ntwx=10, !Every 100 steps, the coordinates will be written to the mdcrd file
  ntwr=100, !Every 100 steps during dynamics, the restart file will be written, ensuring that recovery from a crash will not be so painful. #If ntwr < 0, a unique copy of the file, "restrt_<nstep>", is written every abs(ntwr) steps
  ntr=1, !Turn ON (Cartesian) Restraints
  restraintmask="@CA", !Atoms to be Restrained are specified by a restraintmask
  restraint_wt=5.00, !Force Constant for Restraint, kcal/(mol * A^2)
/
""")

#Heating Stage

with open("heat.in", "w") as file:
    file.write(
"""Type of Simulation Being Done: Heating,
 &cntrl
  ntxo=2, IOUTFM=1, !NetCDF Binary Format.
  imin=0, !MD Simulation
  irest=0, !Do NOT Restart the Simulation; instead, run as a NEW Simulation
  ig=-1, !Pseudo-random number seed is changed with every run.
  ntx=1, !Coordinates, but no Velocities, will be read. Formatted (ASCII) coordinate file is expected
  nstlim=5000000, !Number of MD-steps to be performed. Default 1.
  dt=0.00200, !The time step (psec). Recommended MAXIMUM is .002 if SHAKE is used, or .001 if SHAKE is NOT used
  cut=10.00000, !Cut Off Distance for Non-Bounded Interactions
  igb=0, !No Generalized Born
  ntp=0, !No pressure scaling (Default)
  ntf=2, !Bond Interactions involving H-atoms omitted
  ntc=2, !Bonds involving Hydrogen are Constrained
  ntpr=10000, !Every 100 steps, energy information will be printed in human-readable form to files "mdout" and "mdinfo"
  ntwx=10000, !Every 100 steps, the coordinates will be written to the mdcrd file
  ntwr=10000, !Every 1000 steps during dynamics, the restart file will be written, ensuring that recovery from a crash will not be so painful. #If ntwr < 0, a unique copy of the file, "restrt_<nstep>", is written every abs(ntwr) steps
  ntt=3, !Use Langevin Dynamics with the Collision Frequency GAMA given by gamma_ln,
  gamma_ln=2.00000, !Collision Frequency, ps ^ (-1)
  get_ipython().system('temp0=310.00000, !Reference temperature at which the system is to be kept')
  tempi=0.00000, !Initial Temperature
  ntr=1, !Turn ON (Cartesian) Restraints
  restraintmask="@CA", !Atoms to be Restrained are specified by a restraintmask
  restraint_wt=5.00, !Force Constant for Restraint, kcal/(mol * A^2)
/
""")


init_step = 0
total_step = 500000
step_inc  = 8928
temp_inc  = 5
init_temp = 50.0
norm_temp = 310.0
max_temp  = 320.0

heat_in = "heat.in"
topfile  = "temp_inc.in"
    
TempFile = open(topfile, 'w')
while (init_temp < 320.0):
    TempFile.writelines("&wt type='TEMP0', istep1=" + str(init_step) + ", istep2=" + str(init_step+step_inc) + ", value1=" + str(init_temp) + ", value2=" + str(init_temp+5) +", /" + "\n")
    init_step = init_step + step_inc
    init_temp = init_temp + 5

while (max_temp > 310.0):
    if max_temp == 310.0:
        TempFile.writelines("&wt type='TEMP0', istep1=" + str(init_step) + ", istep2=" + str(total_step) + ", value1=" + str(norm_temp) + ", value2=" + str(norm_temp) +", /" + "\n")
    else:
        TempFile.writelines("&wt type='TEMP0', istep1=" + str(init_step) + ", istep2=" + str(init_step+step_inc) + ", value1=" + str(max_temp) + ", value2=" + str(max_temp-5) +", /" + "\n")
        init_step = init_step + step_inc
        max_temp = max_temp - 5
        
TempFile.writelines("&wt type='END' /")
TempFile.close()

#merge the heating input file and the step by step temp increment
merge_file = "cat " + heat_in + " " + topfile + " > heating.in"
os.system(merge_file)


# NVT Equilibration

with open("equi_NVT.in", "w") as file:
    file.write(
"""Type of Simulation Being Done: NVT Equilibration,
 &cntrl
  ntxo=2, IOUTFM=1, !NetCDF Binary Format.
  imin=0, !MD Simulation
  irest=0, !Do NOT Restart the Simulation; instead, run as a NEW Simulation
  ig=-1, !Pseudo-random number seed is changed with every run.
  ntx=1, !Coordinates, but no Velocities, will be read. Formatted (ASCII) coordinate file is expected
  nstlim=5000000, !Number of MD-steps to be performed. Default 1.
  dt=0.00200, !The time step (psec). Recommended MAXIMUM is .002 if SHAKE is used, or .001 if SHAKE is NOT used
  cut=10.00000, !Cut Off Distance for Non-Bounded Interactions
  igb=0, !No Generalized Born
  ntp=0, !No pressure scaling (Default)
  ntf=2, !Bond Interactions involving H-atoms omitted
  ntc=2, !Bonds involving Hydrogen are Constrained
  ntpr=10000, !Every 100 steps, energy information will be printed in human-readable form to files "mdout" and "mdinfo"
  ntwx=50000, !Every 100 steps, the coordinates will be written to the mdcrd file
  ntwr=10000, !Every 1000 steps during dynamics, the restart file will be written, ensuring that recovery from a crash will not be so painful. #If ntwr < 0, a unique copy of the file, "restrt_<nstep>", is written every abs(ntwr) steps
  ntt=3, !Use Langevin Dynamics with the Collision Frequency GAMA given by gamma_ln,
  gamma_ln=2.00000, !Collision Frequency, ps ^ (-1)
  temp0=310.00000, !Reference temperature at which the system is to be kept
  tempi=310.00000, !Initial Temperature
  ntr=1, !Turn ON (Cartesian) Restraints
  restraintmask="@CA", !Atoms to be Restrained are specified by a restraintmask
  restraint_wt=2.00, !Force Constant for Restraint, kcal/(mol * A^2)
/
 """)

# NPT Equilibration

with open("equi_NPT.in", "w") as file:
    file.write(
"""Type of Simulation Being Done: NPT Equilibration,
 &cntrl
  ntxo=2, IOUTFM=1, !NetCDF Binary Format.
  imin=0, !MD Simulation
  irest=0, !Do NOT Restart the Simulation; instead, run as a NEW Simulation
  ig=-1, !Pseudo-random number seed is changed with every run.
  ntx=5, !Coordinates and Velocities will be read; a formatted (ASCII) coordinate file is expected.
  nstlim=10000000, !Number of MD-steps to be performed. Default 1.
  dt=0.00200, !The time step (psec). Recommended MAXIMUM is .002 if SHAKE is used, or .001 if SHAKE is NOT used
  cut=10.00000, !Cut Off Distance for Non-Bounded Interactions
  igb=0, !No Generalized Born
  ntp=1, !MD with isotropic position scaling
  ntf=2, !Bond Interactions involving H-atoms omitted
  ntc=2, !Bonds involving Hydrogen are Constrained
  ntpr=1000, !Every 1000 steps, energy information will be printed in human-readable form to files "mdout" and "mdinfo"
  ntwx=50000, !Every 10000 steps, the coordinates will be written to the mdcrd file
  ntwr=10000, !Every 10000 steps during dynamics, the restart file will be written, ensuring that recovery from a crash will not be so painful. #If ntwr < 0, a unique copy of the file, "restrt_<nstep>", is written every abs(ntwr) steps
  ntt=3, !Use Langevin Dynamics with the Collision Frequency GAMA given by gamma_ln,
  gamma_ln=2.00000, !Collision Frequency, ps ^ (-1)
  temp0=310.00000, !Reference temperature at which the system is to be kept
  tempi=310.00000, !Initial Temperature
  pres0=1.01300, !Reference Pressure (in units of bars, where 1 bar = 0.987 atm) at which the system is maintained
  ntr=1, !Turn ON (Cartesian) Restraints
  restraintmask="@CA", !Atoms to be Restrained are specified by a restraintmask
  restraint_wt=2.00, !Force Constant for Restraint, kcal/(mol * A^2)
/
""")

# Production Run
with open("md_NPT.in", "w") as file:
    file.write(
"""Type of Simulation Being Done: Production Run,
 &cntrl
  ntxo=2, IOUTFM=1, !NetCDF Binary Format.
  imin=0, !MD Simulation
  irest=0, !Start new Simulation;
  ig=-1, !Pseudo-random number seed is changed with every run.
  ntx=5, !Coordinates and Velocities will be read; a formatted (ASCII) coordinate file is expected.
  nstlim=10000000, !Number of MD-steps to be performed. Default 1.
  dt=0.00200, !The time step (psec). Recommended MAXIMUM is .002 if SHAKE is used, or .001 if SHAKE is NOT used
  cut=10.00000, !Cut Off Distance for Non-Bounded Interactions
  igb=0, !No Generalized Born
  ntp=1, !MD with isotropic position scaling
  ntf=2, !Bond Interactions involving H-atoms omitted
  ntc=2, !Bonds involving Hydrogen are Constrained
  ntpr=1000, !Every 1000 steps, energy information will be printed in human-readable form to files "mdout" and "mdinfo"
  ntwx=50000, !Every 10000 steps, the coordinates will be written to the mdcrd file
  ntwr=10000, !Every 10000 steps during dynamics, the restart file will be written, ensuring that recovery from a crash will not be so painful. #If ntwr < 0, a unique copy of the file, "restrt_<nstep>", is written every abs(ntwr) steps
  ntt=3, !Use Langevin Dynamics with the Collision Frequency GAMA given by gamma_ln,
  gamma_ln=2.00000, !Collision Frequency, ps ^ (-1)
  temp0=310.00000, !Reference temperature at which the system is to be kept
  tempi=310.00000, !Initial Temperature
  pres0=1.01300, !Reference Pressure (in units of bars, where 1 bar = 0.987 atm) at which the system is maintained
/
 """)


# In[46]:


#Submission File
topfile = "cmp.parm7"
coordfile = "cmp.crd"
min_prefix = "min"
heat_prefix = "heating"
equiNVT_prefix = "equi_NVT"
equiNPT_prefix = "equi_NPT"
prod_prefix = "md_NPT"
runfile = "run.sh"

GPU_CARD_ID ="1"
AMBERPATH = "/mnt/Tsunami_HHD/opt/amber20"

RunFile = open(runfile,'w')
RunFile.writelines("export CUDA_VISIBLE_DEVICES=" + GPU_CARD_ID + "\n")
RunFile.writelines("source " + AMBERPATH + "/amber.sh" + "\n\n")
RunFile.writelines("AMBERPATH=" + AMBERPATH + "/bin" + "\n\n")

RunFile.writelines("$AMBERPATH/pmemd.cuda_SPFP -O -i " + min_prefix + "0.in -p " +  topfile + " -c " + coordfile + " -o " + min_prefix + "0.out -r " + min_prefix + "0.rst -ref " + coordfile + " -x " + min_prefix + "0.nc -inf " + min_prefix + "0.info" + "\n\n")
RunFile.writelines("$AMBERPATH/pmemd.cuda_SPFP -O -i " + min_prefix + "1.in -p " +  topfile + " -c " + min_prefix + "0.rst -o " + min_prefix + "1.out -r " + min_prefix + "1.rst -ref " + min_prefix + "0.rst -x " + min_prefix + "1.nc -inf " + min_prefix + "1.info" + "\n\n")
RunFile.writelines("$AMBERPATH/pmemd.cuda_SPFP -O -i " + min_prefix + "2.in -p " +  topfile + " -c " + min_prefix + "1.rst -o " + min_prefix + "2.out -r " + min_prefix + "2.rst -ref " + min_prefix + "1.rst -x " + min_prefix + "2.nc -inf " + min_prefix + "2.info" + "\n\n")
RunFile.writelines("$AMBERPATH/pmemd.cuda_SPFP -O -i " + heat_prefix + ".in -p " +  topfile + " -c " + min_prefix + "2.rst -o " + heat_prefix + ".out -r " + heat_prefix + ".rst -ref " + min_prefix + "2.rst -x " + heat_prefix + ".nc -inf " + heat_prefix + ".info" + "\n\n")
RunFile.writelines("$AMBERPATH/pmemd.cuda_SPFP -O -i " + equiNVT_prefix + ".in -p " +  topfile + " -c " + heat_prefix + ".rst -o " + equiNVT_prefix + ".out -r " + equiNVT_prefix + ".rst -ref " + heat_prefix + ".rst -x " + equiNVT_prefix + ".nc -inf " + equiNVT_prefix + ".info" + "\n\n")
RunFile.writelines("$AMBERPATH/pmemd.cuda_SPFP -O -i " + equiNPT_prefix + ".in -p " +  topfile + " -c " + equiNVT_prefix + ".rst -o " + equiNPT_prefix + ".out -r " + equiNPT_prefix + ".rst -ref " + equiNVT_prefix + ".rst -x " + equiNPT_prefix + ".nc -inf " + equiNPT_prefix + ".info" + "\n\n")
RunFile.writelines("$AMBERPATH/pmemd.cuda_SPFP -O -i " + prod_prefix + ".in -p " +  topfile + " -c " + prod_prefix + ".rst -o " + prod_prefix + ".out -r " + prod_prefix + ".rst -ref " + equiNPT_prefix + ".rst -x " + prod_prefix + ".nc -inf " + prod_prefix + ".info" + "\n\n")
RunFile.close()






