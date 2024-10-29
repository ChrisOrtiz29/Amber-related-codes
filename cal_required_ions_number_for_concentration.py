#!/usr/bin/env python
import os
import sys
from math import ceil, floor

# import load_settings

ion_concentration = float(sys.argv[1])
N_water = int(sys.argv[2])
N_neu_Na = int(sys.argv[3])
N_neu_Cl = int(sys.argv[4])


def count_ion(water_number, ion_concentration):
    # water_number     : Number of water molecules
    # ion_concentration: Ion concentration (mol/L)

    water_number = float(water_number)
    ion_concentration = float(ion_concentration)

    """ Constant definition """

    # Avogadro's Number (Number of molecules or atoms per mol)
    Avogadro = 6.02214129e23

    # Water density (g/cm^3, human body temperature)
    WaterDensity = 0.994

    # Mass of one molecule of water (g/mol)
    MassWater = 18

    # Atomic mass unit (kg/u)
    # Amu = 1.660538921E-27

    """ Constant definition """

    # Mass of all water (g)
    total_water_mass = (water_number / Avogadro) * MassWater

    # Number of liters of water (L)
    liters_of_water = total_water_mass / WaterDensity * 1e-3

    # Number of ions (mol)
    mols_of_ions = liters_of_water * ion_concentration

    # Number of ions molecules
    number_of_ions = mols_of_ions * Avogadro
    ion_number = int(round(number_of_ions))

    return ion_number


N_ion = count_ion(N_water, ion_concentration)
T_neu = "Na+" if N_neu_Na > N_neu_Cl else "Cl-"  # type of ion for neutraliazation
N_neu = N_neu_Na if N_neu_Na > N_neu_Cl else N_neu_Cl

if N_ion < N_neu:
    print(
        "The number of ions for concentraiton is smaller than the number of ions for neutralization. No extra ions required!"
    )
    print("0 0")

else:
    if T_neu == "Na+":
        N_total_Na = N_ion + int(ceil(N_neu / 2))
        N_total_Cl = N_ion - int(floor(N_neu / 2))
    elif T_neu == "Cl-":
        N_total_Na = N_ion - int(floor(N_neu / 2))
        N_total_Cl = N_ion + int(ceil(N_neu / 2))
    else:
        print('"cal_required_ions_number_for_concentration.py" fail!')

        sys.exit(-1)

print(N_total_Na, N_total_Cl)
