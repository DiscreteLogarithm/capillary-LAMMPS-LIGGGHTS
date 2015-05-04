# capillary-LAMMPS-LIGGGHTS

Simple capillary force model for granular matter in LAMMPS/LIGGGHTS

## Installation

First apply the patch:

    patch -p1 < patch.capillary

Then enable this package by:

    make yes-user-capillary

Then compile again as usual.

## Uninstallation

Use following command:

    make no-user-capillary

Then compile again

