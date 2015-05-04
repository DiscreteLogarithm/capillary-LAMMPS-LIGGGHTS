# capillary-LAMMPS-LIGGGHTS

Simple capillary force model for granular matter in LAMMPS/LIGGGHTS

## Installation

* Put the patch file in your top-level LAMMPS directory, where the `LICENSE` and `README` files are.

* Apply the patch by typing the following command in your top-level LAMMPS directory:

        patch -p1 < patch.capillary

* Then include this package by typing the following command in `src` directory:

        make yes-user-capillary

* Then build again as usual.

## Uninstallation

* Use following command in `src` directory:

        make no-user-capillary

* Then build again

