# true_graphene

A Transparent, Reproducible, Usable by Others, Extensible (TRUE) simulation workflow of a graphene
slit pore.

## Notes
A small bug was discovered in the number density calculation which caused the averaging to be incorrect.
The bug has since been fixed and the changes are
reflected in the reference data.  Therefore the figure generated from
this notebook will slightly differ than the number density figure in the paper.  However, the new number density data is not
qualititatively different and does not change any of the major
conclusions of the paper.

## How is This Workflow TRUE?
All of the files needed to build and run this simulation are contained within this workflow.
Additionally, the Python functions and classes to initialize and parametrize the system are freely
available on GitHub, making the workflow transparent.  Because all of the files and code are freely
available, any user should be able to run this simulation and get similar results, making this
workflow reproducible.

Any user of this worklow is
also free to change any of the system or simulation parameters.  For instance, the temperature this
simulation is run at can easily be changed within the GROMACS mdp files contained in this repository.
Thus making the workflow usable by others and extensible.

![graphene_snapshot](https://user-images.githubusercontent.com/25011342/70189374-c994ff00-16b8-11ea-827e-3e6b7576359e.png)

## Requirements
There are several software requirements to successfully run this workflow:
- [mBuild](https://github.com/mosdef-hub/mbuild)
- [Foyer](https://github.com/mosdef-hub/foyer)
- [Matplotlib](https://matplotlib.org)
- [NumPy](https://numpy.org)
- [Pore-Builder](https://github.com/rmatsum836/Pore-Builder)
- [Jupyter](https://jupyter.org)
