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
This is a transparent workflow, as all files (scripts, simulation input files, force field information) are
contained within this repository and available for anyone to use.  `true_graphene` is reproducible because
the explicit steps to set up and run the simulation are contained within `graphene-pore.ipynb`.  This
workflow is extensible, as the Python objects used to initialize the graphene slit pore system are
sufficiently documented for anyone to modify the simulation parameters.

![graphene_snapshot](https://user-images.githubusercontent.com/25011342/70189374-c994ff00-16b8-11ea-827e-3e6b7576359e.png)

## Requirements
There are several software requirements to successfully run this workflow:
- [mBuild](https://github.com/mosdef-hub/mbuild)
- [Foyer](https://github.com/mosdef-hub/foyer)
- [Matplotlib](https://matplotlib.org)
- [NumPy](https://numpy.org)
- [Pore-Builder](https://github.com/rmatsum836/Pore-Builder)
- [Jupyter](https://jupyter.org)
