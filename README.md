# Classical-physics simulations of laser alignment

Classical-mechanics simulations of the rotational dynamics of rigid molecules
and nanoparticles in laser-electric fields.


## Installation

Change your working directory to the root directory of the package and run
```
pip install .
```

For development work `pip install --editable .` might be a useful alternative.


## Documentation

See the online documentation of the latest release on
[GitHub Pages](https://cfel-cmi.github.io/CMIclassirot).

For a start, you can also inspect and run the calculation inputs provided in
`examples`:
```
cmiclassirot -o example.h5 <example>
cmiclassirot-plot example.h5
```

See also the [Release Notes](Release%20Notes.md) and the
[Contributors file](Contributors.md).

Documentation of the development version is available on
[DESY's gitlab pages](https://cmi.pages.desy.de/internal/software/cmiclassirot).


## Citation

When you use this program in scientific work, please cite it as Muhamed Amin,
Jean-Michel Hartmann, Amit K. Samanta, Jochen Küpper, "Laser-induced alignment
of nanoparticles and macromolecules for single-particle-imaging applications",
[arXiv:2306.05870](https://arxiv.org/abs/2306.05870)




<!-- Put Emacs local variables into HTML comment
Local Variables:
coding: utf-8
fill-column: 80
End:
-->
