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

See the inline documentation in the code or see the manual generated using `tbd`
or[online at rtd](https://classirot.readthedocs.org).

For a start, you can run the calculation inputs provided in `examples`:
```
cmiclassirot -o example.h5 <example>
cmiclassirot-plot example.h5
```

See also [Release Notes.md] and [Contributors.md].

For registered developers, updated [documentation of development
versions](https://cmi.pages.desy.de/internal/software/cmiclassirot) is available
at DESY.


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
