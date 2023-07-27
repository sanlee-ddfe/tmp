# OpenLB - Open Source Lattice Boltzmann Code

The OpenLB project is a C++ package for the implementation of lattice Boltzmann
methods adressing a vast range of tansport problems.

## Dependencies

The only mandatory external dependency of OpenLB is GNU Make and a C++ compiler
with C++17 support. This includes all reasonably recent versions of GCC, Clang
and Intel C++.

GPU support depends on Nvidia CUDA 11 or later.

## Installation

1. Adjust the `config.mk` file to fit your local compiler environment
   (examples for some common configurations are available in `config/`
2. Build the embedded dependencies and core library using `make`
3. Switch to any of the examples and compile it using `make`

## Documentation

A comprehensive user guide is available in `doc/userGuide`.

Papers featuring OpenLB are collected at [1].

Up-to-date Doxygen documentation can be generated via `make doxygen`
or accessed at [2].

[1]: https://www.openlb.net/articles/
[2]: https://www.openlb.net/DoxyGen/html/index.html

## Community

The OpenLB forum [3] is an open discussion board for all aspects of LBM and OpenLB.
Feel free to post any problems, questions, suggestions or contributions.

You can also reach us via mail at info@openlb.net

There is a yearly one-week Spring School [4] where you can learn about LBM and
OpenLB directly from the developer team and invited guest lecturers.

A list of all present and past contributors is available at [5].

[3]: https://www.openlb.net/forum/
[4]: https://www.openlb.net/spring-school-2022/
[5]: https://www.openlb.net/authors/

## How to cite OpenLB

Standardized citation metadata is available in the `CITATION.cff` file.

To cite OpenLB in general instead of a specific release we suggest:

```
@article{OpenLB2020,
  title = {OpenLB - Open Source Lattice Boltzmann Code},
  year = {2020},
  issn = {08981221},
  doi = {10.1016/j.camwa.2020.04.033},
  journal = {Computers \& Mathematics with Applications},
  author = {Krause, Mathias J. and Kummerl{\"a}nder, Adrian and Avis, Samuel J. and Kusumaatmaja, Halim and Dapelo, Davide and Klemens, Fabian and Gaedtke, Maximilian and Hafen, Nicolas and Mink, Albert and Trunk, Robin and Marquardt, Jan E. and Maier, Marie-Luise and Haussmann, Marc and Simonis, Stephan},
}
```

This article is available as open access [6].

[6]: https://doi.org/10.1016/j.camwa.2020.04.033

## Code Formatting

The basic formatting rules are described by the `.editorconfig` file [7]
for automatic application in a wide variety of text editors and IDEs.

[7]: https://editorconfig.org/

## License

OpenLB is provided as open source under the terms of the GNU GPL v2 license.

See `LICENSE` for details.
