# IsoGen
IsoGen is a toolbox for predicting isotope distributions in mass spectrometry data.

It includes both an absolute FFT-based calculation and a neural network prediction.

Pretrained models are included for both peptides and RNA based on either average mass or sequence.

The FFT methods are absolute and are limited only by the accuracy of the data you put in. They are a little faster, especially on larger species.

The NN methods are very accurate and can be faster on smaller species. The primary advantage of these is that they can be retrained on non-standard isotope distributions.

## Installation

Precompiled libraries for Windows and Linux are provided in the "bin" directory.

For other platforms, you can build the library from source using CMake. An example compile script is provided.

The only dependency is FFTW, but we have statically compiled it into the library, so no additional installation is necessary unless you want to build your own.

## Usage

A simple test script is provided. It can be run with:

'''
isogen_test.exe -mass 10000
'''

This will create an isotope distribution printout for a mass of 10 kDa.

However, most people will want to use functions in the DLL or Python wrapper.

## License

IsoGen is released under the BSD 3-Clause License. See LICENSE for details.

PLEASE CITE THIS SOFTWARE IN ANY PUBLICATIONS THAT USE IT (publication to follow).

## Contact

If you have any questions, please email mtmarty@utexas.edu or open a ticket on GitHub.



