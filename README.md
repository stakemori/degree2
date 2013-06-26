# degree2

A Sage package for computation of degree 2 Siegel modular forms

# Installation
Here, we install the package in "~/sage\_packages/".

1. First install [Sage](http://www.sagemath.org/).

2. Create the directory and clone the repository.

    ```sh
    mkdir ~/sage_packages
    cd ~/sage_packages
    git clone https://github.com/stakemori/degree2.git
    ```

3. Put the following lines to "~/.sage/init.sage" and replace the
   string between double quotes in the second line by the absolute
   path of "~/sage\_packages".
   If you are a Mac user, it is "/Users/your\_username/sage_packages".

    ```python
    import sys
    sys.path.append("/absolute/path/to/sage_packages")
    from degree2.all import *
    ```
# Basic Usage

* Siegel-Eisenstein of degree two can be obtained by the function
  `EisensteinSeries_degree_2`. Siegel-Eisenstein series is normalized
  so that the constant term is one.

    ```python
    sage: es4 = EisensteinSeries_degree_2(4, 10)
    sage: es4.prec
    10
    sage: es4.fourier_coefficient(2, 1, 3)
    2903040
    ```
  The last line means that the Fourier coefficient of `es4` at the
  half integral symmetric matrix ![alt text](./images/mat1.png) is 2903040.
  The third lines means that `es4` knows Fourier coefficients of
  Siegel-Eisenstein series of weight 4
  at half integral matrices whose trace are less than or equal to 10.
