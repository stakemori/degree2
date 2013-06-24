# degree2

A Sage package for computation of degree 2 Siegel modular forms

# Installation
Here, we install the package in "~/sage\_packages/".

1. First install [Sage](http://www.sagemath.org/).

2. Create the directory and clone the repository.

    ```sh
    mkdir ~/sage_packages
    cd ~/sage_packages
    git clone git://github.com/stakemori/degree2.git
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
