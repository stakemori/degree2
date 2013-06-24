# degree2

A Sage package for computation of degree 2 Siegel modular forms

# Installation
Here, we install the pakcage in "~/sage_packages/".

1. First install [Sage](http://www.sagemath.org/).

2. Create the directory and clone the repository.

    ```sh
    mkdir ~/sage_packages
    cd ~/sage_packages
    git clone git://github.com/stakemori/degree2.git
    ```

3. Put the following lines to "~/.sage/init.sage".

    ```python
    import sys
    sys.path.append("/home/your_user_name/sage_packages")
    from degree2.all import *
    ```
