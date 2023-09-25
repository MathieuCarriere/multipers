# multipers

**Warning: this package is deprecated and no longer maintained. I strongly suggest users to download and run [the MMA package](https://github.com/DavidLapous/multipers-signed-measure) for analyzing multiparameter persistence modules. In particular, the multiparameter persistence summaries of the article "Multiparameter persistence image for topological machine learning" (NeurIPS 2020), available [here](https://papers.nips.cc/paper/2020/hash/fdff71fcab656abfbefaabecab1a7f6d-Abstract.html) can be obtained with `SimplexTree.persistence_approximation(...).get_barcodes(...).to_multipers()` with the MMA API.**

Package for computing multiparameter persistence summaries. This is the code accompanying the article "Multiparameter persistence image for topological machine learning" (NeurIPS 2020), available [here](https://papers.nips.cc/paper/2020/hash/fdff71fcab656abfbefaabecab1a7f6d-Abstract.html). The notebooks contain code for reproducing experiments (`experiments.ipynb`) and testing simple cases (`visualization.ipynb`). 

The computation of multiparameter persistence images involves matching bars of different barcodes together. This can be done either by using your own function (e.g., the matching induced by the Wasserstein distance), or by using [vineyards](https://www.mrzv.org/publications/vineyards/socg06/), see the documentation in the source code. In the latter case, you will have to compile the code in the `dionysus_vineyards` folder into a Python library, which you can easily do by running the first cell in `experiments.ipynb` or `visualization.ipynb`. 

All files in the `dionysus_vineyards` folder (except for `dionysus_vineyards.hpp`, `dionysus_vineyards.pyx` and `setup.py`) are extracted from the Dionysus library, and can also be accessed [here](https://hg.mrzv.org/Dionysus/file/tip/include).
A custom implementation of the vineyards algorithm from the authors (which is also in C++) is also available upon request (by sending an email to mathieu.carriere3@gmail.com). 
