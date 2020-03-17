=======
SVBench
=======

SVBench is a library for benchmarking variant calls vs a reference dataset.
Also an interface is provided for training, saving and applying `sklearn <https://scikit-learn.org/stable/>`_, or other
models to custom models to SV call sets.


Installation
------------
Install using::

    $ python setup.py install
    # Or
    $ pip install -r requirements.txt; pip install .

Requires Python>=3.6, packages needed are listed in requirements.txt.


Usage
-----
For a general tutorial on using the API see the ipython notebook; `svbench_tutorial.ipynb`

To apply a model using the command line interface::

    $ svbench my_model.pkl svs.vcf > svs.classified.vcf


Available models
----------------
Models for a few popular SV callers are provided in the SVBench/models folder


API documentation
-----------------
Under construction
