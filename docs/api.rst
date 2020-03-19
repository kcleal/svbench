SVBench - an interface for classifying variant calls
====================================================


SVBench is a python module that makes it easy to parse, benchmark and classify variant calls in .vcf or
.csv format.

This page provides a quick introduction to the command-line interface for applying
models to variant calls, and the API for parsing and manipulating variant files.

To apply a predefined model to your variant calls, svbench can be invoked from the command-line::

    $ svbench model.pkl variants.vcf > variants.classified.vcf

The model.pkl file is a predefined SVBench object which is typically composed of parsing instructions and
a classifier. Internally, the model parses the input file variants.vcf and applies the classifier to each
variant. Typically this will result in each variant being assigned a probability value, or class label,
although the exact output can be flexible, and can include optional filtering, for example.

To inspect the contents of the model


.. autofunction:: svbench.Col
