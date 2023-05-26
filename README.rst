=======
SVBench
=======

SVBench is a library for benchmarking variant calls vs a reference dataset.


Installation
------------
Install using::

    pip install svbench

Or::

    pip install -r requirements.txt; pip install .

Requires Python>=3.6, packages needed are listed in requirements.txt.


CLI Usage
---------

Benchmark the number of query SVs in a reference/truth set. Results are printed to stderr.::

    svbench truthset.vcf query1.vcf


Multiple query vcfs can also be analysed::

    svbench truthset.vcf query1.vcf query2.vcf ...


API Usage
---------
For a general tutorial on using the API see the ipython notebook; `svbench_tutorial.ipynb`
