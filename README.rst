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
::

    Usage: svbench [OPTIONS] REFERENCE_VCF QUERY_VCFS...
    
    Options:
      --include PATH            Include regions .bed file
      --pass-only               Assess only PASS variants
      --slop INTEGER            Add intervals +/- slop around breakpoints
                                [default: 250]
      --min-size-ref INTEGER    Min SV length  [default: 0]
      --min-size-query INTEGER  Min SV length  [default: 30]
      --no-duplicates           Don't quantify duplicate true positives
      --version                 Show the version and exit.
      --help                    Show this message and exit.




Benchmark the number of query SVs in a reference/truth set. Results are printed to stderr.::

    svbench truthset.vcf query1.vcf


Multiple query vcfs can also be analysed::

    svbench truthset.vcf query1.vcf query2.vcf ...


API Usage
---------
For a general tutorial on using the API see the ipython notebook; `svbench_tutorial.ipynb`
