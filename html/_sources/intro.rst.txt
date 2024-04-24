Introduction
------------

pidp_tools is a collection of functions and classes that provide visualizations and analysis tools for GlueX data. This package is designed for use with `PID Playground <https://duberii.github.io/pid-playground/>`_, which is an educational resource designed to teach physics students about the particle identification process. As such, these tools may or may not be useful in experimental analysis. This package is designed to be used in Google Colaboratory, and using it on your local machine is not recommended.


Installation
------------

pidp_tools is hosted on PyPI, so you can install the latest version (in Google Colab) as follows::

    !pip install pidp_tools

After this, you can import the library. It is recommended to import the library as follows, so you don't need to worry about the prefixes::
    
    from pidp_tools import *

This way, you don't need to worry about whether a function is part of the analysis module or the visualization module.


