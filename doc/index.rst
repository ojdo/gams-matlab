GAMS and how to use it from MATLAB
==================================

:Maintainer: Johannes Dorfner, <johannes.dorfner@tum.de>
:Organization: `Institute for Renewable and Sustainable Energy Systems`_,
               Technische Universität München
:Version: |version|
:Date: |today|
:Copyright:
  The corresponding code is licensed under the `GNU General Public License 3.0`_.
  This documentation is licensed under a `Creative Commons Attribution 4.0 
  International`_ license.

Introduction
------------

This manual consists of two major sections:

The first section is a short **GAMS introduction** with a motivation for
an external data handling component. The basic language elements *sets,
parameters, variables* and *equations* are briefly explained, so that
readers with experience in mathematical modelling should be able to write new 
models.

The second section then introduces GAMS.m, a class written for
MATLAB that provides static functions that allow creating GAMS data
structures, reading and writing GDX files. Finally, the concepts
*entities* and *timeseries* are introduced and how they can be used for
rapid model development.

Contents
--------

.. toctree::
   :numbered: 
   :maxdepth: 2
   
   gams-intro
   gams-m

   
Download
--------

Get ``GAMS.m`` and the usage examples from its `GitHub repository`_.
   

.. _GitHub repository: https://github.com/ojdo/gams-matlab
.. _GNU General Public License 3.0: http://www.gnu.org/licenses/gpl-3.0
.. _Creative Commons Attribution 4.0 International: http://creativecommons.org/licenses/by/4.0/
.. _Institute for Renewable and Sustainable Energy Systems: http://www.ens.ei.tum.de/

