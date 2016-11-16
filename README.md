# GAMS-MATLAB

[![Documentation Status](http://readthedocs.org/projects/gams-matlab/badge/?version=latest)](http://gams-matlab.readthedocs.io/en/latest/?badge=latest)

This repository contains a small, but helpful collction of MATLAB functions that allow **reading** output, **writing** input and **transforming** data to facilitate working with GAMS models.

## Dependencies

1. MATLAB (all 201x versions should work) 
2. GAMS (must include the [GDXMRW](http://www.gams.com/help/topic/gams.doc/tools/gdxmrw/index.html) interface, providing the `rgdx` and `wgdx` functions) 

## Installation

1. Put file `GAMS.m` anywhere into your MATLAB path, e.g. in `C:\Users\%username%\Documents\MATLAB`
2. Verify that the main directory of your local GAMS installation is in the MATLAB path (e.g. `C:\GAMS\win64\24.5`, depending on your installed version). To test this, you can type `help rgdx` within MATLAB. This should yield a help message, not an error.

## Usage

### MATLAB

    % read input data from any source or - for sake of an example -
    % create it in MATLAB programmatically
    cities = GAMS.set('cities', {'New York City' 'Los Angeles' 'Washington D.C.'})
    population = GAMS.param('population', [8046 3884 659]*1e3, cities.uels);

    % write both to a GDX file
    GAMS.putGDX('input.gdx', cities, population);

### GAMS

    $gdxin input.gdx
    Set   cities;
    Param pop(cities);
    $load cities pop
    [... remaining model code ...]
       
### MATLAB

    % after running the model using 'gams model.gmx -gdx=output.gdx'
    growth = GAMS.getGDX('output.gdx', 'growth');

    % plot it
    bar(growth.val);
    ylabel('Annual population growth (%)');
    xticklabels(growth.uels{1}); % matching labels automatically

    % or export to a spreadsheet
    GAMS.putXLS('growth.xlsx', growth

## Documentation

A pretty verbose technical documentation is [hosted on ReadTheDocs](http://gams-matlab.readthedocs.io/en/latest/), including a short syntax primer on GAMS.
  
## Copyright

Copyright (C) 2011-2014, 2016 Johannes Dorfner

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
