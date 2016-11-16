# GAMS-MATLAB

This repository contains a small, but helpful collction of MATLAB functions that allow **reading** output, **writing** input and **transforming** data to facilitate working with GAMS models.

## Dependencies

1. MATLAB (all 201x versions should work) 
2. GAMS (must include the [GDXMRW][1] interface, providing the `rgdx` and `wgdx` functions) 

## Installation

1. Put file `GAMS.m` anywhere into your MATLAB path, e.g. in `C:\Users\%username%\Documents\MATLAB`
2. Verify that the main directory of your local GAMS installation is in the MATLAB path (e.g. `C:\GAMS\win64\24.5`, depending on your installed version). To test this, you can type `help rgdx` within MATLAB. This should yield a help message, not an error.

## Further reading 

For now, refer to the [documentation](doc/index.rst) in this directory and the corresponding code examples in the subdirectories.

  
## Copyright

Copyright (C) 2011-2013, 2016 Johannes Dorfner

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


[1]: http://www.gams.com/help/topic/gams.doc/tools/gdxmrw/index.html
