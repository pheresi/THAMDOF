# Time History Analysis of MDOF systems: THAMDOF

THAMDOF is an educational Matlab GUI for performing simplified nonlinear time history analyses (THA) of multiple-degree-of-freedom (MDOF) systems. The solution algorithm consists in Newmark's constant acceleration method, with Newton-Raphson iterations (Chopra 4th Edition, Table 16.3.3).
Here you will find the source code and a standalone executable version.

# Simplified model

THAMDOF uses a shear-building lumped-mass model of the MDOF system, with masses lumped at the floor levels. Each story lateral hysteretic model consists in a tri-linear model. P-Δ effects can be included.

The user inputs the model information through a .csv file. An example of this file is provided for reference.

<img src="Figures/BldgModel.JPG" width="250" title="Shear-building lumped-mass model"/> <img src="Figures/HystModel.JPG" width="500" title="Hysteretic model"/> 

# Ground motions

The user inputs the set of ground motion records to be analyzed and their scale factor through a .csv file. An example of this file is provided for reference. The name of each ground motion is at the first row, the number of points at the second row, the time-step (dt) at the third row, the scale factor at the fourth row, and the acceleration time series (in [g]) is from the fifth row to the end.

If the number of points specified in the fourth row is greater than the number of acceleration points in the time series, zeros are padded at the end of the record to match the specified number of points. This is particularly useful, for example, for correctly estimating residual deformations. On the other hand, if the number of points specified in the fourth row is lower than the number of acceleration points in the time series, the record is trimmed to match the specified number of points.

THAMDOF has two options for the user:
1. Analyzing every ground motion at the specified scale factor
2. Scale every ground motion in order to obtain the scale factor that results in collapse (in this option, scale factors speficied at the input file are neglected)

# Outputs

THAMDOF provides visualization of several responses. Also, it is capable of producing a video of the structural response under a given ground motion.

The user can export different outputs from THAMDOF:
1. Max and Min values of displacement, interstory displacement, interstory drift ratio (IDR), residual interstory displacement, residual IDR, absolute velocity, interstory velocity, total acceleration, and story restoring force of each floor, during every ground motion
2. Time history series of IDR and story restoring force 
3. Scale factors that produce collapse
Results are exported in a .xlsx file.

# Other information

I developed this tool initially on 2015 for the course CEE385-Performance Based Earthquake Engineering (as its TA), thought by Professor Miranda at Stanford. Since then, every year I've implemented some minor modifications. If you have any recommendation for new features (especially a new name, please... hehe), don't hesitate in sharing them with us. 

I'm currently working on a Python version of this tool as well.

This tool is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
