This project includes the codes (C++) used for the simulation of Fig4--Fig9 
in the manuscript
"Epigenetic state inheritance drivers drug-tolerant persister-induced resistance 
in solid tumors: A stochastic agent-based model". 
The C++ codes for model simulation and data generation.

The C++ codes include:
BCTool.cpp, System.cpp, CStemCell.cpp, Cell.cpp,Random.cpp and the corresponding header files

In parfiles, par.dat includes the system parameter values, 
drugschedules.dat includes different drug adaministration schedules used in Fig4-Fig9
md-**.in is used as input files of simulation.


The program has been test in Window 10 and Linux (ubuntu 22.04.4)

Usage:

windows:
1. Use IDE such VS code  or Dev C++ to generate the execute file BCTool.exe
2. in the command windown, BCTool  parfiles/md-**.in  to generate the simulation data.
md-**.in denote the files in parfiles with different schedules, such as md-control.in, md-sdrug.in

windows:
1. g++ -o BCTool BCTool.cpp to generate the execute file BCTool.exe
2. ./BCTool  parfiles/md-***.in  to generate the simulation data.
md-**.in denote the files in parfiles with different schedules,