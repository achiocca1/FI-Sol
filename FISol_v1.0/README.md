ALL THE FILES PRESENT IN THE FOLDER REFER TO A READY-TO-USE MODEL.

JUST TRY TO RUN MAIN_FI.m IN MATLAB ENVIRONMENT TO TEST THE CODE.

************************
# DESCRIPTION OF THE FILES 
************************

- **ANSYS_Post_Process_Results.mac** = APDL command necessary for a correct post-processing of results in Ansys (it can be used only in Ansys, if any other FEM software is used you need to correctly post-process results as the ones contained in the file RESULTS.csv and COORD.csv)

- **MAIN_FI.m** = main Matlab function to be run

- **RESULTS.csv** and **COORD.csv** = files containing all the nodal results of stress and strain from the finite element analysis

*********************************
# NECESSARY STEPS TO RUN THE SCRIPT 
*********************************
0) **Generate the stress and strain tensors result files**

>*Only in the case of Ansys Workbench and APDL you can use the ready-to-use code provided (ANSYS_Post_Process_Results.mac). If another FEM software is used you have to export stress and strain results by creating .csv files structured as the ones contained in **RESULTS.csv** and **COORD.csv** *
>- 0.1) Create a Named Selection called "Nodes_Circ" and select all the nodes that have to be evaluated through the critical plane method
>- 0.2) Paste and copy the ANSYS_Post_Process_Results.mac in an APDL command in the solution environment of Ansys Workbench. The code will automatically generate .csv files containing all the stress and strain results at each load step.

1) **Open MAIN_FI.m**
2) **Set in PARAMETERS the quantities you are interested in:**

- LoadSteps             % Considered loadsteps for Findley critical plane factor calculation
- kFin                  % Material constant of Findley critical plane factor
- R                     % Load ratio (it is necessary to be defined only when one load step is involved in the simulation)

3) **RUN the script**
4) **ENJOY!**
