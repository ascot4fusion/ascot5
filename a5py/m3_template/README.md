#The directories content

## Micro

The 'micro' model - provides a template of the MPI Python code that can interact with MUSCLE3 library

## Macro

The 'macro' model - very simple logic, added to allow the micro model tests. 
The model launches the micro model several times, sending to it the input IDSes and receives the out IDS 
at every iteration

## Workflow
yMMLS file that defines micr-macro models connections and their implementations. 

**WARNING!**
Remember to correct the code parameter file path to point to a proper locations!
(The MUSCLE3 0.8.0 cannot use system variables in 'settings' part of the yMMSL file.)

## Scripts
`set-env.sh`, `run.sh`:  Auxiliary scripts to set up the environment (IMAS/3.42.0/foss) and run the workflow.
