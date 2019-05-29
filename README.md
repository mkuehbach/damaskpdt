# damaskpdt
Additional MPI/OpenMP parallelized post-processing functionalities for DAMASK (the Düsseldorf Advanced MAterial Simulation Kit) to characterize the spatial distribution of state variable values as a function of the distance of each material point to interfaces.

# Dependencies:
C/C++ compiler, C2011 spec capable, e.g. GNU or Intel compiler  
Message Passing Library (MPI), MPIThreadFunneled capable, e.g. Intel MPI, MPICH  
Intel Math Kernel Library v2018.4  
Boost C/C++ header libraries v1.66  
Hierarchical data format library HDF5 v1.10.2  

# Prepare code compilation:
1.) Check existence and correct installation of all dependencies  
2.) Establish your shell environment variable (e.g. via module environment)  
3.) Modify upper section of the CMakeLists.txt to change path of your specific library locations.  
4.) Go into the damaskpdt build folder  

# Compile the code:
cmake -DCMAKE_BUILD_TYPE=Release ..  
make  

# Prepare an analysis:
0.) damaskpdt uses relative paths.  
1.) Make sure a spectralOut file and associated meta files *.outputConstitutive, *.outputHomogenization, *.outputCrystallite exist.  
2.) Set quantities and analysis tasks in the *.xml settings files. An example is in the build folder.  
3.) Set your environment for OpenMP multithreading:  
export OMP_NUM_THREADS=x  
export MKL_NUM_THREADS=1  
export OMP_PLACES=cores  
x should be the number of threads you want to use.  
Do not use two threads per hyperthreading core. It makes no sense in most cases.
Feel free to contact me for advise.  
Linear algebra operations on material points are threaded with OpenMP already.  
Therefore, do not use more than one thread for the Intel MKL library, as otherwise too many threads/tasks created.  

# Run an analysis:
export simid=y  
damaskpdt $simid Settings.xml MySpectralOut 1>DAMASKPDT.SimID.$simid.STDOUT.txt 2>DAMASKPDT.SimID.$simid.STDERR.txt  

y is here an unsigned integer (32bit) to distinguish results from different runs.  
Mind that multiple runs with the same ID and in the same folder will be overwritten silently!  
Settings.xml is the name of your specific settings file.  
MySpectralOut is the name of your spectralOut file.  
Do not add the extension *.spectralOut file; it will be added automatically!  

# Problems, questions:
Feel free to contact me: https://www.bigmax.mpg.de/39151/bigmax-software-engineering-consultant  

# Example dataset:
An example dataset is available on Zenodo with the following DOI: 10.5281/zenodo.1249282  




