A repository of tools to perform analysis of (e, e'pi) events using CLAS12 data.  To use, source:
```
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module load clas12
```
To compile:
```
mkdir build
cd build
cmake ../
make
```

TO DO:
- Update paths to data directory in classes
- Figure out cut value storage scheme (update Auxiliary class to make cleaner)
- Reimpliment various analysis macros from old analysis
