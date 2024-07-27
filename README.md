A repository of tools to perform analysis of (e, e'pi) events using CLAS12 data.  To use, source:
```
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module load clas12
module switch clas12root/1.8.4 clas12root/dev
```
To compile:
```
mkdir build
cd build
cmake ../
make
```

TO DO:
- Update final skimmer to take multiple files via 'reader' class
- Reimpliment various analysis/correction macros from old analysis
