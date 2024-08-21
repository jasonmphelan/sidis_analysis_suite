A repository of tools to perform analysis of (e, e'pi) events using CLAS12 data.  To use, source:
```
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module load clas12
module switch clas12root/1.8.4 clas12root/dev
```
To compile:

Set environment variables SIDIS_DATA_PATH and SIDIS_HIST_PATH (used as output directory) and run
(will update)

```
mkdir build
cd build
cmake ../
make
```

TO DO:
- Update final skimmer to do 3d matching
- Update skimmers to handle MC and generator
- Reimpliment various analysis/correction macros from old analysis
