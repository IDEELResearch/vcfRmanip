# vcfR2manip
  [![Travis build status](https://travis-ci.org/IDEELResearch/vcfR2manip.svg?branch=master)](https://travis-ci.org/IDEELResearch/vcfRmanip). 
 [![Codecov test coverage](https://codecov.io/gh/IDEELResearch/vcfR2manip/branch/master/graph/badge.svg)](https://codecov.io/gh/IDEELResearch/vcfRmanip?branch=master)

# Overview 
This package includes a few useful functions for manipulating VCFs. It is not meant to be comprehensive and at times, will not be memory efficient. Many of these functions are inspired by personal work (that I have found to be redundant between projects) and are meant to be "quick" and "small" manipulations.

This package will undergo version control for manuscript posterity. Major and minor version changes will be recorded.

### Notes

  1.  Some of these functions were inherited from `NFBtools` (now depreciated)
  2. This package has several redundancies with `vcfdo`. The key difference is that this package is meant to be interactive within the Rstudio interface versus the command line. However, the Rstudio environment and many of these functions are memory limited. This will be an issue with large vcf files, where `vcfdo` may be preferable. **Redundancies include**:   
  	1. LD pruning  
  	2. WSAF calculations (storage difference)   
  	3. `dab` calculation from PMC4786412 .  
