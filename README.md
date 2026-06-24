# Gums taxonomy egypt
Scripts used in the paper:
New insights from MALDI-TOF MS and multivariate data analysis on the botanical origin of polysaccharide-basedpaint binders in ancient Egypt 

Clara Granzotto(1), Amra Aksamija(1,2), Gerjen H. Tinnevelt(3), Viktoriia Turkina(3), Ken Sutherland(1) 

(1) Art institute of Chicago, Department of Conservation and Science, 11 S. Michigin Ave, Chicago, IL 60603, USA.
(2) Center for Scientific Studies in Arts, Northwestern University, Tech building, 2145 Sheridan Road, Evanston, IL 60208, USA.
(3) Radboud University, Institute for Molecules and Materials,(Analytical Chemistry & Chemometrics), P.O. Box 9010, 6500 GL, Nijmegen, The Netherlands


## Description
Matlab scripts. 
Imports Maldi spectra from .txt, performs pre-processing, peak picking, pick alignment, calculates similarity matrix and performs PCA.


## Usage
Run the Analysis.m file. The script requires blanks in the current folder and additional folders with spectra. Reads only .txt files. From line 74 the analysis of Acacia was performed and creates the figures of the publication.

### R version (no coding required)
An R port with a point-and-click Shiny app is available in the `r/` folder. See **[README_R.md](README_R.md)** for setup and usage:

```r
source("r/install_packages.R")  # once
source("r/launch_app.R")        # opens the analysis app
```

Or run from the command line: `Rscript r/run_analysis.R`

## Support
gerjen.tinnevelt@ru.nl

## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
Open to contribution and collobaration. 

## Authors and acknowledgment
G.H. Tinnevelt
C. Granzotto, A Aksamija,V. Turkina, K. Sutherland

## Toolboxes used:
Francesco Savorani (2024). icoshift - interval correlation optimized shifting (for MATLAB V. 2014b and above) (https://www.mathworks.com/matlabcentral/fileexchange/29359-icoshift-interval-correlation-optimized-shifting-for-matlab-v-2014b-and-above), MATLAB Central File Exchange. Retrieved March 22, 2024.
Jeff (2024). simbin(mat1,mat2,type,mask) (https://www.mathworks.com/matlabcentral/fileexchange/55190-simbin-mat1-mat2-type-mask), MATLAB Central File Exchange. Retrieved March 22, 2024.

## License
CC-BY 4.0

## Project status
Paper submission.