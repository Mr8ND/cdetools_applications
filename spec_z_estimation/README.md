Spec-z estimation: functional input (Section 4.3)
===

Standard prediction methods, such as random forests, do not typically 
fare well with functional covariates, simply treating them as unordered 
vectorial data and ignoring the functional structure.
However, often there are substantial benefits to explicitly taking advantage 
of that structure as in fRFCDE.  
We compare the performance of a vectorial implementation of RFCDE with fRFCDE 
and Flexcode-Series (Lee & Izbicki 2016) for a spectroscopic sample from the 
Sloan Digital Sky Survey (Alan et al. 2015), 
where the input x is a high-resolution spectrum of a galaxy, 
and the response the galaxy's redshift z.

<br>

## Reproducing Results

* Setup the current folder as working directory;
* We use the SDSS spectra, applying Richards et al. 2009 and Gaussian noise. See the paper 
for details about the data;
* Once the data are available, say in the `data` folder, run `spec_z_preprocess.R` to apply
Gaussian noise and split between training and testing;
* Run `lensing_posterior.R`, setting the `data_file` variable to reflect location of the
data in use;
* Code to obtain the figure, CDE Loss and also HPD coverage histograms is included
in `lensing_posterior.R`.


