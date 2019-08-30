Likelihood-free cosmological inference: multivariate response. (Section 4.2)
===

For multiple response components we want to model the often complicated dependencies 
between these components; 
this is in contrast to approaches which model each component separately, 
implicitly introducing an assumption of conditional independence. 
We use an example of LFI for simulated weak lensing shear data to show how NNKCDE 
and RFCDE can capture more challenging bivariate distributions with curved structures; 
in this toy example y represents cosmological parameters (\Omega_M, \sigma_8) 
in the LambdaCDM model, and x represents (coarsely binned) weak lensing 
shear correlation functions.

<br>

## Reproducing Results

* Setup the current folder as working directory;
* Data are from the GALSIM toolkit (https://github.com/GalSim-developers/GalSim) - 
we simulate 25,000 points as total number of simulations available;
* Run `lensing_posterior.R`, setting the `data_file` variable to reflect location of the
data in use;
* Code to obtain the figure, CDE Loss and also HPD coverage histograms is included
in `lensing_posterior.R`.


