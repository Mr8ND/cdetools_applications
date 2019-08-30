Photo-z estimation: univariate response with multivariate input (Section 4.1)
===

This is the standard prediction setting in which all four methods apply. 
In Section 4.1, we apply the methods to the Teddy photometric redshift data 
(retrieved at https://github.com/COINtoolbox/photoz_catalogues), and illustrate 
the need for loss functions to properly assess PDF estimates of redshift z 
given photometric colors x.

<br>

## Reproducing Results

* Setup the current folder as working directory;
* Download Teddy photometric redshift data A and B 
(https://github.com/COINtoolbox/photoz_catalogues) and place them in a folder, say `Teddy`;
* Run `process_data.R`, setting the `datadir` folder as the folder you have downloaded
the Teddy data to. At this point you have setup the training and test data in the `data`
folder;
* Each method is run with the `method_*.R` file, with the exception of DeepCDE, which is
run using Python and the Tensorflow implementation. Each file will output a `.hdf5`
containing:
    - `y_grid`: the response grid over which the CDE has been calculated;
    - `y_true`: the true values of the response;
    - `cde`: the predicted conditional density estimates.
    
* Run `results.R` to reproduce the figures in the paper, along with other interesting
visualizations and photo-z metrics from Dahlen et al. (2013).


