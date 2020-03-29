# CDE Methods Applications: Photo-z, LFI and Spec-Z Estimation

Repository with code for examples (Section 4) of the "Conditional Density Estimation 
Tools in Python and R with Applications to Photometric Redshifts and Likelihood-Free 
Cosmological Inference" paper. <br>
Folders are arranged in the same way as the example section in the paper:
 - `photometric_redshift_teddy` is the code for photo-z estimation with 
 [TEDDY](https://github.com/COINtoolbox/photoz_catalogues) catalogue A and B (section 4.1);
 - `lfi_cosmological_inference` is the code for LFI for \Omega_M and \sigma_8 parameters
 examples using [Galsim](https://github.com/GalSim-developers/GalSim) toolkit (section 4.2);
 - `spec_z_estimation` is the code for the perturbed [SDSS 6](http://classic.sdss.org/dr6/) spectra and 
 spec-z classification (section 4.3).
 
To **reproduce the results** please see the README in the respective folder.

<br>

The methods can be found at the following repositories:
- [NNKCDE](https://github.com/tpospisi/nnkcde)
- [RFCDE](https://github.com/tpospisi/RFCDE)
- [Flexcode (R)](https://github.com/rizbicki/FlexCoDE) and [Flexcode (Python)](https://github.com/tpospisi/FlexCode)
- [DeepCDE](https://github.com/tpospisi/DeepCDE)
- [cdetools](https://github.com/tpospisi/cdetools)


Citation
===
@article{dalmasso2020cdetools,
       author = {{Dalmasso}, N. and {Pospisil}, T. and {Lee}, A.~B. and {Izbicki}, R. and
         {Freeman}, P.~E. and {Malz}, A.~I.},
        title = "{Conditional density estimation tools in python and R with applications to photometric redshifts and likelihood-free cosmological inference}",
      journal = {Astronomy and Computing},
         year = 2020,
        month = jan,
       volume = {30},
          eid = {100362},
        pages = {100362},
          doi = {10.1016/j.ascom.2019.100362}
}