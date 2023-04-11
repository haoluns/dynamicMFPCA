# MDFPCA (Dynamic Prediction for Sparse Time-Varying Images using Functional Principal Component Analysis)
R codes for implementing the dynamic Prediction for sparse time-varying images using functional principal component analysis.
# Description
Our work is motivated by predicting the progression of Alzheimer's disease (AD) based on a series of longitudinally observed brain scan images.  Existing works on dynamic prediction for AD focus primarily on extracting predictive information from multivariate time-varying biomarker values or brain imaging data at the baseline; whereas in practice, the subject's brain scan image represented by a multi-dimensional data matrix is collected at each follow-up visit. It is of great interest to predict the progression of AD directly from a series of longitudinally observed images. We propose a novel multi-dimensional functional principal component analysis based on alternating regression on tensor-product B-spline, which circumvents the computational difficulty of doing eigendecomposition, and offers the flexibility of accommodating sparsely and irregularly observed image series. We then use the functional principal component scores as features in the Cox proportional hazards model. We further develop a dynamic prediction framework to provide a personalized prediction that can be updated as new images are collected. Our method extracts visibly interpretable images of the functional principal components and offers an accurate prediction of the conversion to AD. We examine the effectiveness of our method via simulation studies and illustrate its application on the Alzheimer's Disease Neuroimaging Initiative data.


The repository includes two functions:
* ```func_mdfpca.r```: The R source code for fitting the FPCA for sparse time-varying images, which contains the following key functions.
```rscript
first_FPC_2d_image(beta1, observed, timepoints1, timepoints2, 
basis1, basis2, threshold, minit) 
second_FPC_conditional_2d_image(beta1, pc_index, observed, timepoints1, timepoints2, 
basis1, basis2, betalist, threshold, minit)
```
* ```demo_simu.r```: A demo script for a simulation study of the FPCA for sparse time-varying images.
* ```simu.RData```: The first three FPCs for the simulation study. Can be downloaded with ```download.file```.

