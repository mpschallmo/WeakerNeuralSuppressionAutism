# WeakerNeuralSuppressionAutism
Computational modeling code (MATLAB) for the paper "Weaker Neural Suppression in Autism" by Schallmo and colleagues,
available on bioRxiv: doi.org/10.1101/645846

1. System requirements
- Requires MATLAB (The MathWorks), and the Statistics and Machine Learning Toolbox
- Also requires makeNeuralImage.m & convolveImage.m, which are included in the directory
- Tested on version R2016a

2. Installation guide
- Clone the repository
- If you wish to run the code outside of the directory, add the directory to your
  MATLAB paths
- Time to install: < 1 min

3. Instructions for use
- To run, call the function "Normalization_Models" at the command line
- Inputs: which_version - integer, acceptable values are:
    1 - Weaker normalization model, depicted in Figure 5A, B, & C
    2 - Larger excitatory spatial filter model, from Figure 5D, E, & F
    3 - Narrower spatial top-down modulation model, from Figure 5G, H, & I
    4 - Extra-narrow top-down modulation, from Supp. Figure 4A, B, & C
    5 - Ultra-narrow top-down modulation, from Supp. Figure 4D, E, & F
- Outputs: model - structure, with fields:
    p = model parameters
    stim = stimulus parameters
    resp = model peak response (arbitrary units), used to calculate 
           duration thresholds
    thresh = model duration threshold
    SI = model size indices
    variable_param = structure, with field = parameter that varies in the
                     chosen model version
    which_version = which model version (input, as above)
    dimension_key = key for dimension size of resp & thresh
- Time to run: < 1 min

This code is made available under the following license:
Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) 
https://creativecommons.org/licenses/by-nc/4.0/

We wish to thank Geoffrey M. Boynton for providing MATLAB functions for the model.
