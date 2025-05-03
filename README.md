# OrthogonalizedSyntheticControl

Contains R code to for Synthetic Control Estimators which estimates the weights on control units by using outcomes of units excluded from the set of controls as instruments.

GMM-SCE.R gives a function to estimate the control weights and model selection method for splitting units into instruments and controls based on "A Method of Moments Approach to Asymptotically Unbiased Synthetic Controls" by Fry (2024). 

OrthogonalizedSCE.R has function for an orthogonalized estimate of the average treatment effect on the treated unit.

RegularizedEstimate.R has functions to estimate the control weights and the weights on the moment conditions for OrthogonalizedSCE.R.

SeriesHAC.R has functions to estimate variance using Orthonormal Series HAC estimator for OrthogonalizedSCE.R.
