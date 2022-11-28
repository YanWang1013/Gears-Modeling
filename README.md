## GEARS

GEARS (Global parameter Estimation with Automated Regularization via Sampling) is a MATLAB toolbox for parameter estimation in nonlinear dynamic models composed of deterministic ordinary differential equations (ODEs). GEARS is based on the combination of three main strategies: 
(i) global optimization (to avoid convergence to local solutions), 
(ii) reduction of the search space (i.e. tighter bounds on parameters), 
(iii) regularized estimation, a strategy used to handle overfitting (i.e. fitting the noise rather than the signal). 
As a result, GEARS can allow the user to avoid both underfitting and overfitting problems while requiring minimum supervision from the user. These capabilities are especially useful when calibrating ODE-based models with highly nonlinear and flexible dynamics (e.g. models of biological oscillators).

## Features

<ul>
<li>Unsupervised parameter estimation </li>
<ul>
<li>Efficient Global optimisation</li>
<li>Automated regularisation tuning</li>
<li>Parameter bounding</li>
</ul>
<li>Post-fit analyses</li>
<ul>
<li>Normalised root mean square error (NRMSE)</li>
<li>R<sup>2</sup> test</li>
<li>chi<sup>2</sup> test</li>
<li>Parameter uncertainty (using Fisher information)</li>
<li>Correlation matrix (using Fisher information)</li>
<li>Active bounds</li>
</ul>
<li>Post-fit plotting</li>
<ul>
<li>Trajectories for Fitting and cross-validation</li>
<li>Parameter space samples</li>
<li>Visualisation of parameter bound reduction</li>
<li>Residuals</li>
<li>Predictions vs measurements</li>
<li>Trajectory uncertainty (via Fisher information)</li>
<li>Convergence curves</li>
</ul>
<li>Results reports</li>
<ul>
<li>html markup of results</li>
<li>xls markup of results</li>
<li>Combined pdf report of all plotted Figures</li>
</ul>
</ul>

## Requirements

<p> GEARS runs on Matlab R2015b or later and is multi-platform (Windows and Linux). Both the optimisation and symbolic mathematics Matlab toolboxes are required to run GEARS. </p>

<p>GEARS requires that the <a href="http://icb-dcm.github.io/AMICI/">AMICI package</a> has been correctly installed.</p>

<p>Optionally, users can use <a href="https://www.ghostscript.com">Ghostscript</a> for the exportation of figure reports.</p>

![gear](https://user-images.githubusercontent.com/118581887/204173806-4ae4850a-209d-4ec2-8225-60bc6792c381.png)




