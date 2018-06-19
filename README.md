# Matteo's electrophysiology

Scripts and helper functions for my electrophysiology data analysis. Not guaranteed to work on your machine, or at all.

## Organization
Things are in general organized according to their use, e.g. 'oscillations' for things which deal with LFP oscillations, or 'RFs' for receptive field mapping, etc. Note that not all code is guaranteed to be functional or fully-formed, as some (several) were aborted early on. 

## Dependencies
Probably the entire [spikes repository](https://github.com/cortex-lab/spikes "to it's GitHub page"), or at least very many functions therein. In addition, some parts require these:

+ Statistics and Machine Learning Toolbox
+ [Chronux](http://chronux.org/)
+ [npy-matlab](https://github.com/kwikteam/npy-matlab) (mainly `readNPY` and `writeNPY`)
+ [Rigbox](https://github.com/cortex-lab/Rigbox.git) (specifically +dat and `loadVar`)
+ [alyx-matlab](https://github.com/cortex-lab/alyx-matlab.git) (specifically +alf)
+ Optimization toolbox (only for RF fitting)
+ Signal Processing toolbox (only for cerebellum code)
+ [Circular Statistics toolbox](https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-) (only for cerebellum code)

In addition, all of the code for loading data and such .

## Hints
Some functions in this repository which were useful when exploring the data, in descending order of utility:

+ **plotlinked** Interactive plot (or imagesc) of some function for the selected point in a given scatter plot. Uses the mouse or Data Cursor ![alt text](https://uk.mathworks.com/help/matlab/ref/datacursortool.png "Mathworks") to select a point on the scatter plot, and show data for that point. Used e.g. to see how RF relates to various measured features, or to confirm that autocorrelogram shape matches with measured quantities.

+ **scatmat** Scatter-plot matrix with histograms along the diagonal, with the option to plot density rather than points, and to specify axis labels. Useful for looking at all the ephys features.

+ **quickLoadData** Load structures from a specific dataset of the given 'db' structure. Used to e.g. load the ephys information for the probe of a certain neuron being analyzed.

+ **samelims** A script which changes the current axis limits to be the same for x and y; different from the command-line `axis square`, which only makes the scales the same. 

+ **viewTogether** Plot two signals on upper and lower subplots, with linked x axes.