# Matteo's electrophysiology

Scripts and helper functions for my electrophysiology data analysis. Not guaranteed to work on your machine, or at all.

## Organization
Things are in general organized according to their use, e.g. 'oscillations' for things which deal with LFP oscillations, or 'RFs' for receptive field mapping, etc. Note that not all code is functional or fully-formed, as some (several) were aborted early on. Anything general purpose.

## Dependencies
Probably the entire [spikes repository](https://github.com/cortex-lab/spikes "to it's GitHub page"), or at least very many functions therein. In addition, some parts require these:

+ Statistics and Machine Learning toolbox
+ [npy-matlab](https://github.com/kwikteam/npy-matlab) (mainly `readNPY` and `writeNPY`)
+ [Rigbox](https://github.com/cortex-lab/Rigbox.git) (specifically +dat and `loadVar`)
+ [alyx-matlab](https://github.com/cortex-lab/alyx-matlab.git) (specifically +alf)
+ [Chronux](http://chronux.org/) (only for e-phys features)
+ Optimization toolbox (only for RF fitting)
+ Signal Processing toolbox (only for cerebellum code)
+ [Circular Statistics toolbox](https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-) (only for cerebellum code)

In addition, all of the code for loading data assumes the standard (as of June 2018) cortex lab heirarchy for data storage.

## Hints
If you haven't, check out [this handy guide](https://github.com/cortex-lab/neuropixels/wiki/Other_analysis_methods) for looking through Neuropixels data.

On top of that, here are some functions which were useful when exploring my data, in descending order of utility:

+ **plotlinked** Interactive plot (or imagesc) of some function for the selected point in a given scatter plot. Uses the mouse or Data Cursor ![alt text](https://uk.mathworks.com/help/matlab/ref/datacursortool.png "Mathworks") to select a point on the scatter plot, and show data for that point. Used e.g. to see how RF relates to various measured features, or to confirm that autocorrelogram shape matches with measured quantities.

+ **scatmat** Scatter-plot matrix with histograms along the diagonal, with the option to plot density rather than points, and to specify axis labels. Useful for looking at all the ephys features.

+ **quickLoadData** Load structures from a specific dataset of the given 'db' structure. Used to e.g. load the ephys information for the probe of a certain neuron being analyzed.

+ **samelims** A script which changes the current axis limits to be the same for x and y; different from the command-line `axis square`, which only makes the scales the same. 

+ **viewTogether** Plot two signals on upper and lower subplots, with linked x axes.

+ **viewRaw_bin.py** A super wobbly app for viewing raw neuropixels data in spikeGLX format. Made it before I realized there was a way to do this with phy. It's run on the commond line from the ephys directory by doing either:<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`python viewRaw_bin.py`<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;which will let you navigate to the *.ap.bin file, or<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`python viewRaw_bin.py Path\To\Your\File\foo.imec.ap.bin`<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbspwhich loads it directly. From there nagivate with A/D for back/forth in time, and up/down arrow keys for scale.<br/>


