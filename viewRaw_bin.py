import sys
from PyQt5 import QtGui
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QApplication
import pyqtgraph as pg
import numpy as np
import h5py as h5
import scipy.signal as sig
import scipy.io as scio

class currentData:
	""" The object which stores the data currently in view, and contains methods for loading and manipulating it.
	Needs to be handled from within the 'mainWindow' class."""
	
	def __init__(self, fname, intType = 'int16', edge_buffer = 5000, ref_chan = None):
		""" Initialise the current data according to default parameters. First reads the meta-data associated
		with the *.bin loaded from the GUI ('fname'), assuming it came from spike GLX. Stores data as nChan x nSamp arrays. 
		Also stores a 'time' array of the same size."""
		self.fname = fname
		self.intType = intType 
		self.e_buff = edge_buffer # don't know why this exists
		self.ref_chan = None
		
		print('\nLoading ... ') # first read metadata
		
		if 'CAR' in fname:
			meta = fname[:-8] + '.meta'
		else:
			meta = fname[:-4] + '.meta'
			
		nChans, fs, prb, filetime, filesize = self.readmeta(binname = self.fname, metaname = meta)
		if  prb == 'imec': # default parameters for probe types
			winTime = 0.2
			chan2plot = int(nChans)
			plot_chan = np.arange(chan2plot)
			scl = 38.5/(200*(2**15))
			space = 0.001
		elif prb == 'nidq':
			winTime = 1
			chan2plot = int(nChans) - 2
			plot_chan = np.arange(chan2plot)
			scl = 2.5/(200*(2**15))
			space = 0.0007
		else:
			sys.exit('\num?')
		
		self.winTime = winTime # lots of info about what we're lookin at
		self.plot_chan = plot_chan
		self.nChans = int(nChans)
		self.prb = prb
		self.fs = fs
		self.nSamp = int(self.fs*self.winTime)
		self.space = space
		self.scl = scl
		
		phil = open(self.fname,'rb')
		phil.seek(0,0)
		datatype = np.dtype(str(self.nChans)+'<i2')
		tempdat = np.fromfile(phil,dtype = datatype,count = self.nSamp)
		b1 = phil.tell()
		bytesize = (b1)/self.nSamp # the size in bytes of data size nSamp
		phil.close()
		
		if self.prb == 'imec': # if it's imec data, load the LFP separately
			# lfname = fname[:-6] + '.lf.bin'
			# templeaf = self.getleaf(lfname)
			print('ooo, fancy')
		else: # otherwise, highpass filter the data
			now = pg.ptime.time()
			print('Filtering ... ') # can take some time depending on the data
			filtdat = np.array(self.getfilts(np.transpose(tempdat)))
			print('... took %0.2f seconds' % (pg.ptime.time() - now))
			self.filt = filtdat
		
		t = np.linspace(0,self.winTime,self.nSamp)
		
		self.time = np.tile(t,(chan2plot,1))
		self.raw = np.transpose(tempdat) # the raw data
		self.byteSize = bytesize
		self.fileTime = filetime
		self.fileSize = filesize
		self.fileSamps = int(filetime*fs)
		self.moved = True
		self.where = 0
		
	def updateData(self,newpos = None):	
		phil = open(self.fname,'rb')
		datatype = np.dtype(str(self.nChans)+'<i2')
		
		if newpos == None:
			self.moved = False
			newpos = self.where
		else:
			self.moved = True
			newtime = newpos / self.fs
			t = np.linspace(newtime,newtime+self.winTime,self.nSamp)
			self.time = np.tile(t,(len(self.plot_chan),1))
	
		newdat = self.loadbin(self.fname, pos = newpos, leap = self.nSamp, bytesize = self.byteSize)
		if self.prb == 'imec':
			ind = self.fname.index('ap')
			lfname = self.fname[:ind] + 'lf.bin'
			lfp = self.getleaf(lfname)
	
			self.raw = np.transpose(newdat) + lfp
			self.filt = np.transpose(newdat)
		else:
			newfilt = np.array(self.getfilts(np.transpose(newdat)))
			self.raw = np.transpose(newdat)
			self.filt = newfilt
			
		self.where = newpos
		
	def readmeta(self, binname, metaname = None):
		""" Given a *.bin file name, will read the metadata from the corresponding *.meta structure.
		Assumes that this structure is coming from SpikeGLX."""
		
		if metaname == None:
			metaName = binname[:-4] + '.meta'
		else:
			metaName = metaname
		with open(metaName,'r') as phil:
			for lin in phil: # add contitions here
				ind = lin.index('=') + 1 # assumes no space after the equals
				if 'nSavedChans' in lin:
					nChans = lin[ind:-1]
				elif 'SampRate' in lin:
					fs = int(lin[ind:-1])
				elif 'typeThis' in lin:
					prb = lin[ind:-1]
				elif 'fileTimeSecs' in lin:
					fileTime = int(float(lin[ind:-1]))
				elif 'fileSizeBytes' in lin:
					fileSize = int(float(lin[ind:-1]))
					
		return nChans, fs, prb, fileTime, fileSize
	
	def loadbin(self,binname, pos, leap, bytesize = None):
		"""Load binary file. Must give the filename, desired position, and how many samples to take. Specify the size of each
		sample in bytes ('bytesize'), or it will calculate from the file size and number of samples."""
		datatype = np.dtype(str(self.nChans)+'<i2')
		if bytesize == None:
			bs = int(self.fileTime*self.fs)/self.fileSize
			bs = int(bs)
		else:
			bs = bytesize
		
		phil = open(binname,'rb')
		bytepos = int(phil.tell() + pos*bs)
		phil.seek(bytepos,0)
		dat = np.fromfile(phil, dtype = datatype, count = leap)

		phil.close()
		
		return dat
		
	def getleaf(self, lfname):
		""" Given the name of an *.imec.lf.bin file, will load it, up-sample to Fs of ap band, and store it as 'leaf'.
		Position ('pos') should be given in units of raw data samples. 
		Asumes that the number of samples to take is already specified and stored as an attribute of self."""
		_, lffs, _, _, lfsize = self.readmeta(lfname)
		leafbyte = int(self.fileTime * lffs)/lfsize
		leafsamp = int(self.nSamp * (lffs/self.fs))
		rawlfp = self.loadbin(lfname, pos = self.where, leap = leafsamp, bytesize = leafbyte)
		rawlfp = np.transpose(rawlfp)
		
		x = np.arange(self.nSamp) # see SciPy documentation
		xp = np.arange(0,self.nSamp,int(self.fs/lffs))
		leaf = []
		for i, chan in enumerate(rawlfp):
			leaf.append(np.interp(x,xp,rawlfp[i]))
		
		leaf = np.array(leaf)		
		return leaf	

	def getfilts(self, rawdat, fstop = 600, fpass = 700):
		fnyq = self.fs / 2.
		desired = (0, 0, 1, 1)
		bands = (0, fstop, fpass, fnyq) # see SciPy documentation 
		filtcoefs = sig.firls(1001, bands, desired, nyq=fnyq)
		filtdat = sig.filtfilt(filtcoefs,[1],rawdat,axis = -1)
		return filtdat
		
	# def findSpikes(self): # work in progress
		# # Uses the loaded clusters to make a logical indexing of spike times
		# # The vector is length nClu x nSamp, with a 30-sample window highlighting the time of each spike, 
		# # indexed according to the main channel of that cluster (i.e. [0 0 0 ... 0 3 3 3 3 ... 3 3 3 0 0 0 ...]).
		
		# nSpks = self.spikes.shape[0]
		# kern = np.ones((1,51))
		# kern[:21] = 0
		
		# current_window = self.spikes < self.time[-1] and self.spikes > self.time[0]
		# newspikes = np.zeros((nSpks, self.nSamp))
		# for iClu, cluster in enumerate(spikes):
			# newspikes[cluster[current_window]] = sites[iClu]
			
class mainWindow(QtGui.QMainWindow): 

	def __init__(self): 
		# Make the UI
		QtGui.QMainWindow.__init__(self)
		print(len(sys.argv))
		if 	len(sys.argv) > 1:
			if '.bin' in sys.argv[1]:
				raw_file = sys.argv[1]
		else:
			raw_file, _ = QtGui.QFileDialog.getOpenFileName(self, 'Open Raw') 
			if len(raw_file) == 0:
				sys.exit('\noh, okay, then')
	
		data = currentData(fname = raw_file)
		
		self.view = pg.GraphicsLayoutWidget()
		self.setCentralWidget(self.view)
		
		newFile = QtGui.QAction("&Load data", self) # open new bin from within the app
		newFile.setShortcut("Ctrl+O")
		newFile.setStatusTip('Open bin file')
		newFile.triggered.connect(self.new_bin)
		openSpikes = QtGui.QAction("&Load clusters", self) # so the I can see spikes
		openSpikes.setShortcut("Ctrl+Q")
		openSpikes.setStatusTip('Open File')
		openSpikes.triggered.connect(self.spike_open)
		getChanMap = QtGui.QAction("&Load channel map", self) # to add channel mapping
		getChanMap.setShortcut("Ctrl+M")
		getChanMap.setStatusTip('Open File')
		getChanMap.triggered.connect(self.cmap_open)
		
		changeWind = QtGui.QAction("&Change window", self) # for changing the time to plot
		changeWind.setShortcut("Shift+T")
		changeWind.setStatusTip('Enter value')
		changeWind.triggered.connect(self.new_winsize)
		changeChan = QtGui.QAction("&Change channels", self) # for changing the channels to plot
		changeChan.setShortcut("Shift+Q")
		changeChan.setStatusTip('Enter value')
		changeChan.triggered.connect(self.new_plotchan)
		
		mainMenu = self.menuBar()
		
		fileMenu = mainMenu.addMenu('&File')
		fileMenu.addAction(newFile)
		fileMenu.addAction(openSpikes)
		fileMenu.addAction(getChanMap)
		
		windMenu = mainMenu.addMenu('&Window')
		windMenu.addAction(changeWind)
		windMenu.addAction(changeChan)
        
		# Plot 
		self.plotWindow = self.view.addPlot()
		self.plotWindow.disableAutoRange()
		
		self.ds = 2
		self.data = data
		self.whichdata = 'raw'
		self.istart = 0
		self.ichan = 0
		
		self.updatePlot()
		print('\nShowing: \n%s' % raw_file)

	def updatePlot(self):
		if self.whichdata == 'raw':
			dat = self.data.raw[self.data.plot_chan,:self.data.nSamp]
		elif self.whichdata == 'filt':
			dat = self.data.filt[self.data.plot_chan,:self.data.nSamp]
		
		dat = np.multiply(dat,self.data.scl,casting = 'unsafe')
		dat = np.add(dat,(np.arange(len(self.data.plot_chan))[:,np.newaxis]*self.data.space),casting = 'unsafe')
		
		self.plotWindow.disableAutoRange()
		self.plotWindow.clear()
		# now = pg.ptime.time()
		
		lines = MultiLine(self.data.time, dat, ds = self.ds)
		self.plotWindow.addItem(lines) # plot it
		for onee,chan in enumerate(dat): # add channel labels
			labels = pg.TextItem(str(self.data.plot_chan[onee]), 'w',anchor = (1,1), border = pg.mkPen(0.5),fill = pg.mkBrush(0.0))
			labels.setPos(float(self.data.time[0,0]),float(chan[0]))
			self.plotWindow.addItem(labels, ignoreBounds = True)
			
		# print('Updated in %0.2f seconds' % (pg.ptime.time() - now))
		if self.data.moved:
			self.plotWindow.autoRange() 
	
	def keyPressEvent(self, e):
		if e.key() == Qt.Key_W:
			e.ignore()
			
		elif e.key() == Qt.Key_S:
			if (self.istart + self.data.plot_chan) < self.data.nChans:
				self.istart += self.data.winTime
			else:
				self.istart = self.data.nChans - self.data.nSamp
			self.data.updateData(self.istart)
			self.updatePlot()

		elif e.key() == Qt.Key_A:
			if self.istart >= self.data.nSamp:
				self.istart -= self.data.nSamp
			else:
				self.istart = 0
			self.data.updateData(self.istart)
			self.updatePlot()
			
		elif e.key() == Qt.Key_D:
			if (self.istart + self.data.nSamp) < self.data.fileSamps:
				self.istart += self.data.nSamp
			else:
				self.istart = self.data.fileSamps - self.data.nSamp
			self.data.updateData(self.istart)
			self.updatePlot()
			
		elif e.key() == Qt.Key_F:
			if self.whichdata == 'raw':
				#if self.data.prb != 'imec':
				self.whichdata = 'filt'
				#else:
				#	print('\nHey, don\'t do that.')
			elif self.whichdata == 'filt':
				self.whichdata = 'raw'
			self.data.moved = False
			self.updatePlot()
		
		elif e.key() == Qt.Key_Z: # get rid of downsampling
			self.ds = 1
			self.data.moved = False
			self.updatePlot()
			
		elif e.key() == Qt.Key_Left: # increase the downsampling
			if self.ds > 1:
				self.ds -= 1
			else:
				e.ignore()
			self.data.moved = False
			self.updatePlot()
			
		elif e.key() == Qt.Key_Right:
			self.ds += 1
			self.data.moved = False
			self.updatePlot()
		
		elif e.key() == Qt.Key_Up: #change scale
			if e.modifiers() & Qt.ShiftModifier:
				self.data.space *= 1.25 # change spacing if shift + pressing
			else:
				self.data.scl *= 1.5
			self.data.moved = False
			self.updatePlot()
			
		elif e.key() == Qt.Key_Down:
			if e.modifiers() & Qt.ShiftModifier:
				self.data.space *= 0.75
			else:
				self.data.scl *= 0.666
			self.data.moved = False
			self.updatePlot()
			
		elif e.key() == Qt.Key_Q:
			if e.modifiers() & Qt.ShiftModifier:
				self.new_plotchan()
			
		else:
			e.ignore()
	
	def new_bin(self):
		raw_file, _ = QtGui.QFileDialog.getOpenFileName(self, 'Open Raw') 
		if len(raw_file)>0:
			data = currentData(fname = raw_file)
			self.data = data
			self.updatePlot()
			print('\nShowing: \n%s' % raw_file)
		else:
			print('oh, okay, then')
	
	def new_winsize(self):
		num, ok = QtGui.QInputDialog.getDouble(self,"Select window size","enter a number (in seconds)")
		if ok:
			self.data.nSamp = int(self.data.fs*num) 
			self.data.winTime = num
			self.data.updateData()
			self.updatePlot()
			self.plotWindow.autoRange()
	
	def new_plotchan(self):
		"""Convert string input to numpy array."""
		input_chans, ok = QtGui.QInputDialog.getText(self,"Specify channels","enter new channels (e.g: 1, 3:6, 10)")
		if ok:
			chans = [] # ugly it but works
			if ',' in input_chans:
				input_chans = input_chans.split(',')
			else:
				input_chans = input_chans.split(' ')
				
			for val in input_chans: 
				if ':' in val:
					c = val.index(':')
					if len(val) == 1:
						chans = range(0, self.data.nChans)
						break
					else:
						first = int(val[:c])
						last = int(val[c+1:]) + 1
						[chans.append(i) for i in range(first,last)]
				else:
					numb = int(val)
					chans.append(numb)
			
			chans = np.array(chans)
			self.data.plot_chan = chans.flatten()
			self.data.updateData()
			self.updatePlot()
			self.plotWindow.autoRange()
	
	def spike_open(self): # work in progress
		"""Open a *_jrc.mat file and extract the spike times, cluster identity, and main channel.
		Adds these to the current data structure. Check this in case the structure names change."""
		
		name, _ = QtGui.QFileDialog.getOpenFileName(self, 'Open File') # this loads a tuple normally
		if len(name) > 0:
			spks_mat = scio.loadmat(name)
			try:
				spikes = spks_mat['t_spk'][0]
				sites  = spks_mat['t_spk'][1]
			except:
				print('\nHey!')
				print('\n "' + name + '" isn\'t a cluster file!')
				return
			
			self.data.spikes = spikes
			self.data.sites = sites
			
	def cmap_open(self): # work in progress
		"""Open a *.prb or *.mat file and add the channel map to current data structure."""
		# name = QtGui.QFileDialog.getOpenFileName(self, 'Open File')
		
		# if len(name) > 2:
			# with open(name,'r') as phil:
				# for lin in phil:
					# if 'channels' in lin:
						# cmap = []
						# if '...' in lin:
							# ind = lin.index('[') + 1
							# fin = lin.index('...')
							# vect = lin[ind:fin]
							# cmap.append(int(i) for i in vect)
						# elif ';' in lin:
							# ind = lin.index('[') + 1
							# fin = lin.index(']')
							# vect = lin[ind:fin]
							# cmap.append(int(i) for i in vect)
							# break
		
		# else:
		
		# self.cmap = cmap
			
class MultiLine(pg.QtGui.QGraphicsPathItem):
	# Skeleton courtesy of Luke Campagnola
	
	def __init__(self, x, y, slice = 'all', ds = 1, col = 'w'):
		"""x and y are 2D arrays of shape (Nplots, Nsamples).
		Slice is a logical vector saying which values to plot. If given, must be a vector of length nSamp.
		By default, it is 'all', and all samples are plotted.
		'ds' is the downsampling factor. At the moment, it just subsamples, but could be made to do something
		better, like take the mean/max in each bin."""
		
		if slice == 'all':
			connect = np.ones(x.shape, dtype=bool)
			connect = connect[:,::ds] # downsample
			connect[:,-1] = 0 # don't draw the segment between each trace
		else:
			connect = np.zeros(x.shape, dtype=bool)
			connect[:,slice] = 1 # draw segments only between specified points
			connect = connect[:,::ds] # downsample
			connect[:,-1] = 0 # just to make sure
			
		x = x[:,::ds]
		y = y[:,::ds]
		self.path = pg.arrayToQPath(x.flatten(), y.flatten(), connect.flatten())
		pg.QtGui.QGraphicsPathItem.__init__(self, self.path)
		self.setPen(pg.mkPen(col))
		
	def shape(self): # override because QGraphicsPathItem.shape is too expensive.
		return pg.QtGui.QGraphicsItem.shape(self)
	
	def boundingRect(self):
		return self.path.boundingRect()

		
def main(): 
	app = QtGui.QApplication(sys.argv) 
	
	win = mainWindow()
	win.show()
	
	sys.exit(app.exec_())
	
main()