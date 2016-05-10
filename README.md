# DTWave_cluster
Unsupervised clustering using artificial neural network and dynamic time warping


##### Parameters Description Default

'htkconfigdir', ''
'verbose', 0
'outputfile', fullfile(wavdir,'dtwave_cluster.csv'),
'wintime', 0.3
'hoptime', 0.1,...
'numcep', 13
'lifterexp', -22
'sumpower', 1
'preemph', 0.97
'dither', 0
'minfreq', 300
'maxfreq', 20000
'nbands', 26
'bwidth', 1.0
'dcttype', 3
'fbtype', 'htkmel'
'usecmp', 0
'modelorder', 0
'broaden', 0
'epoch', 4
'LR', [0.5 0.05]
'NS', [2 1]
'NC', [4 2]
'thexpan', -1
'gama', 0.95
'ce', 0
'nrepeat', 5
'depth', 0



## References:

A technical description of the method is available in:
Louis Ranjard, Howard A. Ross, Unsupervised bird song syllable classification using evolving neural networks, Journal of the Acoustical Society of America 123(6):4358-4368, 2008



--
Copyright (c) 2008-2016, Louis Ranjard
