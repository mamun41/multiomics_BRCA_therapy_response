from matplotlib import rcParams
import matplotlib.pyplot as plt
from aux_performance import join_outputs
from aux_multiseed_performance import plotMultiPCR
from aux_plotting import applyPlotStyle


rcut=1
applyPlotStyle()

from aux_multiseed_performance import getROC

train_parent_folder = '../trained_models/results_submission_20200916_100417'

testfiles = {}
testfiles['neg'] = '../validation_results/results_her2neg'
testfiles['pos'] = '../validation_results/results_her2pos'

output_folder = 'graphs-files'

getROC(testfiles,rcut,'agnost','pCR',['diablo427','chemo'],['DIABLO','ORIGINAL'],output_folder=output_folder)