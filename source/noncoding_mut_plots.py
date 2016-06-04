#!python

import cPickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from operator import itemgetter

## Titles for plots
titles = {'mse': 'Mean Squared Error (MSE)',
					'r2': 'R-square Score (R2)',
					'muae': 'Mean Absolute Error',
					'mdae': 'Median Absolute Error',
					'expvar': 'Explained Variance Score',
					'imp': 'Feature Importance'}

## Marker for plots
markers = {'rf': 'go--', 
					 'sgd': 'ro--', 
					 'nb': 'bd--', 
					 'ada': 'ms--', 
					 'gradient': 'yv--',
					 'poisson': 'c*--',
					 'ann':'kp--',
					 'sae':'c+--'}

## Labels for plots
labels = {'rf': 'Random Forest', 
	        'ada': 'AdaBoosting',
	        'sgd': 'Stochastic Gradient Descent', 
	        #'svr' : 'SVR',
	        'gradient': 'Gradient Boosting',
	        'poisson': 'Poisson Regression',
	        'nb': 'Negative Binomial',
	        'ann': 'Neural Network',
	        'sae': 'Stacked AutoEncoder'}

features = ['tfbs_cnt',
						'tfbs_max_sc',
						'tfbs_avg_sc',
						'atf3_cnt',
						'atf3_max_sc',
						'atf3_avg_sc',
						'cebpb_cnt',
						'cebpb_max_sc',
						'cebpb_avg_sc',
						'cebpd_cnt',
						'cebpd_max_sc',
						'cebpd_avg_sc',
						'creb1_cnt',
						'creb1_max_sc',
						'creb1_avg_sc',
						'egr1_cnt',
						'egr1_max_sc',
						'egr1_avg_sc',
						'ets1_cnt',
						'ets1_max_sc',
						'ets1_avg_sc',
						'maff_cnt',
						'maff_max_sc',
						'maff_avg_sc',
						'dhs_src_cnt',
						'dhs_max_sc',
						'gerp_sc',
						'tss_dist',
						'gc_per']

#----------------------------------------------------------
# Put all metrics value in outoput of regr model into mem
#----------------------------------------------------------
def appendToMaps(predMaps,regressor,metricMaps,metrics):
	trainMetrics = metricMaps['train']
	testMetrics = metricMaps['test']

	for metric in metrics:
		if trainMetrics.has_key(metric):
			if metric != 'imp':
				(predMaps[metric])[regressor] = [trainMetrics[metric],testMetrics[metric]]
			else:
				(predMaps[metric])[regressor] = trainMetrics[metric]
	#//end for
	return predMaps
#//end appendToMaps

#------------------------------------------------------------
# Plot each metrics over differnt parameters in all models
#------------------------------------------------------------
def myplots(title,metric,markers,labels,regressors,params,vvMets,figname):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.title(title)

	i = 0
	first = 0
	indx = []
	if metric == 'imp':
		for vmet in vvMets:
			if len(vmet) > 0:
				# First we need to extract the max importance value of each feature
				# over different model parameters
				if regressors[i] not in ['poisson','nb']:
					vmet_max = np.amax(np.array(vmet),axis=0)
					if sum(vmet_max) != 0:
						#vmet_ranks = np.argsort(np.argsort(np.array(vmet),axis=1),axis=1)
						#max_vmet_ranks = np.amax(vmet_ranks,axis=0)
						#plt.plot(max_vmet_ranks,markers[regressors[i]],label=labels[regressors[i]])

						if first == 0:
							indx = vmet_max.argsort()[::-1]
							#features_new = itemgetter(*indx)(features)
							vmet_max = itemgetter(*indx)(vmet_max)
							first = 1
						else:
							vmet_max = itemgetter(*indx)(vmet_max)

						plt.plot(vmet_max,markers[regressors[i]],label=labels[regressors[i]])
				else:
					vmet_percent = vmet/float(sum(vmet))
					#print len(vmet),len(vmet_percent),len(indx)
					vmet_percent = itemgetter(*indx)(vmet_percent)
					plt.plot(vmet_percent,markers[regressors[i]],label=labels[regressors[i]])
			else:
				i = i + 1
				continue
			i = i + 1
	else:
		# plot other metrics: mse,r2
		for vmet in vvMets:
			plt.plot(vmet,markers[regressors[i]], label = labels[regressors[i]])
			i = i + 1
	#//end if

	if metric != 'imp':
		plt.ylabel(metric)
		plt.xlabel("Parameters-K")
		plt.xticks(np.arange(0,len(params['rf']),1.0))
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	else:
		## feature importances
		plt.xlabel('Features')
		plt.ylabel('Importance Score')
		#plt.xticks(range(len(features)),features,rotation=45,fontsize=6)
		plt.xticks(range(len(features)),itemgetter(*indx)(features),rotation=45,fontsize=6)
		plt.legend(loc='best')

	## Save figure to image file
	plt.savefig(figname,dpi=300,bbox_inches='tight')
	plt.show()
#//end function myplots

#-----------------------------------------------------------------------
# Plot all metrics value over different regressors and 
# different tunning parameters
#-----------------------------------------------------------------------
def plotMetrics(predMaps, regressors, params, metrics, outpath, msepath, out_suffix):	
	v_test_metrics = {}
	for metric in metrics:
		v_train_metric = []
		valid_regressors = []
		for regressor in regressors:
			if predMaps[metric].has_key(regressor):
				valid_regressors.append(regressor)
				if metric != 'imp':
					v_train_metric.append((predMaps[metric])[regressor][0])
					if not v_test_metrics.has_key(metric):
						v_test_metrics[metric] = [(predMaps[metric])[regressor][1]]
					else:
						v_test_metrics[metric].append((predMaps[metric])[regressor][1])
				else:
					v_train_metric.append((predMaps[metric])[regressor])
				#//end if
		#//end for

		title = titles[metric]
		figname = outpath+'nc.plots.train.sampled.'+out_suffix+"."+metric+'.png'
		if metric != 'imp':
			title = 'Validation-'+titles[metric]

		myplots(title,metric,markers,labels,valid_regressors,params,v_train_metric,figname)
	#//end for

	fignm = msepath+'nc.plots.test.sampled.'+out_suffix+'.png'
	outtsv= msepath+'nc.mse.out.test.sampled.'+out_suffix
	outfhdl=open(outtsv,'w')
	#myplots('Testing MSE','mse',markers,labels,regressors,params,v_test_metrics,fignm)
	print "regressors:",regressors
	outfhdl.write('\t'.join(regressors))
	outfhdl.write('\n')

	for metric in ['mse','muae','mdae','r2','expvar']:
		print metric,":",v_test_metrics[metric]
		if metric == 'mse':
			outfhdl.write('\t'.join(str(i) for i in v_test_metrics[metric]))
			outfhdl.write('\n')

	outfhdl.close()
#//end func plotMetrics


#----------------------------------------------------------
# Main function
#----------------------------------------------------------
def main():
	# Names format of the files generated by model learning
	tsv_suffix = "0.6"
	png_suffix = "0.6.1"

	# nc.mse.out.313.<regressor>.<param>
	#nm_prefix = 'nc.mse.out.316.1M.0.1M'
	nm_prefix = 'nc.mse.out.sampled.'+png_suffix
	
	## Output directory
	in_path = "/home/tanjinxu/Project/noncoding/sklearn_out/dumps/p60/"
	out_path = "/home/tanjinxu/Project/noncoding/sklearn_out/dumps/p60/figures/"
	mse_path = "/home/tanjinxu/Project/noncoding/sklearn_out/figures/"

	#metrics = ['mse','r2','ps','sp', 'imp', 'muae', 'mdae', 'expvar']
	metrics = ['mse','muae','mdae','r2','expvar','imp']

	params = {'rf': [30,40,50,100,150,200,300],           #number of trees
						'ada': [50,100,150,200,250,300,350],        #number of decision stumps
						'knn': [5,7,9,10,15,20,30],                 #number of neighbours
						'sgd': [5, 10, 20, 50, 100, 150, 200],      #number of iterations
						'svr': [0.001,0.005,0.01,0.05,0.1,0.5,1.0], #penalty of the error
						'gradient': [50,100,120,150,200,250,300]    #number of boosting steps
					}

	#regressors = ['poisson','nb']
	regressors = ['rf','ada','gradient','sgd','poisson','nb','ann', 'sae']

	mseMap = {}   # Mean Squared Error
	r2Map  = {}   # R-Squared score
	psMap  = {}   # Pearson-coefficient
	spMap  = {}   # Spearman-coefficient
	muaeMap = {}  # Mean Absolute Error
	mdaeMap = {}  # Median Absolute Error
	expsMap = {}  # Explained Variance Score
	impMap  = {}  # feature-importance

	predMaps = {'mse':mseMap,
							'r2':r2Map,
							'muae':muaeMap,
							'mdae':mdaeMap,
							'expvar':expsMap,
							'imp':impMap}

	for regressor in regressors:
		fname = in_path + nm_prefix + "." + regressor
		if regressor == 'ann' or regressor == 'sae':
			annRes = cPickle.load(open(fname,'rb'))
			expvars = []
			if annRes.has_key('exps'):
				expvars = annRes['exps']
			elif annRes.has_key('expvar'):
				expvars = annRes['expvar']

			annparms = annRes['neurons']
			train_metrics = {'mse': annRes['val_loss']}
			test_metrics = {'mse': annRes['mse'],
											'r2': annRes['r2'],
											'muae':annRes['muae'],
											'mdae':annRes['mdae'],
											'expvar': expvars}
			metricsMap = {'train':train_metrics, 'test':test_metrics}
		elif regressor == 'poisson' or regressor == 'nb':
			csv_error = in_path+"nc.res.raw."+regressor+".err."+tsv_suffix+".tsv"
			csv_imp = in_path+"nc.res.raw."+regressor+".imp."+tsv_suffix+".tsv"
			df_error = pd.read_csv(csv_error,
														 header = 0,
														 sep = '\t',
														 index_col = 0)
			df_imp = pd.read_csv(csv_imp,
														header = 0,
														index_col = 0,
														sep = '\t'
													)
			imp_array = df_imp.as_matrix()
			error_array = df_error.as_matrix()
			#print imp_array.shape,imp_array
			#print error_array
			
			train_metrics = {'muae':[error_array[0,0]]*7,
											 'mse': [error_array[1,0]]*7,
											 'meae': [error_array[2,0]]*7,
											 'expvar': [error_array[3,0]]*7,
									 		 'r2': [error_array[4,0]]*7,
									 		 'imp': imp_array[:,0]}
			test_metrics = {'muae': error_array[10,0],
											'mse': error_array[11,0],
											'meae': error_array[12,0],
											'expvar': error_array[13,0],
											'r2': error_array[14,0]}

			metricsMap = {'train':train_metrics, 'test':test_metrics}
		else:
			metricsMap = cPickle.load(open(fname,'rb'))
		#//end if

		predMaps = appendToMaps(predMaps,regressor,metricsMap,metrics)
	#//end for

	#print predMaps['imp']
	plotMetrics(predMaps,regressors,params,metrics,out_path,mse_path,png_suffix)
#//end main function

if __name__ == "__main__":
	main()
