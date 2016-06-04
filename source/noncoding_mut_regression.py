#!python

import math
import numpy as np
import cPickle
import scipy
import datetime
import random
import gc
from time import time

from sklearn.grid_search import GridSearchCV
from sklearn import linear_model
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import RandomForestRegressor,AdaBoostRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.svm.classes import SVR
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import mean_squared_error,r2_score
from sklearn.metrics import mean_absolute_error, median_absolute_error
from sklearn.metrics import explained_variance_score
from sklearn.cross_validation import KFold

#-------------------------------------------------------------
# Split data into training set and testing/validation set or
# just sampling part of the data
#-------------------------------------------------------------
def sample_split_data(x,y,labels,split_rate):
	(nrows,ncols) = x.shape
	mcords = [str.split(cord,'-') for cord in labels]
	chr_cords_list = [cord[0]+str(int(cord[1])/10**6) for cord in mcords]

	uniq_list_val = list(set(chr_cords_list))
	uniq_first_val = random.sample(uniq_list_val,int(split_rate*len(uniq_list_val)))
	uniq_second_val = list(set(uniq_list_val)-set(uniq_first_val))

	first_indices = [i for i,item in enumerate(chr_cords_list) if item in uniq_first_val]
	second_indices = list(set(range(len(chr_cords_list))) - set(first_indices))

	first_x = x[first_indices,]
	first_y = y[first_indices,]
	second_x = x[second_indices,]
	second_y = y[second_indices,]
	return first_x,first_y,second_x,second_y
#//end sample_split_data

#-------------------------------------------------------------------
# K-Fold Cross Validation
#-------------------------------------------------------------------
def crossval_predict(estimator,k_fold,train_x,train_y,train_labels,has_importances):
	## Predicted feature importances in training
	train_avg_imp = np.zeros(train_x.shape[1])

	## validation loss and errors
	valid_avg_mse = 0.0
	valid_avg_r2 = 0.0
	valid_avg_sp = np.array([0,0])
	valid_avg_ps = np.array([0,0])
	valid_avg_muae = 0.0
	valid_avg_mdae = 0.0
	valid_avg_expvar = 0.0

	for i in range(k_fold):
		## Split the train data into train_subset and validation_subset
		## based on the hash-value of <chrom_id>_<mutation_coordinate>
		## Note: originally I was trying to just split randomly with the k-fold
		x_train,y_train,x_valid,y_valid = sample_split_data(train_x,train_y,train_labels,1-1/float(k_fold))

		x_train = np.array(x_train,dtype=np.float32)
		y_train = np.array(y_train,dtype=np.float32)
		x_valid = np.array(x_valid,dtype=np.float32)
		y_valid = np.array(y_valid,dtype=np.float32)

		print("train:{}, validation:{}".format(x_train.shape,x_valid.shape))

		## Fit the training data set
		reg_fit = estimator.fit(x_train,y_train)
		if has_importances:
			train_avg_imp = train_avg_imp + reg_fit.feature_importances_

		## Fit the model on the validation dataset
		y_hat_valid = reg_fit.predict(x_valid)

		mse_valid = mean_squared_error(y_valid,y_hat_valid,multioutput='uniform_average')
		r2_valid = r2_score(y_valid,y_hat_valid)
		pearson_valid = scipy.stats.pearsonr(y_valid,y_hat_valid)[0]
		spearman_valid = scipy.stats.spearmanr(y_valid,y_hat_valid)
		muae_valid = mean_absolute_error(y_valid,y_hat_valid)
		mdae_valid = median_absolute_error(y_valid, y_hat_valid)
		expvar_valid = explained_variance_score(y_valid, y_hat_valid)

		## Debugging
		print("Valid y_range:[{},{}], yhat_range:[{},{}]".format(
			     min(y_valid),max(y_valid),round(min(y_hat_valid),2),round(max(y_hat_valid),2)))
		print("Validation MSE: {}, R2:{}, PS:{}, MUAE:{}, MDAE:{}, EXPS:{}".format(
					mse_valid,
					r2_valid,
					pearson_valid,
					muae_valid,
					mdae_valid,
					expvar_valid))

		valid_avg_mse += mse_valid
		valid_avg_r2 += r2_valid
		valid_avg_sp = valid_avg_sp + np.array(spearman_valid)
		valid_avg_ps = valid_avg_ps + np.array(pearson_valid)
		valid_avg_muae += muae_valid
		valid_avg_mdae += mdae_valid
		valid_avg_expvar += expvar_valid
	#//end for

	preds = {'reg':estimator,
					 'mse':valid_avg_mse/k_fold,
					 'r2':valid_avg_r2/k_fold,
					 'sp':valid_avg_sp/k_fold,
					 'ps':valid_avg_ps/k_fold,
					 'muae':valid_avg_muae/k_fold,
					 'mdae':valid_avg_mdae/k_fold,
					 'expvar':valid_avg_expvar/k_fold,
					 'imp':train_avg_imp/k_fold
					 }
	gc.collect()
	return preds
#//end crossval_predict

#-----------------------------------------------------------------------------------------
# Compute the Mean Squared Errors (MSE) 
# Input:
#   regressor - regression model name
#   param - parameters used in the regression model
#   has_importances - the model output feature importances
#   train_data - training set
#   test_data  - testing set
#   n_jobs - number of threads in parallel
# Output:
#   mse_y_hat - MSE of the predicted result
#------------------------------------------------------------------------------------------
def compute_mse(regressor, params, has_importances, train_data, test_data, n_jobs, k_fold):
	##
	## Note: the last three columns (samplefreq, rawcounts, rname) are not features, instead
	##   [samplefreq] - the main response variable (default)
	##   [rawcounts] - the optional response variable (only raw counts, no sample freq is counted)
	##   [rname] - hash labels for "hash" sampling when splitting training data into 
	##   training_subset and validation subset in cross-validation step
	##
	## In our regression recurrence model, the sample frequency is an impotant indicator because
	## in recent research paper by Weihold and Melton, the "driver"/regulatory hot noncoding
	## mutations are frequent because of the sample frequency.  
	##

	## Use sample frequency as recurrence response variable
	train_x = train_data[:, 0:-3]
	train_y = train_data[:, -3]
	train_labels = train_data[:,-1]
	test_x = test_data[:, 0:-3]
	test_y = test_data[:, -3]

	## Standarize/normalize the data
	## SVR does NOT do any scaling of the data
	train_x = (train_x - np.mean(train_x)) / np.std(train_x)
	test_x = (test_x - np.mean(test_x)) / np.std(test_x)

	## Construct the shuffled dataset against the test dataset
	shuffle_x = np.copy(test_x)

	feature_importances = []
	predsList = []

	## Learn the model with different parameters
	for param in params[regressor]:
		preds = {} #//learned model metrics for each param

		# Record the training start time
		start_time = time()
		print("train:{},{}; test:{},{}".format(
				train_x.shape,len(train_y),test_x.shape,len(test_y)))
		print("start training model {}: param {}, time:{}".format(
					regressor,param,
					datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")))

		## Check regressor type and initialize the regressor object
		if(regressor == 'sgd'):
			# Does not hurt to try a linear regression model
			reg = linear_model.SGDRegressor(loss = 'squared_loss',
												 n_iter = param,
												 penalty = 'l2',
												 l1_ratio = 0.15,
												 epsilon = 0.1,
												 eta0 = 0.01,
												 verbose = 0)
		elif(regressor == 'rf'):
			# random forest regressor
			reg = RandomForestRegressor(n_estimators = param, 
																	max_features = train_x.shape[1], 
																	criterion = 'mse', 
																	oob_score = True,
																	verbose = 0, 
																	min_samples_split = 2,
																	warm_start=False,
																	n_jobs = n_jobs)
		elif(regressor == 'ada'):
			# Ada Boosting with DecisionTree Stump
			rng = np.random.RandomState(1000)
			reg = AdaBoostRegressor(DecisionTreeRegressor(
																criterion = 'mse',
																max_features = 'auto',
																max_depth = 4,
																min_samples_split = 2
																),
															n_estimators  = param, #number of DT stumps
															learning_rate = 0.1,
															loss = 'square', #loss function used for weights update
															random_state = rng)
		elif(regressor == 'knn'):
			# K-Nearest Neighbor(KNN) Regressor
			reg = KNeighborsRegressor(n_neighbors = param, 
			  	                      weights = 'uniform', 
				                        algorithm = 'auto', 
			    	                    n_jobs = n_jobs)
		elif(regressor == 'svr'):
			# Support Vector Regressor(SVR)
			reg = SVR(kernel = 'rbf', 
			        	gamma = 'auto', 
			        	max_iter = -1, 
			        	epsilon = 0.1, 
			        	verbose = True,
			        	shrinking = True,
			        	cache_size = 6000,
			        	C = param) #<= penalty of the error 
		elif(regressor == 'gradient'):
			# Gradient Boosting Regression
			l_params = {'n_estimators': param,
									'max_depth': 4,
									'min_samples_split': 2,
									'learning_rate': 0.1,
									'loss': 'ls',
									'subsample':0.5}
			reg = GradientBoostingRegressor(**l_params)
		else:
			raise Exception("No regressor set.")
		#//end regressor initialization

		## Training with 10-Fold Cross-Validation
		preds = crossval_predict(reg,k_fold,train_x,train_y,train_labels,has_importances)
		predsList.append(preds)

		print("end training model {}: param {}, seconds:{}".format(
		    	regressor,param, str(round(time()-start_time,3))))
	#//end for param

	## So far we already tried out all parameters in learning and for each
	## parameter, we collected the model evaluation metrics including:
	##  - Mean Squared Error (MSE)
	##  - Mean Absolute Error (MuAE)
	##  - Median Absolute Error (MdAE)
	##  - Explained Variance Score (ExpVar)
	##  - R-Square (R2) Correlation
	##  - Spearman Correlation Coefficient
	##  - Pearson Correlation Coefficient
	##  - ...
	## Now we need to choose the parameter that fits the regression model
	## with the minimum MSE
	mse = [pred['mse'] for pred in predsList]
	reg_best = (predsList[mse.index(min(mse))])['reg']

	## Debugging
	print("Best param: {}, Valid MSE: {}".format(mse.index(min(mse)),[round(i,3) for i in mse]))

	## Extract the training metrics
	train_metrics = {'model': reg_best,
									 'best': (params[regressor])[mse.index(min(mse))],
									 'mse': [pred['mse'] for pred in predsList],
									 'r2': [pred['r2'] for pred in predsList],
									 'sp': [pred['sp'] for pred in predsList],
									 'ps': [pred['ps'] for pred in predsList],
									 'muae':[pred['muae'] for pred in predsList],
									 'mdae':[pred['mdae'] for pred in predsList],
									 'expvar':[pred['expvar'] for pred in predsList],
									 'imp': [pred['imp'] for pred in predsList]
								}

	## Apply the whole training dataset to the best learned model
	reg_best_fit = reg_best.fit(train_x,train_y)

	## Here we may need to save the model object to disk
	## in case for use of predictions of new features data
	## other than testing data

	## Deploy the model to the testing set
	y_hat_test = reg_best_fit.predict(test_x)

	## Compute the testing metrics
	test_metrics = {'mse': mean_squared_error(test_y,y_hat_test),
									'r2': r2_score(test_y,y_hat_test,multioutput='uniform_average'),
									'sp': np.array(scipy.stats.spearmanr(test_y,y_hat_test)),
									'ps': np.array(scipy.stats.pearsonr(test_y,y_hat_test)),
									'muae': mean_absolute_error(test_y,y_hat_test),
									'mdae': median_absolute_error(test_y, y_hat_test),
									'expvar': explained_variance_score(test_y, y_hat_test)
								}

	## Shuffle the features in test dataset and predict against the shuffled data
	## using the best learned model
	shuffle_mse = []
	shuffle_r2 = []
	shuffle_sp = []
	shuffle_ps = []
	shuffle_muae = []
	shuffle_mdae = []
	shuffle_expvar = []
	for i in range(5):
		shuffle_x = np.copy(test_x)
		for j in range(test_x.shape[1]):
			np.random.shuffle(shuffle_x[:,j])

		y_hat_shuffle = reg_best_fit.predict(shuffle_x)
		shuffle_mse.append(mean_squared_error(test_y,y_hat_shuffle))
		shuffle_r2.append(r2_score(test_y,y_hat_shuffle,multioutput='uniform_average'))
		shuffle_sp.append(np.array(scipy.stats.spearmanr(test_y,y_hat_shuffle)))
		shuffle_ps.append(np.array(scipy.stats.pearsonr(test_y,y_hat_shuffle)))
		shuffle_muae.append(mean_absolute_error(test_y,y_hat_shuffle))
		shuffle_mdae.append(median_absolute_error(test_y, y_hat_shuffle))
		shuffle_expvar.append(explained_variance_score(test_y, y_hat_shuffle))
	#//end for

	shuffle_metrics = {'mse': shuffle_mse,
										 'r2': shuffle_r2,
										 'sp': shuffle_sp,
										 'ps': shuffle_ps,
										 'muae': shuffle_muae,
										 'mdae': shuffle_mdae,
										 'expvar':shuffle_expvar
									}

	## Debugging
	print("Test y_range:[{},{}], yhat_range:[{},{}]".format(
		    min(test_y),max(test_y),round(min(y_hat_test),2),round(max(y_hat_test),2)))
	print("Test MSE:{}, r2:{}, ps:{}, MuAE:{}, MdAE:{}, ExpVar:{}".format(
				test_metrics['mse'], 
				test_metrics['r2'], 
				test_metrics['ps'],
				test_metrics['muae'],
				test_metrics['mdae'],
				test_metrics['expvar']))

	## Return the train/test metrics
	metrics = {'train':train_metrics, 'test':test_metrics, 'shuffle': shuffle_metrics}

	## Garbage collection to free memory
	predsList = []
	gc.collect()

	return metrics
#//end function compute_mse


#------------------------------------------------------------------------------------
# main function
#------------------------------------------------------------------------------------
def main():
	in_path = "/home/tanjinxu/Project/noncoding/data/for_python/"
	out_path = "/home/tanjinxu/Project/noncoding/sklearn_out/dumps/filtered/"
	
	in_str =  in_path + "nc.data.sampled.filtered"
	out_pref = out_path + "nc.res.freq.filtered"
	day_str = "421"

	# Number of threads in parallel
	n_jobs = 5
	k_fold = 5 #K-fold cross-validation

	#sample_rates = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
	sample_rates = [0.3,0.4,0.5]
	sample_iters = 5

	regressors = ['gradient']
	#regressors = ['rf','ada','gradient']

	## Indicators which model outputs variable importances, maily GLM and Trees Models
	importances = {'sgd':False,'rf':True,'ada':True,'knn':False,'svr':False,'gradient':True}

	## Parameters to be tunned for each regression model to get best performance
	## using cross-validation method
	'''
	params = {'sgd': [5, 10, 20, 50, 100, 150, 200],   #number of iterations
						'rf': [100,200,300,400,500,600,700],     #number of trees
						'ada':[50,100,150,200,250,300,350],      #number of trees
						'knn': [5,7,9,10,15,20,30],              #number of neighbours
						'svr': [1.0,5.0,1e1,5e1],                #penalties
						'gradient':[50,100,200,250,300,400,500]  #number of boosting iterations
					}
	'''
	params = {'sgd': [5, 50, 100, 200],    
						'rf': [100,300,500,700],     
						'ada':[50,100,200,300],      
						'knn': [5,7,9,10,15,20,30],             
						'svr': [1.0,5.0,1e1,5e1],               
						'gradient':[50,100,200,400]  
					}

	##
	## Note: since the original whole dataset size is huge (~4.5M) and the project time
	## is limitted, and a huge memory and computing resource is required (e.g. randomForest),
	## here we uniformly sample a subset from the whole dataset, as statistically there
	## should be no big difference. However, to verify the data is uniformly distributed
	## (i.e. no preference on some speical mutation events), we tried different sampling
	## rate, for each sampling rate, we run couples of times. 
	##
	for rate in sample_rates:
		for i in range(sample_iters):
			verstr = "."+str(rate)+"."+str(i+1)+"."
			in_data_file = in_str + verstr + day_str + ".save"
			data = cPickle.load(open(in_data_file,'rb'))
			train_data = data['train_data']
			test_data = data['test_data']

			## Try all proposed regression models
			for regressor in regressors:
				## Learn the regression model with different parameter configurations
				## and compute the errors for model performance evaluation, including
				## validation errors, test errors, and prediction errors on scrambled
				## features (shuffled test data)
				metrics = compute_mse(regressor,
				      		      	    params,
															importances[regressor],
															train_data, 
															test_data,
															n_jobs,
															k_fold)

				## Dump the important in-memory objects including the fitted best models for 
				## later/future reference, which could be done simply by loading everything
				## in the dumpped file into memory 
				dumpf_name = out_path+"nc.model.res.sampled"+ verstr + regressor
				cPickle.dump(metrics, open(dumpf_name,'wb'), protocol=cPickle.HIGHEST_PROTOCOL)

				out_str = out_pref + verstr + regressor
				out_train_err_file = out_str + ".train.err.tsv"
				out_train_imp_file = out_str + ".train.imp.tsv"
				out_test_err_file = out_str + ".test.err.tsv"

				## Also write results to tsv/csv for easier plotting
				## For training ==>
				## 1. params: [...]
				## 2. bestpar: xxx
				## 3. fit.mse:[...]
				## 4. fit.mean_abs_err: [...]
				## 5. fit.median_abs_err: [...]
				## 6. fit.r2_score: [...]
				## 7. fit.exp_score: [...]
				train_metrics = metrics['train']
				result_arr = np.asarray([params[regressor], 
																np.tile(train_metrics['best'],len(params)),
																train_metrics['mse'],
																train_metrics['muae'],
																train_metrics['mdae'],
																train_metrics['r2'],
																train_metrics['expvar']])
				row_names = np.array(['params','best',
															'train.mse','train.mean_ae','train.median_ae',
															'train.r2','train.expvar'])
				np.savetxt(out_train_err_file, np.c_[row_names,result_arr], delimiter="\t", fmt="%s")
				np.savetxt(out_train_imp_file, np.asarray(train_metrics['imp']), delimiter="\t", fmt="%s")

				test_metrics = metrics["test"]
				shuff_metrics = metrics["shuffle"]
				test_row_names = np.array(['test.mse','test.mean_ae','test.median_ae','test.r2','test.expvar'])
				shuff_result_arr = np.asarray([shuff_metrics['mse'],
																		 shuff_metrics['muae'],
																		 shuff_metrics['mdae'],
																		 shuff_metrics['r2'],
																		 shuff_metrics['expvar']])
				test_result_arr = np.asarray([test_metrics['mse'],
																			test_metrics['muae'],
																			test_metrics['mdae'],
																			test_metrics['r2'],
																			test_metrics['expvar']])
				final_test_result = np.c_[test_row_names,test_result_arr]
				np.savetxt(out_test_err_file, np.c_[final_test_result,shuff_result_arr],delimiter="\t", fmt="%s")
			#//end for
		#//end for
	#//end for
#//end main

if __name__ == "__main__":
	main()
