import cPickle
import numpy as np
import theano
import theano.tensor as T
import random
from keras.layers import containers
import keras
from keras.models import Sequential
from keras.layers.core import Dense, Flatten, AutoEncoder, Dropout
from keras.optimizers import SGD,Adam
from keras.utils import np_utils
from keras.callbacks import EarlyStopping

import scipy
from sklearn.metrics import mean_squared_error,r2_score
from sklearn.metrics import mean_absolute_error, median_absolute_error
from sklearn.metrics import explained_variance_score

random.seed(1)
np.random.seed(1)

def writeHist(history, fname):
    fhdl = open(fname,'w')
    for key in ['acc','loss','val_acc','val_loss']:
        fhdl.write('\t'.join(str(round(i,4)) for i in history[key]))
        fhdl.write('\n')
    fhdl.close()
#//end func 

def learnANN(infile, hidden_sizes, dump_out, train_out, test_out):
	data = cPickle.load(open(infile,'rb'))
	train_data = data['train_data']
	test_data = data['test_data']

	## Need normalization/standarization
	## Note: last three columns are non features: 
	## <sample_freq>, <raw_counts>, <rnames>
	## <sample_freq> - sample frequency for recurrence
	## <raw_counts> - raw counts for recurrence
	X_train = train_data[:,0:-3]
	y_train = train_data[:,-3]
	X_test = test_data[:,0:-3]
	y_test = test_data[:,-3]
	X_train = (X_train - np.mean(X_train)) / np.std(X_train)
	X_test = (X_test - np.mean(X_test)) / np.std(X_test)

	print("train:{}, test:{}".format(X_train.shape,X_test.shape))
	print y_train.shape,y_test.shape

	nb_epoch = 10
	batch_size = 50

	'''
	train_subset_y_cat = np_utils.to_categorical(train_subset_y, nb_labels)
	dev_y_cat = np_utils.to_categorical(dev_y, nb_labels)
	test_y_cat = np_utils.to_categorical(test_y, nb_labels)
	'''
	hist_list = []
	model_list = []
	avg_val_loss_list = []

	for hidden_size in hidden_sizes:
		model = Sequential()
		model.add(Dense(hidden_size, input_dim=X_train.shape[1], activation='tanh'))
		#model.add(Dropout(0.5))
		#model.add(Dense(nb_labels, activation='softmax'))
		#model.add(Dense(1,activation='linear'))
		model.add(Dense(1,input_dim=hidden_size))

		adam = Adam(lr=0.01, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
		model.compile(loss='mean_squared_error', optimizer=adam)

		print('Start training')
		hist = model.fit(X_train, y_train,batch_size=batch_size, nb_epoch=nb_epoch,
										show_accuracy=False, verbose=True, validation_split=0.1,
										validation_data=None,shuffle=True)
	
		#print(hist.history)
		#writeHist(hist.history,'noncoding.ann.base.neuros.'+str(hidden_size)+'.hist')
		hist_list.append(hist.history)
		model_list.append(model)
	#//end for

	i = 0
	j = 0
	min_loss_avg = float('inf')
	for hist in hist_list:
		#print hist
		val_loss_avg = np.mean(hist['val_loss'])
		avg_val_loss_list.append(val_loss_avg)

		if val_loss_avg < min_loss_avg:
			j = i
			min_loss_avg = val_loss_avg
		i += 1
	#//end for

	print("Min average validation loss:{} with hidden size:{}".format(val_loss_avg,hidden_sizes[j]))
	#score = model.evaluate(X_test, y_test, show_accuracy=False, verbose=0)
	#print ('Test score:', score)

	## Also dump the train/validation loss to tsv file
	np.savetxt(train_out, np.asarray([hidden_sizes,avg_val_loss_list]), delimiter="\t", fmt="%s")

	## Choose the best model of the minimum average loss
	model = model_list[j]
	best_index = j
	y_predicted = model.predict(X_test)[:,0]
	print("y_test range:({},{}), y_hat range:({},{})".format(
				min(y_test),max(y_test),min(y_predicted),max(y_predicted)))

	test_mse = mean_squared_error(y_test,y_predicted)
	test_r2 = r2_score(y_test,y_predicted)
	test_spearman = scipy.stats.spearmanr(y_test,y_predicted)
	test_pearson = scipy.stats.pearsonr(y_test,y_predicted)
	test_muae = mean_absolute_error(y_test,y_predicted)
	test_mdae = median_absolute_error(y_test,y_predicted)
	test_exps = explained_variance_score(y_test,y_predicted)

	## Also predict against the shuffled features data of test set
	shuffle_mse = []
	shuffle_r2 = []
	shuffle_sp = []
	shuffle_ps = []
	shuffle_muae = []
	shuffle_mdae = []
	shuffle_expvar = []

	for i in range(5):
		shuffle_x = np.copy(X_test)
		for j in range(X_test.shape[1]):
			np.random.shuffle(shuffle_x[:,j])
		#//end for

		y_hat_shuffle = model.predict(shuffle_x)
		shuffle_mse.append(mean_squared_error(y_test,y_hat_shuffle))
		shuffle_r2.append(r2_score(y_test,y_hat_shuffle,multioutput='uniform_average'))
		#shuffle_sp.append(np.array(scipy.stats.spearmanr(y_test,y_hat_shuffle)))
		#shuffle_ps.append(np.array(scipy.stats.pearsonr(y_test,y_hat_shuffle)))
		shuffle_muae.append(mean_absolute_error(y_test,y_hat_shuffle))
		shuffle_mdae.append(median_absolute_error(y_test, y_hat_shuffle))
		shuffle_expvar.append(explained_variance_score(y_test, y_hat_shuffle))
	#//end for

	print("MSE:{},R2:{},Spearman:{},Pearson:{}".format(test_mse,test_r2,test_spearman,test_pearson))

	dump_data = {'neurons':hidden_sizes,
							 'model': model,
							 'best':hidden_sizes[best_index],
							 'val_loss':avg_val_loss_list,
							 'mse':test_mse,
							 'r2':test_r2,
							 'spearman':test_spearman,
							 'pearson':test_pearson,
							 'muae':test_muae,
						 	 'mdae':test_mdae,
						 	 'expvar':test_exps
						}
	cPickle.dump(dump_data, open(dump_out,'wb'))

	## Also write to csv/tsv
	test_row_names = np.array(['test.mse','test.mean_ae','test.median_ae','test.r2','test.expvar'])
	test_result_arr = np.asarray([test_mse, test_muae, test_mdae, test_r2, test_exps])
	shuff_result_arr = np.asarray([shuffle_mse, shuffle_muae, shuffle_mdae, shuffle_r2, shuffle_expvar])
	final_test_result = np.c_[test_row_names,test_result_arr]
	np.savetxt(test_out, np.c_[final_test_result,shuff_result_arr],delimiter="\t", fmt="%s")
#//end function learnANN


#----------------------------------------------------------------------------
# main function
#----------------------------------------------------------------------------
def main():
	#hidden_sizes = [5,7,10,13,18,21,25]
	hidden_sizes = [2, 3, 5, 7, 9, 10, 12]

	#sample_rates = [0.3,0.4,0.5,0.7,0.8]
	sample_rates = [0.3, 0.4, 0.5]
	sample_iters = 5

	in_path = "/home/ramseylab/tj/noncoding_sae/"
	in_str = in_path + "data/nc.data.sampled.filtered"
	out_pref = in_path + "output/nc.res.freq.filtered"
	day_str = "421"

	for rate in sample_rates:
		for i in range(sample_iters):
			verstr = "."+str(rate)+"."+str(i+1)+"."
			out_str = out_pref + verstr + "ann"
			out_train_err_file = out_str + ".train.err.tsv"
			out_test_err_file = out_str + ".test.err.tsv"

			in_data_file = in_str + verstr + day_str + ".save"
			dump_out = in_path + "dumps/nc.model.res.sampled"+verstr+"ann"

			## Apply neural network regression
			learnANN(in_data_file, hidden_sizes, dump_out, out_train_err_file, out_test_err_file)
		#//enf for
	#//end for
#//end main

if __name__ == "__main__":
	main()