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

import scipy
from sklearn.metrics import mean_squared_error,r2_score
from sklearn.metrics import mean_absolute_error, median_absolute_error
from sklearn.metrics import explained_variance_score

random.seed(3)
np.random.seed(3)


def writeHist(history, fname):
    fhdl = open(fname,'w')
    for key in ['acc','loss','val_acc','val_loss']:
        fhdl.write('\t'.join(str(round(i,4)) for i in history[key]))
        fhdl.write('\n')
    fhdl.close()
#//end func 

def learnSAE(infile, hidden_pairwise_layers, dump_out, train_out, test_out):
	data = cPickle.load(open(infile,'rb'))
	train_data = data['train_data']
	test_data = data['test_data']

	## Need normalization/standarization
	X_train = train_data[:,0:-3]
	y_train = train_data[:,-3]
	X_test = test_data[:,0:-3]
	y_test = test_data[:,-3]
	X_train = (X_train - np.mean(X_train)) / np.std(X_train)
	X_test = (X_test - np.mean(X_test)) / np.std(X_test)

	print("train:{}, test:{}".format(X_train.shape,X_test.shape))
	print y_train.shape,y_test.shape

	nb_epoch_pretraining = 10
	batch_size_pretraining = 50
	batch_size = 50
	nb_epoch=10

	#hidden_pairwise_layers = [(20,10),(20,5),(25,10),(25,5),(15,10),(15,5)]
	#nb_hidden_layers = [X_train.shape[1],20,10]
	#X_train_tmp = np.copy(X_train)

	hist_list = []
	model_list = []
	avg_val_loss_list = []

	for hidden_pair in hidden_pairwise_layers:
		# Layer-wise batch_size_pretraining
		encoders = []
		decoders = []

		X_train_tmp = np.copy(X_train)
		nb_hidden_layers = [X_train.shape[1],hidden_pair[0],hidden_pair[1]]

		for i, (n_in, n_out) in enumerate(zip(nb_hidden_layers[:-1], nb_hidden_layers[1:]), start=1):
			print('Training the layer {}: Input {} -> Output {}'.format(i, n_in, n_out))
			# Create AE and training
			ae = Sequential()

			if n_out >= 20:
				encoder = containers.Sequential([Dense(n_out, input_dim=n_in, activation='tanh'), Dropout(0.5)])
			else:
				encoder = containers.Sequential([Dense(n_out, input_dim=n_in, activation='tanh')])

			decoder = containers.Sequential([Dense(n_in, input_dim=n_out, activation='tanh')])  

			ae.add(AutoEncoder(encoder=encoder, decoder=decoder, output_reconstruction=False))
			sgd = SGD(lr=2, decay=1e-6, momentum=0.0, nesterov=True)
			ae.compile(loss='mse', optimizer='adam')
			ae.fit(X_train_tmp, X_train_tmp, batch_size=batch_size_pretraining, nb_epoch=nb_epoch_pretraining, verbose = True, shuffle=True)

			# Store trainined weight and update training data
			encoders.append(ae.layers[0].encoder)
			decoders.append(ae.layers[0].decoder)

			X_train_tmp = ae.predict(X_train_tmp)
		#//end for 

		################
		#End to End Autoencoder training    
		if len(nb_hidden_layers) > 2:
			full_encoder = containers.Sequential()
			for encoder in encoders:
				full_encoder.add(encoder)

			full_decoder = containers.Sequential()
			for decoder in reversed(decoders):
				full_decoder.add(decoder)

			full_ae = Sequential()
			full_ae.add(AutoEncoder(encoder=full_encoder, decoder=full_decoder, output_reconstruction=False))    
			full_ae.compile(loss='mse', optimizer='adam')

			print "Pretraining of full AE"
			full_ae.fit(X_train, X_train, batch_size=batch_size_pretraining, nb_epoch=nb_epoch_pretraining, verbose = True, shuffle=True)

		## Putting all together
		model = Sequential()
		for encoder in encoders:
			model.add(encoder)

		model.add(Dense(1,input_dim=nb_hidden_layers[-1]))
		adam = Adam(lr=0.01, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
		model.compile(loss='mean_squared_error', optimizer=adam)

		print('Start training')
		hist = model.fit(X_train, y_train,batch_size=batch_size, nb_epoch=nb_epoch,
										show_accuracy=False, verbose=True, validation_split=0.1,
										validation_data=None,shuffle=True)

		hist_list.append(hist.history)
		model_list.append(model)
	#//end for

	i = 0
	j = 0
	min_loss_avg = float('inf')
	for hist in hist_list:
		#print hist
		val_loss_avg = min(hist['val_loss'])
		avg_val_loss_list.append(val_loss_avg)

		if val_loss_avg < min_loss_avg:
			j = i
			min_loss_avg = val_loss_avg
		i += 1

	print("Min average validation loss:{} with hidden size:{}".format(min_loss_avg,hidden_pairwise_layers[j]))

	## Also dump the train/validation loss to tsv file
	np.savetxt(train_out, np.asarray([hidden_pairwise_layers,avg_val_loss_list]), delimiter="\t", fmt="%s")

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

	print("MSE:{},R2:{},Spearman:{},Pearson:{}, Exps:{}".format(test_mse,test_r2,test_spearman,test_pearson,test_exps))

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

	dump_data = {'neurons':hidden_pairwise_layers,
							 'model':model,
							 'best':hidden_pairwise_layers[best_index],
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
#//end function learSAE


#----------------------------------------------------------------------------
# main function
#----------------------------------------------------------------------------
def main():
	#hidden_pairwise_layers = [(20,10),(20,5),(25,10),(25,5),(15,10),(15,5)]
	hidden_pairwise_layers =  [(12,10),(10,8),(12,5),(10,5),(8,4),(6,3)]

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
			out_str = out_pref + verstr + "sae"
			out_train_err_file = out_str + ".train.err.tsv"
			out_test_err_file = out_str + ".test.err.tsv"

			in_data_file = in_str + verstr + day_str + ".save"
			dump_out = in_path + "dumps/nc.model.res.sampled"+verstr+"sae"

			## Apply neural network regression
			learnSAE(in_data_file, hidden_pairwise_layers, dump_out, out_train_err_file, out_test_err_file)
		#//enf for
	#//end for
#//end main

if __name__ == "__main__":
	main()