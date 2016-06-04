#!python

import rpy2.robjects as robjects
import pandas.rpy.common as com
import pandas as pd
import numpy as np
import cPickle
import random

## load .RData and converts to pd.DataFrame
robj = robjects.r.load('train_test.rda')
# iterate over datasets the file
'''
for sets in robj:
    myRData = com.load_data(sets)
    # convert to DataFrame
    if not isinstance(myRData, pd.DataFrame):
        myRData = pd.DataFrame(myRData)
        print(myRData.shape)

## save pd.DataFrame to R dataframe
df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6], 'C':[7,8,9]},
                index=["one", "two", "three"])
r_dataframe = com.convert_to_r_dataframe(df)
'''

myRData = com.load_data(robj[0])
train_data = np.array(pd.DataFrame(myRData['train']))
test_data = np.array(pd.DataFrame(myRData['test']))

## Sample the data set (~2.5M -> sampl_size)
train_sampl_size = 1000000
test_sampl_size = 100000
train_data_subset = train_data[random.sample(range(train_data.shape[0]),train_sampl_size),]
test_data_subset = test_data[random.sample(range(test_data.shape[0]),test_sampl_size),]

print('train:{}, test:{}'.format(train_data.shape, test_data.shape))
print('train_subset:{},test_subset:{}'.format(train_data_subset.shape,test_data_subset.shape))

dump_data = {'train_data': train_data,
             'test_data': test_data}

dump_data_subset = {'train_data': train_data_subset,
										'test_data': test_data_subset}

cPickle.dump(dump_data,open('noncoding.data.py.save','wb'),protocol=cPickle.HIGHEST_PROTOCOL)
cPickle.dump(dump_data_subset,open(noncoding.data.py.1m.0.1m.save,'wb'),protocol=cPickle.HIGHEST_PROTOCOL)