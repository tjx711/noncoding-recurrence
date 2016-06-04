#!python

import cPickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ticks = ["50,30","100,30",'100,50','150,100','180,100','200,100','200,150']
valid_err = [0.13172, 0.13188, 0.13172, 0.13172, 0.13172, 0.13172, 0.13172]
test_err = [0.13817,0.13899,0.13817,0.13817,0.13817,0.13817,0.13817]

fig = plt.figure()
ax = fig.add_subplot(110)
plt.title("Accuracy vs. Hidden Layer Size")

plt.plot(1-np.array(valid_err), 'ms--', label='valid_acc')
plt.plot(1-np.array(test_err), 'bd--', label='test_acc')

plt.xlabel('hidden layer size')
plt.ylabel('accuracy')
plt.xticks(range(len(ticks)),ticks,rotation=45,fontsize=10)
plt.legend(loc='best')
plt.show()

plt.savefig("dbn_accur.jpg",dpi=300,bbox_inches='tight')