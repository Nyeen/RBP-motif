import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import CNN_utils as ut
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, Dropout, Conv2D, MaxPool2D, Flatten
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
# from tensorflow.keras.regularizers import L2
from tensorflow.keras.constraints import MaxNorm

from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix


# We know the cut-off nucleotide lengths are 3 to 22
flanking_sequence = 0 # True or False
flanking_sequence_length = 30
cut_off = 22 + flanking_sequence*flanking_sequence_length
mode = 'pstack'




if flanking_sequence:
    try:
        dfp = pd.read_csv('positive_wf_'+str(flanking_sequence_length)+'.csv')
    except:
        dfp = pd.read_csv('positive.csv')
        dfp = ut.add_flanking_sequence(dfp, length = flanking_sequence_length)
        dfp.to_csv('positive_wf_'+str(flanking_sequence_length)+'.csv', index = False)
        
    try:
        dfn = pd.read_csv('negative_wf_'+str(flanking_sequence_length)+'.csv')
    except:
        dfn = pd.read_csv('negative.csv')
        dfn = ut.add_flanking_sequence(dfn, length = flanking_sequence_length)
        dfn.to_csv('negative_wf_'+str(flanking_sequence_length)+'.csv', index = False)

else:
    dfp = pd.read_csv('positive.csv')
    dfn = pd.read_csv('negative.csv')




x_positive, y_positive = ut.get_data_label(dfp, label = 1, cut_off = cut_off, mode = mode)
x_negative, y_negative = ut.get_data_label(dfn, label = 0, cut_off = cut_off, mode = mode)

x_positive = x_positive.astype('float32')
x_negative = x_negative.astype('float32')
random.Random(8).shuffle(x_positive)
random.Random(8).shuffle(x_negative)

positive = list(zip(x_positive, y_positive))
negative = list(zip(x_negative, y_negative))


x_train, y_train, x_test, y_test = ut.prepare_train_test(positive, negative, split = 0.2)
if len(np.shape(x_train)) < 4:
    x_train = x_train.reshape((x_train.shape[0], x_train.shape[1], x_train.shape[2], 1))
    x_test = x_test.reshape((x_test.shape[0], x_test.shape[1], x_test.shape[2], 1))
    



# building a linear stack of layers with the sequential model
model = Sequential()
# convolutional layer1
model.add(Conv2D(32, kernel_size=(4,4), strides=(1,1), padding='same', activation='relu',
                 input_shape=np.shape(x_train)[1:])) #
# model.add(MaxPool2D(pool_size=(2,2), strides = 2))
model.add(Dropout(0.5))
# convolutional layer2
model.add(Conv2D(256, kernel_size=(4,4), strides=(1,1), padding='same', activation='relu'))
# model.add(MaxPool2D(pool_size=(2,2), strides = 2))
model.add(Dropout(0.5))
# convolutional layer3
model.add(Conv2D(160, kernel_size=(4,4), strides=(1,1), padding='same', activation='relu'))
# model.add(MaxPool2D(pool_size=(2,2), strides = 2))
model.add(Dropout(0.5))
# convolutional layer4
model.add(Conv2D(192, kernel_size=(3,3), strides=(1,1), padding='same', activation='relu'))
# model.add(MaxPool2D(pool_size=(2,2), strides = 2))
model.add(Dropout(0.5))
# convolutional layer5
model.add(Conv2D(32, kernel_size=(4,4), strides=(1,1), padding='same', activation='relu'))
# model.add(MaxPool2D(pool_size=(2,2), strides = 2))
model.add(Dropout(0.5))
# flatten output of conv
model.add(Flatten())
# Dense layer1
model.add(Dense(512, activation='relu', kernel_constraint = MaxNorm(2)))
model.add(Dropout(0.5))
# Dense layer2
model.add(Dense(128, activation='relu', kernel_constraint = MaxNorm(2)))
model.add(Dropout(0.5))
# output layer
model.add(Dense(1, activation='sigmoid'))

# compiling the sequential model
opt = Adam(learning_rate=0.0001)
model.compile(loss='binary_crossentropy', metrics=['accuracy'], optimizer=opt)

stopping = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=30)

checkpoint = ModelCheckpoint("best_model",
                             monitor = 'val_loss',
                             mode = 'min',
                             save_best_only=True)

# training the model for 200 epochs
history = model.fit(x_train, y_train, batch_size=256, epochs=200,
                    validation_data = (x_test, y_test),
                    callbacks = [stopping, checkpoint])#validation_split=0.2)





train_loss = history.history['loss']
val_loss = history.history['val_loss']
print('Validation loss:\t', min(val_loss))
# print('Validation accuracy:\t', history.history['accuracy'][np.argmin(val_loss)])

loaded_model = load_model('best_model')
y_prob = loaded_model.predict(x_test, verbose=0)
y_prob = y_prob[:, 0]
y_pred = np.round(y_prob).astype('int32')

# accuracy: (tp + tn) / (p + n)
accuracy = accuracy_score(y_test, y_pred)
print('Accuracy: %f' % accuracy)
# precision tp / (tp + fp)
precision = precision_score(y_test, y_pred)
print('Precision: %f' % precision)
# recall: tp / (tp + fn)
recall = recall_score(y_test, y_pred)
print('Recall: %f' % recall)
# f1: 2 tp / (2 tp + fp + fn)
f1 = f1_score(y_test, y_pred)
print('F1 score: %f' % f1)
 
# kappa
kappa = cohen_kappa_score(y_test, y_pred)
print('Cohens kappa: %f' % kappa)
# ROC AUC
auc = roc_auc_score(y_test, y_prob)
print('ROC AUC: %f' % auc)
# confusion matrix
matrix = confusion_matrix(y_test, y_pred)
print(matrix)



plt.figure(figsize = (10,7))
plt.plot(train_loss)
plt.plot(val_loss)
plt.grid()
plt.show()

