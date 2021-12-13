#!/usr/bin/python
# coding = utf-8

import numpy as np
import matplotlib.pyplot as plt
from keras.layers import Input, Dense
from keras.models import Model
from keras.utils import np_utils
from keras import optimizers
import os

os.environ['CUDA_VISIBLE_DEVICES'] = '/gpu:1'

if __name__ == '__main__':
    # ISIh = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    ISIh = [0.1, 0.3, 0.5, 0.7, 0.9]
    SNR = np.arange(16, 17, 2)
    # SNR = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    # load train data
    # NFeature = 2
    NClass = 16
    # NTrain = 10000
    for idx in range(len(ISIh)):
        for ii in range(len(SNR)):
            train_data = np.loadtxt('./train/traindata' + str(ISIh[idx]) + '-' + str(SNR[ii]) + 'dB.csv', delimiter=',',
                                    dtype=float)
            valid_data = np.loadtxt('./train/validdata' + str(ISIh[idx]) + '-' + str(SNR[ii]) + 'dB.csv', delimiter=',',
                                    dtype=float)
            ydata = train_data[:, -1]
            xdata = np.delete(train_data, -1, axis=1)
            index = np.arange(ydata.shape[0])
            np.random.shuffle(index)
            xdata = xdata[index, :]  # 行数打乱 列数不变
            ydata = ydata[index]
            ydata = np_utils.to_categorical(ydata, NClass)  # 这里是one-hot编码，在某种类别下为1 ，其他的值为0

            yVal = valid_data[:, -1]
            yVal = np_utils.to_categorical(yVal, NClass)
            xVal = np.delete(valid_data, -1, axis=1)  # 读入验证数据集到数组中

            # 模型开始训练 是keras 平台的网络构建的方式
            _in_ = Input(shape=(2,))
            ot = Dense(100, activation='sigmoid')(_in_)
            # ot = Dense(32, activation='relu')(ot)
            _out_ = Dense(NClass, activation='softmax')(ot)
            model = Model(_in_, _out_)
            adamax = optimizers.adamax(lr=0.01)
            model.compile(loss='categorical_crossentropy',  # 定义loss 是交叉熵函数
                          optimizer='adamax',  # 优化器是自适应 优化器
                          metrics=['categorical_accuracy'])  # 将训练的准确率 保存在矩阵中
            # 保存下来训练过程中数据
            history = model.fit(xdata, ydata,
                                epochs=300,
                                batch_size=16,
                                validation_data=(xVal, yVal),
                                shuffle=False)
            print("evaluate the model - train_set:")
            scores = model.evaluate(xdata, ydata)  # 这里是对训练进行评估 其实是训练的数据都保存在里面
            print('loss: ', str(scores[0]))
            print("%s: %.2f%%" % (model.metrics_names[1], scores[1] * 100))
            model.summary()
            epoch = len(history.history['loss'])
            model.save('./inteference-2000/dnn_model' + str(ISIh[idx]) + '-' + str(SNR[ii]) + 'dB.h5')

            plt.plot(range(epoch), history.history['loss'], label='loss')
            # 以下是 画出 训练的loss值 和 准确率变化图
            plt.xlabel("epoch")
            plt.ylabel("loss")
            plt.title("The loss with the epoch h=" + str(ISIh[idx]) + '-' + str(SNR[ii]) + 'dB')
            plt.show()
