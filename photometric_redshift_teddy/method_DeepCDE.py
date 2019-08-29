import tensorflow as tf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import deepcde
import sklearn.metrics
import h5py

run_name = "teddy-deepcde-5"

# Parameters
learning_rate = 1e-2
initial_bias = 0.1
weight_sd = 0.1
n_epochs = 10000
n_hidden = 32
n_log = 25
n_basis = 31
logdir = ".tboard"
logfile = "%s/%s" % (logdir, run_name)

datadir = "data/"

with h5py.File(datadir + "processed.hdf5", "r") as f:
    x_train = f["x_train"][:].T
    y_train = f["y_train"][:].T
    x_test = f["x_test"][:].T
    y_test = f["y_test"][:].T

n_grid = 200
y_grid = np.linspace(0, 1, n_grid)

y_min = np.min(y_grid)
y_max = np.max(y_grid)

# Project to 0, 1
y_train = deepcde.box_transform(y_train, y_min, y_max)
y_test = deepcde.box_transform(y_test, y_min, y_max)
y_grid = deepcde.box_transform(y_grid, y_min, y_max)
shrink_factor = np.prod(y_max - y_min)

# Calculate basis
# basis = deepcde.bases.DiffusionBasis(y_train[:1000, :], n_basis)
# basis.estimate_orthogonalization(basis.evaluate(y_grid))
basis = deepcde.bases.CosineBasis(n_basis)

## Stuff
n_train = x_train.shape[0]
n_feat = x_train.shape[1]
n_test = x_test.shape[0]

# Create basis functions
y_grid_basis = basis.evaluate(y_grid)
y_train_basis = basis.evaluate(y_train)
y_test_basis = basis.evaluate(y_test)

## Define the Graph
tf.set_random_seed(4)

x = tf.placeholder(tf.float32, [None, n_feat])
y_basis = tf.placeholder(tf.float32, [None, y_grid_basis.shape[1] - 1])

# Architecture
def relu_layer(inputs, n_nodes, weight_sd, initial_bias, name="relu"):
    with tf.name_scope(name):
        W = tf.Variable(tf.random_normal([int(inputs.shape[1]), n_nodes],
                                         stddev=weight_sd), name="W")
        b = tf.Variable(tf.fill([n_nodes], initial_bias), name="b")
        outputs = tf.nn.relu(tf.matmul(inputs, W) + b)
        return outputs

# Normalize inputs
x_mu = x_train.mean(axis = 0)
x_sd = x_train.std(axis = 0)
train_dict = {x:(x_train - x_mu) / x_sd, y_basis:y_train_basis[:,1:]}
test_dict = {x:(x_test - x_mu) / x_sd, y_basis:y_test_basis[:,1:]}

marginal_beta = y_train_basis.mean(axis=0)
## Define the model
y1 = relu_layer(x, n_hidden, weight_sd, initial_bias)
y2 = relu_layer(y1, n_hidden*2, weight_sd, initial_bias)
y3 = relu_layer(y2, n_hidden, weight_sd, initial_bias)
beta = deepcde.cde_layer(y3, weight_sd, np.float32(marginal_beta))
loss = deepcde.cde_loss(beta, y_basis, shrink_factor)

optimizer = tf.train.AdamOptimizer(learning_rate).minimize(loss)

init = tf.global_variables_initializer()

## Logging
merged_summary = tf.summary.merge_all()
training_loss = tf.summary.scalar("Train_loss", loss)
test_loss = tf.summary.scalar("Test_loss", loss)

## Train the model
with tf.Session() as sess:
    if logfile is not None:
        write = tf.summary.FileWriter(logfile, sess.graph)
            
        sess.run(init)
        for epoch in range(n_epochs):
            if epoch % 100 == 0:
                print(epoch)
            sess.run(optimizer, feed_dict = train_dict)

            if logfile is not None and epoch % n_log == n_log - 1:
                # Logging
                write.add_summary(sess.run(training_loss, feed_dict=train_dict), epoch)
                write.add_summary(sess.run(test_loss, feed_dict=test_dict), epoch)
                        
        cdes = deepcde.cde_predict(sess, beta, y_min, y_max, y_grid, y_grid_basis,
                                                   test_dict)

        ## Write out to file
        fname = datadir + "DeepCDE.hdf5"
        with h5py.File(fname, 'w') as f:
            f.create_dataset("/y_grid", data=y_grid)
            f.create_dataset("/y_true", data=y_test)
            f.create_dataset("/cde", data=cdes.T)
