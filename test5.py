import tensorflow as tf
from Bio import SeqIO
from Bio.PDB import *
import numpy as np
import sympy


def get_edm(filename):
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(filename, pdir='.', file_format='mmCif')
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(filename, filename + '.cif')
    model = structure.get_models()
    models = list(model)
    for i in range(len(models)):
        print('model', i)
        chains = list(models[i].get_chains())
        modelAtoms = list()
        for j in range(len(chains)):
            residues = list(chains[j].get_residues())
            for k in range(len(residues)):
                atoms = list(residues[k].get_atoms())
                for l in range(len(atoms)):
                    if (len(modelAtoms) < 10) and (atoms[l].get_fullname() == 'CA'):
                        modelAtoms.append(atoms[l])
        matrix = []
        for j in range(len(modelAtoms)):
            matrix.append([0] * len(modelAtoms))
        for j in range(len(modelAtoms)):
            for k in range(len(modelAtoms)):
                matrix[j][k] = (modelAtoms[j] - modelAtoms[k]) * (modelAtoms[j] - modelAtoms[k])
    np.asarray(matrix)
    _, inds = sympy.Matrix(matrix).T.rref()
    m = []
    for ind in inds:
        m.append(matrix[ind])
    return np.asarray(m)


# TODO: opening of fasta database
filename = "6mep.fasta.txt"
structure_name = "6mep"



list_of_fasta_arrs = []
for seq_record in SeqIO.parse(filename, "fasta"):
    list_of_fasta_arrs.append(list(seq_record.seq[:]))
    break
print(list_of_fasta_arrs)
vocab = ['A', ' B', ' C', ' D', ' E', ' F', ' G', ' H', ' I', ' J', ' K', ' L', ' M', ' N', ' O', ' P', ' Q', ' R',
         ' S', ' T', ' U', ' V', ' W', ' X']

t = tf.compat.v1.placeholder(dtype=tf.string, shape=(None,))
matches = tf.stack([tf.equal(t, s) for s in vocab], axis=-1)
onehot = tf.cast(matches, tf.float32)

# TODO: make a loop and change the indexes in list_of_fasta_arrs
with tf.compat.v1.Session() as sess:
    onehot_encoded_sequence = sess.run(onehot, feed_dict={t: list_of_fasta_arrs[0]})
ONE_HOT_ENCODING_LEN = len(onehot_encoded_sequence[0])


def cnn_model_fn(features, labels, mode):
    """Model function for CNN."""
    # Input Layer
    input_layer = tf.reshape(features["x"], [-1, features.shape[0], ONE_HOT_ENCODING_LEN, 1])

    # Convolutional Layer #1
    conv1 = tf.layers.conv2d(
        inputs=input_layer,
        filters=32,
        kernel_size=[5, 5],
        padding="same",
        activation=tf.nn.relu)

    # Pooling Layer #1
    pool1 = tf.layers.max_pooling2d(inputs=conv1, pool_size=[2, 2], strides=2)

    # Convolutional Layer #2 and Pooling Layer #2
    conv2 = tf.layers.conv2d(
        inputs=pool1,
        filters=64,
        kernel_size=[5, 5],
        padding="same",
        activation=tf.nn.relu)

    pool2 = tf.layers.max_pooling2d(inputs=conv2, pool_size=[2, 2], strides=2)

    # # Dense Layer
    # pool2_flat = tf.reshape(pool2, [-1, 7 * 7 * 64])
    # dense = tf.layers.dense(inputs=pool2_flat, units=1024, activation=tf.nn.relu)
    # dropout = tf.layers.dropout(
    #     inputs=dense, rate=0.4, training=mode == tf.estimator.ModeKeys.TRAIN)

    # Logits Layer
    logits = tf.layers.dense(inputs=pool2, units=features.shape[0] * 5)
    print(logits.shape)

    # making a right result to compare with it
    example = get_edm(structure_name)
    mse = tf.reduce_mean(tf.squared_difference(logits, example.reshape(1,)))

    # TODO: loss function and training

    # predictions = {
    #     # Generate predictions (for PREDICT and EVAL mode)
    #     "classes": tf.argmax(input=logits, axis=1),
    #     # Add `softmax_tensor` to the graph. It is used for PREDICT and by the
    #     # `logging_hook`.
    #     "probabilities": tf.nn.softmax(logits, name="softmax_tensor")
    # }
    #
    # if mode == tf.estimator.ModeKeys.PREDICT:
    #     return tf.estimator.EstimatorSpec(mode=mode, predictions=predictions)
    #
    # # Calculate Loss (for both TRAIN and EVAL modes)
    # loss = tf.losses.sparse_softmax_cross_entropy(labels=labels, logits=logits)
    #
    # # Configure the Training Op (for TRAIN mode)
    # if mode == tf.estimator.ModeKeys.TRAIN:
    #     optimizer = tf.train.GradientDescentOptimizer(learning_rate=0.001)
    #     train_op = optimizer.minimize(
    #         loss=loss,
    #         global_step=tf.train.get_global_step())
    #     return tf.estimator.EstimatorSpec(mode=mode, loss=loss, train_op=train_op)
    #
    # # Add evaluation metrics (for EVAL mode)
    # eval_metric_ops = {
    #     "accuracy": tf.metrics.accuracy(
    #         labels=labels, predictions=predictions["classes"])
    # }
    # return tf.estimator.EstimatorSpec(
    #     mode=mode, loss=loss, eval_metric_ops=eval_metric_ops)
