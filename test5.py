import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Conv2D, Flatten, MaxPooling2D
from Bio import SeqIO


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

print(get_edm("6mep"))
print(get_edm("6mep").shape)
# def build_model():
#     model = tf.keras.Sequential()
#     # add model layers
#     model.add(Conv2D(64, kernel_size=3, activation='relu', input_shape=(28, , 1)))
#     model.add(MaxPooling2D(pool_size=(2, 2), strides=None, padding='valid', data_format=None))
#     model.add(Conv2D(32, kernel_size=3, activation='relu'))
#     model.add(MaxPooling2D(pool_size=(2, 2), strides=None, padding='valid', data_format=None))
#     model.add(Flatten())
#
#     optimizer = tf.keras.optimizers.RMSprop(0.001)
#
#     model.compile(loss='mse',
#                   optimizer=optimizer,
#                   metrics=['mae', 'mse'])
#     return model
#
#
# example = get_edm(structure_nam.reduce_mean(tf.squared_difference(logits, example.reshape(1,)))
