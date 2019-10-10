from Bio import SeqIO
import numpy as np
import sympy
import tensorflow as tf

filename = "6mep.fasta.txt"

l = []
for seq_record in SeqIO.parse(filename, "fasta"):
    l.append(list(seq_record.seq[:]))
    break
# print(l)

vocab = ['A', ' B', ' C', ' D', ' E', ' F', ' G', ' H', ' I', ' J', ' K', ' L', ' M', ' N', ' O', ' P', ' Q', ' R',
         ' S', ' T', ' U', ' V', ' W', ' X']

t = tf.compat.v1.placeholder(dtype=tf.string, shape=(None,))
matches = tf.stack([tf.equal(t, s) for s in vocab], axis=-1)
onehot = tf.cast(matches, tf.float32)

with tf.compat.v1.Session() as sess:
    out = sess.run(onehot, feed_dict={t: l[0]})
    # print(out)
    # print(len(out[0]))

mat = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 1, 0], [1, 0, 0, 1]])  # your matrix
_, inds = sympy.Matrix(mat).T.rref()  # to check the rows you need to transpose!
print(inds)
# tf.feature_column.categorical_column_with_vocabulary_file(
#     key=feature_name_from_input_fn,
#     vocabulary_file="code.txt")
# result = tf.one_hot(
#     code,
#     23,
#     on_value=None,
#     off_value=None,
#     axis=None,
#     dtype=None,
#     name=None
# )
# print(result)
# code len = 24 => depth has to be 23
