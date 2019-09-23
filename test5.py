import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Conv2D, Flatten, MaxPooling2D
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

model = tf.keras.Sequential()
# add model layers
model.add(Conv2D(64, kernel_size=3, activation='relu', input_shape=(28, 28, 1)))
model.add(MaxPooling2D(pool_size=(2, 2), strides=None, padding='valid', data_format=None))
model.add(Conv2D(32, kernel_size=3, activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2), strides=None, padding='valid', data_format=None))
model.add(Flatten())
