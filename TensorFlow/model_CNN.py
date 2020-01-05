from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv2D, Dense, Activation, Flatten
from sklearn.model_selection import train_test_split 
import warnings
from Bio.PDB import *
from Bio import SeqIO
import numpy as np
import os
import shutil
import tensorflow as tf

warnings.filterwarnings('ignore')
warnings.filterwarnings('ignore', category=DeprecationWarning)

S_LEN = 35


# TODO: split dataset into training and testing parts
def split_combined(combined):
    train = combined[:1460]
    test = combined[1460:]

    return train, test


def get_matrix(filename):
    a = np.array([])
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(filename, pdir='.', file_format='mmCif')
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure(filename, filename + '.cif')
    except KeyError:
        return a
    model = structure.get_models()
    models = list(model)
    chains = list(models[0].get_chains())
    modelatoms = []
    for j in range(len(chains)):
        residues = list(chains[j].get_residues())
        for k in range(len(residues)):
            atoms = list(residues[k].get_atoms())
            for l in range(len(atoms)):
                if atoms[l].get_fullname() == 'CA':
                    modelatoms.append(atoms[l])
    if len(modelatoms) == 0:
        return a
    else:
        matrix = []
        for j in range(S_LEN):
            matrix.append([0] * S_LEN)
        for j in range(len(modelatoms)):
            for k in range(len(modelatoms)):
                matrix[j][k] = (modelatoms[j] - modelatoms[k]) * (modelatoms[j] - modelatoms[k])
        m = np.asarray(matrix)
        return m


def get_training_data():
    path2base = "C:/Практика Biocad/PDB/resultbase"
    path2python = "C:/Практика Biocad/Python Projects"
    f = open("proteins_names.txt", "r+")
    p_names = f.readlines()
    number_of_proteins = len(p_names)
    f.close()
    data = {}
    for protein_num in range(0, number_of_proteins):
        protein_id = p_names[protein_num].replace('\n', '')
        shutil.copyfile(path2base + "/" + protein_id + ".cif", path2python + "/" + protein_id + ".cif")
        m = get_matrix(protein_id)
        if m.shape != 0:
            os.remove(path2python + "/" + protein_id + ".cif")
            print("call for: " + protein_id)
            f = get_fasta_onehot(get_fasta_onehot(protein_id))
            if f.shape != 0:
                data[f] = m
    dataset = tf.data.Dataset.from_tensor_slices((f, m))
    return dataset


def get_fasta_onehot(protein_id):
    protein_name = str(protein_id) + ".fasta.txt"
    list_of_fasta_arrs = []
    print(type(protein_name))
    print(protein_name)
    a = np.array([])
    try:
        for seq_record in SeqIO.parse(protein_name, "fasta"):
            list_of_fasta_arrs.append(list(seq_record.seq[:]))
            break
    except OSError:
        return a
    for i in range(len(list_of_fasta_arrs[0]) - 1, S_LEN):
        list_of_fasta_arrs[0].append('-1')

    aminoacids_voc = ['-1', 'A', 'R', 'N', 'D', 'B', 'C', 'E', 'Q', 'Z', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
                      'T', 'W', 'Y', 'V']
    numbers = [i for i in range(len(aminoacids_voc))]
    v = dict(zip(aminoacids_voc, numbers))
    indices = []
    for i in list_of_fasta_arrs[0]:
        indices.append(v[i])

    onehot = tf.one_hot(indices, len(aminoacids_voc) - 1, axis=0)
    return onehot


print(get_training_data())

# Get data:

train, target = split_combined(get_training_data())

# Build model

NN_model = Sequential()

# The Input Layer :
NN_model.add(Conv2D(64, kernel_size=3, activation='relu', input_shape=(28, 28, 1)))
NN_model.add(Conv2D(32, kernel_size=3, activation='relu'))
NN_model.add(Flatten())
# The Hidden Layers :

NN_model.add(Dense(32, kernel_initializer='normal', activation='relu'))
NN_model.add(Dense(32, kernel_initializer='normal', activation='relu'))
NN_model.add(Dense(32, kernel_initializer='normal', activation='relu'))

# The Output Layer :
# TODO: add S_LEN here too
NN_model.add(Dense(1125, kernel_initializer='normal', activation='linear'))

# Compile the network :
NN_model.compile(loss='mean_absolute_error', optimizer='adam', metrics=['mean_absolute_error'])
NN_model.summary()

# CallBacks:
checkpoint_name = 'Weights-{epoch:03d}--{val_loss:.5f}.hdf5'
checkpoint = ModelCheckpoint(checkpoint_name, monitor='val_loss', verbose=1, save_best_only=True, mode='auto')
callbacks_list = [checkpoint]



# Learn CNN:
NN_model.fit(train, target, epochs=500, batch_size=32, validation_split=0.2, callbacks=callbacks_list)

