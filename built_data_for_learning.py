import time
from Bio.PDB import *
from Bio import SeqIO
import numpy as np
import sympy
import os
import shutil
import tensorflow as tf
import pandas as pd


# def get_data():
#     # get train data
#     train_data_path = 'train.csv'
#     train = pd.read_csv(train_data_path)
#
#     # get test data
#     test_data_path = 'test.csv'
#     test = pd.read_csv(test_data_path)
#     return train, test

# TODO: split dataset into training and testing parts
def split_combined(combined):
    train = combined[:1460]
    test = combined[1460:]

    return train, test


# combined = tf.data.Dataset.from_tensor_slices(tf.random.uniform([4, 10]))
# train, test = split_combined(combined)


def get_matrix(filename):
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(filename, pdir='.', file_format='mmCif')
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(filename, filename + '.cif')
    model = structure.get_models()
    models = list(model)
    chains = list(models[1].get_chains())
    modelatoms = list()
    for j in range(len(chains)):
        residues = list(chains[j].get_residues())
        for k in range(len(residues)):
            atoms = list(residues[k].get_atoms())
            for l in range(len(atoms)):
                if atoms[l].get_fullname() == 'CA':
                    modelatoms.append(atoms[l])
    matrix = []
    for j in range(len(modelatoms)):
        matrix.append([0] * len(modelatoms))
    for j in range(len(modelatoms)):
        for k in range(len(modelatoms)):
            matrix[j][k] = (modelatoms[j] - modelatoms[k]) * (modelatoms[j] - modelatoms[k])
    matrix_1 = np.asarray(matrix)
    print(matrix_1.shape)
    start_time = time.time()
    _, indexes = sympy.Matrix(matrix_1).T.rref()
    print(time.time() - start_time)
    m = []
    for ind in indexes:
        m.append(matrix_1[ind])
    return m


def get_training_data():
    path2base = "C:/Практика Biocad/PDB/resultbase"
    path2python = "C:/Практика Biocad/Python Projects"
    number_of_proteins = 100
    listdir = os.listdir(path2base)
    for protein_num in range(0, number_of_proteins):
        protein_name = listdir[protein_num]
        shutil.copyfile(path2base + "/" + protein_name, path2python + "/" + protein_name)
        m = get_matrix(protein_name)
        os.remove(path2python + "/" + protein_name)
        data = np.array([0, 0])
        data[0] = get_fasta_onehot(get_fasta_onehot(protein_name))
        data[1] = m
        dataset = tf.data.Dataset.from_tensor_slices(data)
    return dataset


def get_fasta_onehot(filename):
    filename += ".fasta.txt"
    list_of_fasta_arrs = []
    for seq_record in SeqIO.parse(filename, "fasta"):
        list_of_fasta_arrs.append(list(seq_record.seq[:]))
        break
    # print(list_of_fasta_arrs)
    vocab = ['A', ' B', ' C', ' D', ' E', ' F', ' G', ' H', ' I', ' J', ' K', ' L', ' M', ' N', ' O', ' P', ' Q', ' R',
             ' S', ' T', ' U', ' V', ' W', ' X']
    t = tf.compat.v1.placeholder(dtype=tf.string, shape=(None,))
    matches = tf.stack([tf.equal(t, s) for s in vocab], axis=-1)
    onehot = tf.cast(matches, tf.float32)

    # TODO: make a loop and change the indexes in list_of_fasta_arrs
    with tf.compat.v1.Session() as sess:
        onehot_encoded_sequence = sess.run(onehot, feed_dict={t: list_of_fasta_arrs[0]})
    ONE_HOT_ENCODING_LEN = len(onehot_encoded_sequence[0])
    return onehot_encoded_sequence
