from Bio.PDB import *
import numpy as np

pdbl = PDBList()
filename = '6gav'
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
                if atoms[l].get_fullname() == 'CA':
                    modelAtoms.append(atoms[l])

    matrix = []
    for j in range(len(modelAtoms)):
        matrix.append([0] * len(modelAtoms))

    for j in range(len(modelAtoms)):
        for k in range(len(modelAtoms)):
            matrix[j][k] = modelAtoms[j] - modelAtoms[k]
            print(matrix[j][k], end=' ')
        print()
