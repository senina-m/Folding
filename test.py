from Bio.PDB import *

pdbl = PDBList()
filename = '6gam'
pdbl.retrieve_pdb_file(filename, pdir='.', file_format='mmCif')
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure(filename, filename + '.cif')
model = structure.get_models()
models = list(model)
for i in range(len(models)):
    chains = list(models[i].get_chains())
    for j in range(len(chains)):
        residues = list(chains[j].get_residues())
        for k in range(len(residues)):
            atoms = list(residues[k].get_atoms())
            for l in range(len(atoms)):
                if atoms[l].get_fullname() == 'CA':
                    print(atoms[l].get_vector())
