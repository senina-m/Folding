from Bio.PDB import *

pdbl = PDBList()
filename = '6gov'
pdbl.retrieve_pdb_file(filename, pdir='.', file_format='mmCif')
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure(filename, filename + '.cif')
fullMatrix = [["ARG",
               [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "NE"], ["CG", "CD", "NE", "CZ"],
                ["CD", "NE", "CZ", "Nh1"]]],
              ["ASN", [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]]],
              ["ASP", [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]]],
              ["CYS", [["N", "CA", "CB", "SG"]]],
              ["GLN", [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "OE1"]]],
              ["GLU", [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "OE1"]]],
              ["HIS", [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]]],
              ["IlE", [["N", "CA", "CB", "CG1"], ["CA", "CB", "CG1", "CD"]]],
              ["LEU", [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]],
              ["LYS",
               [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "CE"], ["CG", "CD", "CE", "NZ"]]],
              ["MET", [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "SD", "CE"]]],
              ["PHE", [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]],
              ["PRO", [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"]]],
              ["SER", [["N", "CA", "CB", "CG"]]],
              ["THR", [["N", "CA", "CB", "OG1"]]],
              ["TRP", [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]],
              ["TYR", [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]]],
              ["VAL", [["N", "CA", "CB", "CG1"]]]]

result = [[[], [], [], [], []], [[], []], [[], []], [[]], [[], [], []], [[], [], []], [[], []], [[], []], [[], []],
          [[], [], [], []],
          [[], [], []], [[], []], [[], []], [[]], [[]], [[], []], [[], []], [[]]]

for residue in structure.get_residues():
    for i in 18:
        if residue.get_resname == fullMatrix[i][0]:
            for j in range(len(fullMatrix[i][1])):
                vector1 = residue[fullMatrix[i][1][j][0]].get_vector()
                vector2 = residue[fullMatrix[i][1][j][1]].get_vector()
                vector3 = residue[fullMatrix[i][1][j][2]].get_vector()
                vector4 = residue[fullMatrix[i][1][j][3]].get_vector()
                angle = calc_dihedral(vector1, vector2, vector3, vector4)
                result[i][j].append(angle)
