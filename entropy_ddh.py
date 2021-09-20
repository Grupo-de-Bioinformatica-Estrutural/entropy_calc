from Bio import PDB
import math
from subprocess import call as run_terminal
from itertools import chain
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse


def plot_entropy_in_dynamics(entropies, frames):

    plt.plot(frames, entropies)
    plt.xlabel("Number of frames")
    plt.ylabel("Entropy (kJ/mol)")
    plt.title("Entropy calculated by the DDH method")
    plt.savefig("../Entropies_DDH.png",dpi=300)
    plt.close()

def histogram(list_dihedrals, n, min_value, max_value):
    ''' Creates probability histogram of n bins from a list of encountered dihedrals from a simulation
        Input: list of dihedrals, number of bins, max and min value of possible dihedrals
        Output: Probability histrogram in n bins as a list'''



    bin_range = int ((abs(min_value) + abs(max_value)) / n)


    hist_bin_values = []        # List of each histogram bin
    count_hist = []             # List of counts to each bin


    for i in range(min_value, max_value, bin_range):

        hist_bin_values.append(i)
        count_hist.append(0)

    for dihedral in list_dihedrals:

        for i in range(0, len(hist_bin_values)):

            bin_value = hist_bin_values[i]

            try:

                upper_bin_value = hist_bin_values[i + 1]

            # If except, it means the dihedral belongs to the last bin
            except:

                count_hist[i] = count_hist[i] + 1
                continue


            # If dihedral is bigger or equal than this bin and lowe than next bin
            if dihedral >= bin_value and dihedral < upper_bin_value:
                count_hist[i] = count_hist[i] + 1                           # Add one count to this bin
                break



    # Probability histogram as a list
    prob_hist = []

    n_dihedrals = len(list_dihedrals)

    for count in count_hist:

        prob_hist.append(count/float(n_dihedrals))

    return prob_hist

def calc_chi_dihedrals(residue_name, residue):
    ''' Calculates all chi dihedrals possible from polypeptide residues
    Input: residue name, residue struct
    Output: List of chi values from respective aminoacid residue'''

    list_chi_angles = []

    if residue_name == "ARG":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG":

                CG = atom_coords

            elif atom_name == "CD":

                CD = atom_coords

            elif atom_name == "NE":

                NE = atom_coords

            elif atom_name == "CZ":

                CZ = atom_coords

            elif atom_name == "NH1":

                NH1 = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,CG))
        chi2 = math.degrees(PDB.calc_dihedral(CA,CB,CG,CD))
        chi3 = math.degrees(PDB.calc_dihedral(CB,CG,CD,NE))
        chi4 = math.degrees(PDB.calc_dihedral(CG,CD,NE,CZ))
        chi5 = math.degrees(PDB.calc_dihedral(CD,NE,CZ,NH1))

        list_chi_angles.append([chi1,chi2,chi3,chi4,chi5])

    elif residue_name == "ASN":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG":

                CG = atom_coords

            elif atom_name == "OD1":

                OD1 = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,CG))
        chi2 = math.degrees(PDB.calc_dihedral(CA,CB,CG,OD1))

        list_chi_angles.append([chi1,chi2])

    elif residue_name == "ASP":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG":

                CG = atom_coords

            elif atom_name == "OD1":

                OD1 = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,CG))
        chi2 = math.degrees(PDB.calc_dihedral(CA,CB,CG,OD1))

        list_chi_angles.append([chi1,chi2])

    elif residue_name == "CYS":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG":

                CG = atom_coords

            elif atom_name == "SG":

                SG = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,SG))

        list_chi_angles.append([chi1])

    elif residue_name == "GLN":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG":

                CG = atom_coords

            elif atom_name == "CD":

                CD = atom_coords

            elif atom_name == "OE1":

                OE1 = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,CG))
        chi2 = math.degrees(PDB.calc_dihedral(CA,CB,CG,CD))
        chi3 = math.degrees(PDB.calc_dihedral(CB,CG,CD,OE1))

        list_chi_angles.append([chi1,chi2,chi3])

    elif residue_name == "GLU":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG":

                CG = atom_coords

            elif atom_name == "CD":

                CD = atom_coords

            elif atom_name == "OE1":

                OE1 = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,CG))
        chi2 = math.degrees(PDB.calc_dihedral(CA,CB,CG,CD))
        chi3 = math.degrees(PDB.calc_dihedral(CB,CG,CD,OE1))

        list_chi_angles.append([chi1,chi2,chi3])

    elif residue_name == "HIS":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG":

                CG = atom_coords

            elif atom_name == "ND1":

                ND1 = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,CG))
        chi2 = math.degrees(PDB.calc_dihedral(CA,CB,CG,ND1))

        list_chi_angles.append([chi1,chi2])

    elif residue_name == "ILE":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG1":

                CG1 = atom_coords

            elif atom_name == "CD":

                CD = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,CG1))
        chi2 = math.degrees(PDB.calc_dihedral(CA,CB,CG1,CD))

        list_chi_angles.append([chi1,chi2])

    elif residue_name == "LEU":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG":

                CG = atom_coords

            elif atom_name == "CD1":

                CD1 = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,CG))
        chi2 = math.degrees(PDB.calc_dihedral(CA,CB,CG,CD1))

        list_chi_angles.append([chi1,chi2])

    elif residue_name == "LYS":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG":

                CG = atom_coords

            elif atom_name == "CE":

                CE = atom_coords

            elif atom_name == "CD":

                CD = atom_coords

            elif atom_name == "NZ":

                NZ = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,CG))
        chi2 = math.degrees(PDB.calc_dihedral(CA,CB,CG,CD))
        chi3 = math.degrees(PDB.calc_dihedral(CB,CG,CD,CE))
        chi4 = math.degrees(PDB.calc_dihedral(CG,CD,CE,NZ))

        list_chi_angles.append([chi1,chi2,chi3,chi4])

    elif residue_name == "MET":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG":

                CG = atom_coords

            elif atom_name == "SD":

                SD = atom_coords

            elif atom_name == "CE":

                CE = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,CG))
        chi2 = math.degrees(PDB.calc_dihedral(CA,CB,CG,SD))
        chi3 = math.degrees(PDB.calc_dihedral(CB,CG,SD,CE))

        list_chi_angles.append([chi1,chi2,chi3])

    elif residue_name == "PHE":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG":

                CG = atom_coords

            elif atom_name == "CD1":

                CD1 = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,CG))
        chi2 = math.degrees(PDB.calc_dihedral(CA,CB,CG,CD1))

        list_chi_angles.append([chi1,chi2])

    elif residue_name == "PRO":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG":

                CG = atom_coords

            elif atom_name == "CD":

                CD = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,CG))
        chi2 = math.degrees(PDB.calc_dihedral(CA,CB,CG,CD))

        list_chi_angles.append([chi1,chi2])

    elif residue_name == "SER":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "OG":

                OG = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,OG))

        list_chi_angles.append([chi1])

    elif residue_name == "THR":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "OG1":

                OG1 = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,OG1))

        list_chi_angles.append([chi1])

    elif residue_name == "TRP":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG":

                CG = atom_coords

            elif atom_name == "CD1":

                CD1 = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,CG))
        chi2 = math.degrees(PDB.calc_dihedral(CA,CB,CG,CD1))


        list_chi_angles.append([chi1,chi2])

    elif residue_name == "TYR":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG":

                CG = atom_coords

            elif atom_name == "CD1":

                CD1 = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,CG))
        chi2 = math.degrees(PDB.calc_dihedral(CA,CB,CG,CD1))

        list_chi_angles.append([chi1,chi2])

    elif residue_name == "VAL":
        for atom in residue:

            atom_name = atom.get_name()
            atom_coords = atom.get_vector()

            if atom_name == "N":

                N = atom_coords

            elif atom_name == "CA":

                CA = atom_coords

            elif atom_name == "CB":

                CB = atom_coords

            elif atom_name == "CG1":

                CG1 = atom_coords

        chi1 = math.degrees(PDB.calc_dihedral(N,CA,CB,CG1))

        list_chi_angles.append([chi1])



    list_chi_angles = list(chain(*list_chi_angles))

    return list_chi_angles

def probability_sum(prob_hist):
    ''' Somatory of probabilitys of each dihedral in simulation '''
    n = len(prob_hist)
    sum = 0

    for i in range(0,n):
        if prob_hist[i] != 0:
            sum += (prob_hist[i] * (math.log(prob_hist[i])))

    return sum * 8.3145

def calc_dihedrals(structure_name):
    ''' Calculates all polypeptide dihedrals to calculate entropy '''
    structure = PDB.PDBParser().get_structure(structure_name,structure_name)
    list_dihedrals = []


    # Get all phi psi dihedrals
    phi_psi_list = []

    for chain in structure:
        ppb = PDB.PPBuilder().build_peptides(chain)
        for poly in ppb:
            phi_psi_list.append(poly.get_phi_psi_list())


    for chain in phi_psi_list:
        for phi_psi in list(zip(*phi_psi_list)):
            for angle in list(zip(*phi_psi)):
                if angle[0] != None:
                    list_dihedrals.append(math.degrees(float(angle[0])))

    for model in structure:
        for chain in model:
            for residue in chain:
                residue_name = residue.get_resname()
                if residue_name == "ALA" and residue_name == "GLY":

                    continue

                else:

                    residue_chi_angles = calc_chi_dihedrals(residue_name,residue)

                for angle in residue_chi_angles:
                    list_dihedrals.append(angle)

    return list_dihedrals

def sum_indi_entropies(all_entropies):
    ''' Takes entropy from each dihedral and sums into final entropy
        Input: List of each dihedral entropy
        Output: Total entropy'''

    final_entropy = 0


    for entropy in all_entropies:

        final_entropy += entropy

    return final_entropy

def entropy_calculator(list_dihedrals, n, min_value, max_value):


    R = 8.3145
    prob_hist = histogram(list_dihedrals, n, min_value, max_value)

    sum_prob = probability_sum(prob_hist)

    return R/2 - (R * math.log(n)) - (sum_prob)

def main_protein(gmx_prefix, xtc_filename, tpr_filename, n, min_value, max_value):

    try:

        os.mkdir("frames/")

    except:

        raise Exception("A directory 'frames' already exist, please rm your current directory 'frames' or rename it")


    run_terminal("echo Protein Protein | {} trjconv -s {}.tpr -f {}.xtc -sep -skip 10 -o frames/frame_.pdb -pbc mol -center".format(gmx_prefix, xtc_filename, tpr_filename), shell = True)


    os.chdir("frames")
    all_frames = os.listdir(".")

    entropy_for_dihedral = []

    # Plot variables
    entropy_plot = []
    frames_plot = []



    # First case
    frame_dihedrals = calc_dihedrals(all_frames[0])
    all_dihedrals = [[None]] * len(frame_dihedrals)


    with open("../entropy_for_frames.txt","w") as results:

        for i,dihedral in enumerate(frame_dihedrals):
            all_dihedrals[i] = [dihedral]

        all_frames = all_frames[1:]

        for list_dihedrals in all_dihedrals:

            entropy_for_dihedral.append(entropy_calculator(list_dihedrals, n, min_value, max_value))



        # Takes first frame as first points
        skip = 0
        number_frames = 1

        for frame in all_frames:

            skip += 1
            number_frames += 1

            frame_dihedrals = calc_dihedrals(frame)

            for i,dihedral in enumerate(frame_dihedrals):

                all_dihedrals[i].append(dihedral)



            if skip == 500:

                skip = 0

                for i,list_dihedrals in enumerate(all_dihedrals):

                    entropy_for_dihedral[i] = entropy_calculator(list_dihedrals, n, min_value, max_value)


                final_entropy = sum_indi_entropies(entropy_for_dihedral)

                # Plot variables
                entropy_plot.append(final_entropy)
                frames_plot.append(number_frames)


        end_entropy = entropy_plot[-1]

        # Relative entropy
        for i,value in enumerate(entropy_plot):

            entropy_plot[i] = value - end_entropy

        for value_frame,value_entropy in zip(frames_plot,entropy_plot):

            results.write("{}\t{}\n".format(value_frame,value_entropy))

    plot_entropy_in_dynamics(entropy_plot, frames_plot)

def string_has_numbers(string):
    ''' Return true if string has at least one number, false otherwise '''
    return any(ch.isdigit() for ch in string)

def read_index(filename):
    ''' Reads a index of dihedrals of small molecules to make possible the calculation of entropy from these molecules
    And return chosen dihedral to calculate
    In the index file, every line is 4 atoms of a specific dihedral, example:
    1,5,4,3
    2,3,4,6
    '''

    dihedrals = []

    with open(filename,"r") as index:
        for line in index:
            if string_has_numbers(line):
                line_dihedral = line.split(",")
                # Strip possible spaces from string
                for i,atom_number in enumerate(line_dihedral):
                    line_dihedral[i] = int(atom_number)
                if len(line_dihedral) != 4:
                    raise Exception("Your index dihedral file {} has a line \
                                    that doesnt contain 4 atoms, but {} atoms! \
                                    \nFix this line:\n{}\
                                    ".format(filename,len(line_dihedral),line))
                dihedrals.append(line_dihedral)

            elif line != "\n":
                raise Exception("The line:\n{}Its not in the expected format!".format(line))

    return dihedrals

def calc_specific_dihedral(dihedral_atoms,structure):
    ''' Given a structure and index of atoms in dihedral, calculate this atoms index dihedral'''

    atom_index = 0

    dihedral_atom_vectors = []

    for residue in structure.get_residues():

        for atom in residue:

            atom_index += 1

            if atom_index in dihedral_atoms:

                dihedral_atom_vectors.append(atom.get_vector())

    atom1 = dihedral_atom_vectors[0]
    atom2 = dihedral_atom_vectors[1]
    atom3 = dihedral_atom_vectors[2]
    atom4 = dihedral_atom_vectors[3]

    return math.degrees(PDB.calc_dihedral(atom1, atom2, atom3, atom4))

def calc_dihedrals_from_index(dihedrals_atoms_index, struct_filename):
    ''' Calculates all dihedrals chosen in index from structure chosen'''

    all_dihedrals = []

    structure = PDB.PDBParser().get_structure(struct_filename,struct_filename)

    for list_atoms in dihedrals_atoms_index:

        all_dihedrals.append(calc_specific_dihedral(list_atoms, structure))

    return all_dihedrals

def main_small_molecule(gmx_prefix, tpr_filename, xtc_filename, index_filename , n, min_value, max_value):

    dihedrals_atoms_index = read_index(index_filename)

    # Plot variables
    entropy_plot = []
    frames_plot = []

    try:

        os.mkdir("frames/")

    except:

        raise Exception("A directory 'frames' already exist, please rm your current directory 'frames' or rename it")


    run_terminal("echo non-Water non-Water | {} trjconv -s {}.tpr -f {}.xtc -sep -skip 10 -o frames/frame_.pdb -pbc mol -center".format(gmx_prefix, xtc_filename, tpr_filename), shell = True)

    os.chdir("frames")
    all_frames = os.listdir(".")

    entropy_for_dihedral = []

    # Gets all chosen dihedrals

    # First case
    frame_dihedrals = calc_dihedrals_from_index(dihedrals_atoms_index,all_frames[0])
    all_dihedrals = [[None]] * len(frame_dihedrals)


    with open("../entropy_for_frames.txt","w") as results:

        for i,dihedral in enumerate(frame_dihedrals):
            all_dihedrals[i] = [dihedral]

        all_frames = all_frames[1:]

        for list_dihedrals in all_dihedrals:

            entropy_for_dihedral.append(entropy_calculator(list_dihedrals, n, min_value, max_value))



        skip = 0
        number_frames = 1

        for frame in all_frames:

            skip += 1
            number_frames += 1

            frame_dihedrals = calc_dihedrals_from_index(dihedrals_atoms_index, frame)

            for i,dihedral in enumerate(frame_dihedrals):

                all_dihedrals[i].append(dihedral)



            if skip == 500:

                skip = 0

                for i,list_dihedrals in enumerate(all_dihedrals):

                    entropy_for_dihedral[i] = entropy_calculator(list_dihedrals, n, min_value, max_value)


                final_entropy = sum_indi_entropies(entropy_for_dihedral)

                # Plot variables
                entropy_plot.append(final_entropy)
                frames_plot.append(number_frames)


        end_entropy = entropy_plot[-1]

        # Relative entropy
        for i,value in enumerate(entropy_plot):

            entropy_plot[i] = value - end_entropy

        for value_frame,value_entropy in zip(frames_plot,entropy_plot):

            results.write("{}\t{}\n".format(value_frame,value_entropy))

    plot_entropy_in_dynamics(entropy_plot, frames_plot)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-gmx",
        dest="gmx",
        help="gmx executable name. Default: gmx",
        default="gmx"
    )

    parser.add_argument(
        "-xtc",
        dest="xtc",
        help="XTC file from the dynamics."
    )

    parser.add_argument(
        "-tpr",
        dest="tpr",
        help="TPR file from the dynamics."
    )

    parser.add_argument(
        "-n",
        dest="n",
        help="Number of bins in entropy DDH method. Default: 72",
        default=72
    )

    parser.add_argument(
        "-i",
        dest="index",
        help="File with dihedrals atoms. See more info at the code repository.",
        default=72
    )

    parser.add_argument(
        "-min_value",
        dest="min_value",
        help="Minimum value of dihedral value to consider. Default: -180",
        default=-180
    )

    parser.add_argument(
        "-max_value",
        dest="max_value",
        help="Maximum value of dihedral value to consider. Default: 180",
        default=180
    )

    args = parser.parse_args()

    main_small_molecule(
        args.gmx,
        args.tpr,
        args.xtc,
        args.index,
        args.n,
        args.min_value,
        args.max_value
    )

#main_protein("gmx_51X","protein_min2","sgc1_1us-center", 72, -180, 180)
