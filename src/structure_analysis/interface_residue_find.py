#!/usr/bin/env python
"""
This script reads a PDB file containing a protein complex (chains A and B),
identifies interface interactions (such as salt bridges, hydrogen bonds, and hydrophobic contacts)
by searching for specific interaction types rather than simply using a distance cutoff.
It outputs the protein sequences (one-letter code) with residue numbers annotated and 
generates a circular plot on a single circle where the sequences for both chains are concatenated.
A clear, small gap is inserted both between Chain A and Chain B 
and between Chain B and Chain A so that the two chains are clearly separated.
Residues are represented as follows:
   - Residues from Chain A are labeled with an "A" (and come first) and those from
     Chain B with a "B".
   - Each amino acid is colored based on its type:
        • Positive (K, R, H) – blue
        • Negative (D, E) – red
        • Hydrophobic (A, V, L, I, M, F, W, Y) – green
        • Polar (S, T, N, Q) – orange
        • Special (C, G, P) – purple
        • Others – black
   - Residue numbers (with chain ID) are annotated if it is the first residue,
     the last residue, divisible by 5, or if it is an interface residue.
   - Connecting lines are drawn between interacting residue pairs using
     different styles/colors according to the interaction type:
         • Salt bridges: solid red lines.
         • Hydrogen bonds: solid turquoise lines.
         • Hydrophobic interactions: solid gold lines.
     
Usage:
    python "20250123_FKBP12/6_interface_residue_find.py" <pdb_file> [-c <cutoff>] [-a <chainA>] [-b <chainB>] [-first <first_chain_name>] [-second <second_chain_name>]

Example:
    python "20250123_FKBP12/6_interface_residue_find.py" 1a12.pdb -c 5.0 -first "FKBP12" -second "Binder"
"""

from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser, NeighborSearch, Polypeptide
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use("nature")


#########################
# Interaction Functions #
#########################

def is_salt_bridge(resA, resB, cutoff=4.0):
    """
    Determine if resA and resB form a salt bridge.
    A salt bridge is defined here as an interaction between an acidic residue (D,E)
    and a basic residue (K,R,H) in which at least one acidic side-chain oxygen (OD1,OD2,OE1,OE2)
    is within the specified cutoff (default 4.0 Å) of one basic side-chain nitrogen.
    """
    acidic_set = set("DE")
    basic_set = set("KRH")
    resA_type = seq1(resA.resname)
    resB_type = seq1(resB.resname)
    # Check that one residue is acidic and the other is basic.
    if (resA_type in acidic_set and resB_type in basic_set) or (resA_type in basic_set and resB_type in acidic_set):
        acidic_atoms = []
        basic_atoms = []
        # Determine which is acidic.
        if resA_type in acidic_set:
            for atom in resA:
                if atom.get_name() in ['OD1', 'OD2', 'OE1', 'OE2']:
                    acidic_atoms.append(atom)
            for atom in resB:
                if resB_type == "K" and atom.get_name() == "NZ":
                    basic_atoms.append(atom)
                elif resB_type == "R" and atom.get_name() in ["NH1", "NH2", "NE"]:
                    basic_atoms.append(atom)
                elif resB_type == "H" and atom.get_name() in ["ND1", "NE2"]:
                    basic_atoms.append(atom)
        else:
            for atom in resB:
                if atom.get_name() in ['OD1', 'OD2', 'OE1', 'OE2']:
                    acidic_atoms.append(atom)
            for atom in resA:
                if resA_type == "K" and atom.get_name() == "NZ":
                    basic_atoms.append(atom)
                elif resA_type == "R" and atom.get_name() in ["NH1", "NH2", "NE"]:
                    basic_atoms.append(atom)
                elif resA_type == "H" and atom.get_name() in ["ND1", "NE2"]:
                    basic_atoms.append(atom)
        min_distance = None
        for a_atom in acidic_atoms:
            for b_atom in basic_atoms:
                distance = a_atom - b_atom
                if min_distance is None or distance < min_distance:
                    min_distance = distance
        if min_distance is not None and min_distance <= cutoff:
            return True, min_distance
    return False, None


def is_hydrophobic_interaction(resA, resB, cutoff=5.0):
    """
    Determine if resA and resB form a hydrophobic interaction.
    Here the interaction is defined between two hydrophobic residues (A, V, L, I, M, F, W, Y)
    if any pair of non-backbone heavy atoms are within the cutoff (default 5.0 Å).
    """
    hydrophobic_set = set("AVLIMFWY")
    resA_type = seq1(resA.resname)
    resB_type = seq1(resB.resname)
    if resA_type in hydrophobic_set and resB_type in hydrophobic_set:
        min_distance = None
        for a_atom in resA:
            if a_atom.get_name() in ['N', 'CA', 'C', 'O']:
                continue
            for b_atom in resB:
                if b_atom.get_name() in ['N', 'CA', 'C', 'O']:
                    continue
                distance = a_atom - b_atom
                if min_distance is None or distance < min_distance:
                    min_distance = distance
        if min_distance is not None and min_distance <= cutoff:
            return True, min_distance
    return False, None


def is_hydrogen_bond(resA, resB, cutoff=3.5):
    """
    Determine if resA and resB form a hydrogen bond.
    This simple implementation uses a distance-based criterion.
    We designate candidate donor atoms and acceptor atoms, and if any donor from one residue is
    within the cutoff distance of any acceptor from the other residue, a hydrogen bond is assumed.

    Candidate donor atoms: ["N", "NE", "ND", "NZ", "NH1", "NH2", "ND1", "NE2", "OG", "OH"]
    Candidate acceptor atoms: ["O", "OD1", "OD2", "OE1", "OE2", "OXT"]
    """
    donor_atoms = ["N", "NE", "ND", "NZ",
                   "NH1", "NH2", "ND1", "NE2", "OG", "OH"]
    acceptor_atoms = ["O", "OD1", "OD2", "OE1", "OE2", "OXT"]

    # Check resA's donor atoms to resB's acceptor atoms.
    for a_atom in resA:
        if a_atom.get_name() in donor_atoms:
            for b_atom in resB:
                if b_atom.get_name() in acceptor_atoms:
                    distance = a_atom - b_atom
                    if distance <= cutoff:
                        return True, distance
    # Check resB's donor atoms to resA's acceptor atoms.
    for b_atom in resB:
        if b_atom.get_name() in donor_atoms:
            for a_atom in resA:
                if a_atom.get_name() in acceptor_atoms:
                    distance = b_atom - a_atom
                    if distance <= cutoff:
                        return True, distance
    return False, None


def get_interaction(resA, resB):
    """
    Check for possible specific interactions between resA and resB.
    Returns a tuple (interaction_type, minimal_distance) if an interaction is found.
    Priority order: salt_bridge > hydrogen_bond > hydrophobic.
    """
    is_sb, d_sb = is_salt_bridge(resA, resB)
    if is_sb:
        return "salt_bridge", d_sb
    is_hbond, d_hbond = is_hydrogen_bond(resA, resB)
    if is_hbond:
        return "hydrogen_bond", d_hbond
    is_hydro, d_hydro = is_hydrophobic_interaction(resA, resB)
    if is_hydro:
        return "hydrophobic", d_hydro
    return None, None


def get_interactions(chain_a, chain_b):
    """
    Iterate over residue pairs from chain_a and chain_b and search for specific interactions.
    Returns a dictionary with keys as tuples (resnum_A, resnum_B) and values as (interaction_type, distance).
    """
    interactions = {}
    for res_a in chain_a.get_residues():
        if not Polypeptide.is_aa(res_a, standard=True):
            continue
        for res_b in chain_b.get_residues():
            if not Polypeptide.is_aa(res_b, standard=True):
                continue
            inter_type, inter_dist = get_interaction(res_a, res_b)
            if inter_type is not None:
                resnum_A = res_a.get_id()[1]
                resnum_B = res_b.get_id()[1]
                key = (resnum_A, resnum_B)
                # Save the interaction; if multiple contacts exist, keep the one with the shortest distance.
                if key not in interactions or (inter_dist is not None and inter_dist < interactions[key][1]):
                    interactions[key] = (inter_type, inter_dist)
    return interactions


def draw_curly_line(ax, xA, yA, xB, yB, linestyle, line_color):
    """
    Draw a curved (Bézier) line between points (xA, yA) and (xB, yB) using
    a quadratic Bézier curve with a control point offset by a perpendicular factor.
    """
    # Calculate the midpoint.
    xm = (xA + xB) / 2.0
    ym = (yA + yB) / 2.0
    # Compute directional vector from A to B.
    dx = xB - xA
    dy = yB - yA
    # Calculate a perpendicular vector.
    perp = np.array([-dy, dx])
    norm = np.linalg.norm(perp)
    if norm != 0:
        perp_unit = perp / norm
    else:
        perp_unit = np.array([0, 0])

    # Determine how far the control point should be from the midpoint.
    distance = np.sqrt(dx**2 + dy**2)
    # Adjust this value to change the curvature (higher value = curlier)
    curvature = 0.15
    offset = curvature * distance

    # Define the control point by offsetting the midpoint perpendicularly.
    control_x = xm + offset * perp_unit[0]
    control_y = ym + offset * perp_unit[1]

    # Generate the Bézier curve (parametric calculation)
    t = np.linspace(0, 1, 100)
    bezier_x = (1-t)**2 * xA + 2*(1-t)*t * control_x + t**2 * xB
    bezier_y = (1-t)**2 * yA + 2*(1-t)*t * control_y + t**2 * yB
    ax.plot(bezier_x, bezier_y, linestyle=linestyle,
            color=line_color, linewidth=2, alpha=1)


#########################
# Main Plotting Section #
#########################


def main():
    # Parse command-line arguments.
    parser = argparse.ArgumentParser(
        description="Identify specific interface interactions (salt bridges, hydrogen bonds, hydrophobic contacts) between two chains, and plot a circular diagram "
                    "with concatenated sequences (with a small gap separating the chains on both sides) plus connecting interaction lines."
    )
    parser.add_argument("pdb_file", type=str,
                        help="Input PDB file with protein complex.")
    parser.add_argument("-c", "--cutoff", type=float, default=5.0,
                        help="Distance cutoff (in Å) for hydrophobic interaction determination (default: 5.0 Å). Note: Salt bridge and hydrogen bond cutoffs are fixed within their functions.")
    # New arguments for specifying chain names/identifiers.
    parser.add_argument("-a", "--chainA", type=str, default="A",
                        help="Chain identifier for Chain A in the PDB file (default: A)")
    parser.add_argument("-b", "--chainB", type=str, default="B",
                        help="Chain identifier for Chain B in the PDB file (default: B)")
    parser.add_argument("-first", "--first_chain_name", type=str, default="",
                        help="Name of first chain (default: empty string)")
    parser.add_argument("-second", "--second_chain_name", type=str, default="",
                        help="Name of second chain (default: empty string)")
    args = parser.parse_args()

    pdb_file = args.pdb_file
    cutoff = args.cutoff

    # Parse the PDB structure (using the first model).
    parser_pdb = PDBParser(QUIET=True)
    structure = parser_pdb.get_structure("protein", pdb_file)
    model = structure[0]

    # Retrieve specified chains.
    try:
        chain_a = model[args.chainA]
    except KeyError:
        sys.exit(f"Error: Chain {args.chainA} not found in the PDB file.")
    try:
        chain_b = model[args.chainB]
    except KeyError:
        sys.exit(f"Error: Chain {args.chainB} not found in the PDB file.")

    # Build NeighborSearch objects for interface detection.
    chain_b_atoms = list(chain_b.get_atoms())
    ns_b = NeighborSearch(chain_b_atoms)
    chain_a_atoms = list(chain_a.get_atoms())
    ns_a = NeighborSearch(chain_a_atoms)

    # Identify interface residues for Chain A.
    interface_a = []
    for res in chain_a.get_residues():
        if not Polypeptide.is_aa(res, standard=True):
            continue
        for atom in res:
            if ns_b.search(atom.coord, cutoff):
                interface_a.append(res)
                break

    # Identify interface residues for Chain B.
    interface_b = []
    for res in chain_b.get_residues():
        if not Polypeptide.is_aa(res, standard=True):
            continue
        for atom in res:
            if ns_a.search(atom.coord, cutoff):
                interface_b.append(res)
                break

    # Build sequence for Chain A.
    interface_ids_a = {res.get_id() for res in interface_a}
    chainA_letters = []       # one-letter codes for Chain A
    # annotated residue numbers (as strings, with chain ID)
    chainA_resnum_list = []
    # mapping: residue number (int) -> local index for Chain A
    chainA_map = {}
    # local indices for interface residues (Chain A)
    interface_positions_A = []
    index = 0
    for res in chain_a:
        if not Polypeptide.is_aa(res, standard=True):
            continue
        try:
            letter = seq1(res.resname)
        except Exception:
            letter = "X"
        chainA_letters.append(letter)
        resnum = res.get_id()[1]
        chainA_resnum_list.append("" + str(resnum))
        chainA_map[resnum] = index
        if res.get_id() in interface_ids_a:
            interface_positions_A.append(index)
        index += 1
    n_A = len(chainA_letters)

    # Build sequence for Chain B.
    interface_ids_b = {res.get_id() for res in interface_b}
    chainB_letters = []       # one-letter codes for Chain B
    # annotated residue numbers (as strings, with chain ID)
    chainB_resnum_list = []
    # mapping: residue number (int) -> local index for Chain B
    chainB_map = {}
    # local indices for interface residues (Chain B)
    interface_positions_B_local = []
    index = 0
    for res in chain_b:
        if not Polypeptide.is_aa(res, standard=True):
            continue
        try:
            letter = seq1(res.resname)
        except Exception:
            letter = "X"
        chainB_letters.append(letter)
        resnum = res.get_id()[1]
        chainB_resnum_list.append("" + str(resnum))
        chainB_map[resnum] = index
        if res.get_id() in interface_ids_b:
            interface_positions_B_local.append(index)
        index += 1
    n_B = len(chainB_letters)

    # (Optional) Print individual sequences and interface positions.
    print("Chain A sequence:")
    print("".join(chainA_letters))
    print("Interface residue positions in Chain A (local index):",
          interface_positions_A)
    print("\nChain B sequence:")
    print("".join(chainB_letters))
    print("Interface residue positions in Chain B (local index):",
          interface_positions_B_local)

    # Instead of using a simple distance cutoff, search for specific interaction types.
    interactions = get_interactions(chain_a, chain_b)
    if interactions:
        print("\nInteraction pairs (residue numbers, interaction type, minimal distance):")
        for key, (inter_type, dist) in sorted(interactions.items(), key=lambda x: x[0][0]):
            print(
                f"Chain A residue {key[0]} interacts with Chain B residue {key[1]} via {inter_type} (distance {dist:.2f} Å)")
    else:
        print("\nNo specific interactions detected.")

    # ---------------------------
    # Circular Plot: Single Circle with Two Gaps (one between chains, and one between chains)
    # ---------------------------
    R = 1.0           # radius of the circle
    gap = 0.1         # angular gap in radians
    total_angle = 2 * np.pi
    total_gap = 2 * gap  # one gap between Chain A and Chain B in each direction
    available_arc = total_angle - total_gap

    # Allocate arcs for each chain proportionally to their number of residues.
    arc_A = available_arc * (n_A / (n_A + n_B))
    arc_B = available_arc * (n_B / (n_A + n_B))

    # Distribute the gaps symmetrically:
    # Center Chain A's arc at 3π/2, so subtract half of arc_A.
    start_A = 3 * np.pi / 2 - arc_A / 2
    angles_A = np.linspace(start_A, start_A + arc_A, n_A, endpoint=False)

    # Center Chain B's arc at π/2, so subtract half of arc_B.
    start_B = np.pi / 2 - arc_B / 2
    angles_B = np.linspace(start_B, start_B + arc_B, n_B, endpoint=False)

    # Create figure.
    fig, ax = plt.subplots(figsize=(12, 12))
    ax.set_aspect('equal')
    ax.axis('off')

    # Helper function for residue color based on amino acid.
    def get_color(letter):
        if letter in "KRH":
            return 'blue'
        elif letter in "DE":
            return 'red'
        elif letter in "AVLIMFWY":
            return 'green'
        elif letter in "STNQ":
            return 'orange'
        elif letter in "CGP":
            return 'purple'
        else:
            return 'black'

    # Plot Chain A residues.
    for i, letter in enumerate(chainA_letters):
        theta = angles_A[i]
        x = R * np.cos(theta)
        y = R * np.sin(theta)
        if i in interface_positions_A:
            fontweight = 'bold'
            fontsize = 12
        else:
            fontweight = 'normal'
            fontsize = 6
        ax.text(x, y, letter, fontsize=fontsize, color=get_color(letter),
                fontweight=fontweight, ha='center', va='center')
        # Annotate residue numbers if first, last, or every 5th.
        if (i == 0) or (i == n_A - 1) or ((i + 1) % 5 == 0):
            x_num = (R + 0.05) * np.cos(theta)
            y_num = (R + 0.05) * np.sin(theta)
            ax.text(x_num, y_num, chainA_resnum_list[i], fontsize=12,
                    color='black', ha='center', va='center')

    # Plot Chain B residues.
    for i, letter in enumerate(chainB_letters):
        theta = angles_B[i]
        x = R * np.cos(theta)
        y = R * np.sin(theta)
        if i in interface_positions_B_local:
            fontweight = 'bold'
            fontsize = 12
        else:
            fontweight = 'normal'
            fontsize = 6
        ax.text(x, y, letter, fontsize=fontsize, color=get_color(letter),
                fontweight=fontweight, ha='center', va='center')
        # Annotate residue numbers if first, last, or every 5th.
        if (i == 0) or (i == n_B - 1) or ((i + 1) % 5 == 0):
            x_num = (R + 0.05) * np.cos(theta)
            y_num = (R + 0.05) * np.sin(theta)
            ax.text(x_num, y_num, chainB_resnum_list[i], fontsize=12,
                    color='black', ha='center', va='center')

    # Add labels for each chain.
    mid_A = start_A + arc_A / 2
    mid_B = start_B + arc_B / 2
    # You can customize these labels or even add new command-line arguments for display names.
    ax.text((R + 0.1) * np.cos(mid_A), (R + 0.1) * np.sin(mid_A), f"{args.first_chain_name}",
            fontsize=14, color='black', ha='center', va='center', fontweight='bold')
    ax.text((R + 0.1) * np.cos(mid_B), (R + 0.1) * np.sin(mid_B), f"{args.second_chain_name}",
            fontsize=14, color='black', ha='center', va='center', fontweight='bold')

    # Draw connecting lines between interacting residues.
    for (resnum_A, resnum_B), (inter_type, dist) in interactions.items():
        if (resnum_A in chainA_map) and (resnum_B in chainB_map):
            i_A = chainA_map[resnum_A]   # Chain A local index
            i_B = chainB_map[resnum_B]   # Chain B local index
            theta_A = angles_A[i_A]
            theta_B = angles_B[i_B]
            xA = (R - 0.02) * np.cos(theta_A)
            yA = (R - 0.02) * np.sin(theta_A)
            xB = (R - 0.02) * np.cos(theta_B)
            yB = (R - 0.02) * np.sin(theta_B)
            if inter_type == "salt_bridge":
                linestyle = '-'
                line_color = '#E94560'  # red
            elif inter_type == "hydrogen_bond":
                linestyle = '-'
                line_color = '#4EBBC1'  # turquoise
            elif inter_type == "hydrophobic":
                linestyle = '-'
                line_color = '#D9AA28'  # gold
            else:
                linestyle = '-'
                line_color = get_color(chainB_letters[i_B])
            draw_curly_line(ax, xA, yA, xB, yB, linestyle, line_color)

    image_file = f"{pdb_file.split('.')[0]}.png"
    plt.savefig(image_file, dpi=300, bbox_inches='tight')
    print(f"\nCircular sequence plot saved as {image_file}")
    plt.show()


if __name__ == "__main__":
    main()
