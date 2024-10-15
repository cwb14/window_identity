#!/usr/bin/env python3

# Run like:
# python calculate_k2p.py <(awk '{for(i=1;i<=NF;i++) if($i ~ /^cg:Z:/){c=substr($i,6);s=0; while(match(c,/([0-9]+)M/,a)){s+=a[1];c=substr(c,RSTART+RLENGTH)}; print $0 "\t" s}}' ./Nipp_to_9311.paf) > Nipp_to_9311_adjust.paf
# The awk line counts the number of Ms in the CIGAR (matches & mismatches) for use as total_aligned_bases in the K2P calculation.

import argparse
import math
import sys
import re

def parse_cs(cs_str):
    """
    Parse the cs:Z: string from PAF to extract operations.
    Returns a list of tuples: (op, value)
    where op is one of ':', '*', '+', '-'
    and value is the associated string or number.
    """
    pattern = re.compile(r'([:*+-])([^:*+-]+)')
    matches = pattern.findall(cs_str)
    return matches

def is_transition(ref, alt):
    """Check if the substitution is a transition."""
    transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
    return (ref, alt) in transitions

def calculate_k2p(transitions, transversions, total_sites):
    """
    Calculate the Kimura 2-Parameter distance.
    K2P = -0.5 * ln(1 - 2p - q) - 0.25 * ln(1 - 2q)
    where p = transitions / total_sites
          q = transversions / total_sites
    """
    if total_sites == 0:
        return 'NA'
    p = transitions / total_sites
    q = transversions / total_sites
    try:
        k2p = -0.5 * math.log(1 - 2*p - q) - 0.25 * math.log(1 - 2*q)
        return f"{k2p:.6f}"
    except ValueError:
        return 'NA'  # Handle math domain error if arguments are out of range

def process_paf_line(line, exclude_adjacent, exclude_near_indel, n_bp):
    """
    Process a single PAF line to calculate K2P.
    """
    fields = line.strip().split('\t')
    # Extract total aligned bases from the last column
    try:
        total_aligned_bases = int(fields[-1])
    except ValueError:
        # If the last column is not an integer, search for it
        total_aligned_bases = int(fields[-1].split()[0])

    # Extract the cs:Z: field
    cs_field = None
    for field in fields[12:]:
        if field.startswith('cs:Z:'):
            cs_field = field[5:]
            break
    if cs_field is None:
        # If no cs:Z: field, cannot proceed
        return line.strip() + '\tNA'

    # Parse the cs:Z: field
    operations = parse_cs(cs_field)

    substitutions = []
    indels = []
    current_pos = 0  # Position in the aligned sequence

    for op, value in operations:
        if op == ':':
            # Match, advance position
            match_length = int(value)
            current_pos += match_length
        elif op == '*':
            # Substitution, advance by 1
            if len(value) == 2:
                ref, alt = value[0], value[1]
                substitutions.append((current_pos, ref.upper(), alt.upper()))
                current_pos += 1
            else:
                # Handle unexpected substitution format
                current_pos += len(value)
        elif op in ('+', '-'):
            # Indel, record position and length
            indel_length = len(value)
            indels.append(current_pos)
            # Indels do not consume positions in the aligned sequence
        else:
            # Unknown operation, skip
            pass

    # Apply filtering
    excluded = set()
    excluded_d_count = 0
    excluded_r_count = 0

    if exclude_adjacent:
        # Sort substitutions by position
        subs_sorted = sorted(substitutions, key=lambda x: x[0])
        for i in range(1, len(subs_sorted)):
            prev = subs_sorted[i-1]
            current = subs_sorted[i]
            if current[0] == prev[0] + 1:
                excluded.add(i-1)
                excluded.add(i)
                excluded_d_count += 1

    if exclude_near_indel:
        # For each substitution, check distance to nearest indel
        indel_positions = sorted(indels)
        for i, sub in enumerate(substitutions):
            sub_pos = sub[0]
            # Find the closest indel position
            if not indel_positions:
                break
            # Binary search could be used here for efficiency
            min_distance = min(abs(sub_pos - indel) for indel in indel_positions)
            if min_distance <= n_bp:
                excluded.add(i)
                excluded_r_count += 1

    # Count transitions and transversions
    transitions_count = 0
    transversions_count = 0
    excluded_count = len(excluded)

    for i, (pos, ref, alt) in enumerate(substitutions):
        if i in excluded:
            continue
        if is_transition(ref, alt):
            transitions_count += 1
        else:
            transversions_count += 1

    # Adjust total aligned bases
    adjusted_total = total_aligned_bases - excluded_count # Its a bit ambiguous to me if the filted variants should be excluded from total_aligned_bases (due to being gaps) or included, due to being mismatch. Here, theyre excluded. The '-d' and '-r' filters are experimental and require further analysis know. I see a case that '-d' are indels and '-r' are suspisious mismatches we're ignoring and so both should be excluded from total_aligned_bases.

    # Calculate K2P
    k2p = calculate_k2p(transitions_count, transversions_count, adjusted_total)

    # Debugging output
#    print(f"Debug info for PAF line:")
#    print(f"  Total transitions: {transitions_count}")
#    print(f"  Total transversions: {transversions_count}")
#    print(f"  Excluded with '-d' (adjacent variants): {excluded_d_count}")
#    print(f"  Excluded with '-r' (near indels): {excluded_r_count}")
#    print(f"  Total aligned bases used for K2P: {adjusted_total}\n")

    return line.strip() + f'\t{k2p}'

def main():
    parser = argparse.ArgumentParser(description="Calculate Kimura 2-Parameter (K2P) distance for PAF alignments with debugging info.")
    parser.add_argument('-d', action='store_true', help='Exclude adjacent variants.')
    parser.add_argument('-r', type=int, metavar='N', help='Exclude variants within N bp of indels.')
    parser.add_argument('paf_file', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Input PAF file (default: stdin)')
    args = parser.parse_args()

    exclude_adjacent = args.d
    exclude_near_indel = False
    n_bp = 0
    if args.r is not None:
        exclude_near_indel = True
        n_bp = args.r

    for line in args.paf_file:
        if line.startswith('#') or line.strip() == '':
            continue  # Skip headers or empty lines
        k2p_line = process_paf_line(line, exclude_adjacent, exclude_near_indel, n_bp)
        print(k2p_line)

if __name__ == "__main__":
    main()
