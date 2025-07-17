#!/usr/bin/env python3
import argparse
import re
import sys

# Not sure if strand flipping is required... i need to look into this.

def parse_cg(cg_str):
    """
    Parse the cg:Z: string into a list of (length, op) tuples.
    For example: "359M1D293M1D56M6I326M..."
    Returns a list of (int, op) tuples.
    """
    pattern = re.compile(r'(\d+)([MID])')
    return [(int(num), op) for num, op in pattern.findall(cg_str)]

def process_paf_line(fields, cutoff):
    """
    Process one PAF alignment line.
    Returns a list of events, where each event is a tuple:
    (indel_type, query_name, q_indel_start, q_indel_end, target_name, t_indel_start, t_indel_end, indel_length)
    
    Coordinates are computed exactly based on the cg operations.
    Both query and target intervals will have the same span equal to the indel length.
    Only events with length >= cutoff are returned.
    
    For minus strand alignments (column 5 == '-'), the cg operations are read in reverse order
    starting from query_end/target_end and subtracting lengths.
    """
    # Required fields from PAF:
    query_name   = fields[0]
    query_start  = int(fields[2])
    query_end    = int(fields[3])
    strand       = fields[4]
    target_name  = fields[5]
    target_start = int(fields[7])
    target_end   = int(fields[8])
    
    # Find the cg:Z: field (optional fields start at index 12)
    cg_field = None
    for field in fields[12:]:
        if field.startswith("cg:Z:"):
            cg_field = field[5:]
            break
    if cg_field is None:
        return []

    # Set up starting positions and operation order based on strand.
    if strand == '-':
        # For minus strand, start from the alignment end and process operations in reverse.
        qpos = query_end
        tpos = target_end
        ops = list(reversed(parse_cg(cg_field)))
    else:
        qpos = query_start
        tpos = target_start
        ops = parse_cg(cg_field)

    events = []
    # Process each operation in order
    for length, op in ops:
        if op == 'M':
            if strand == '-':
                qpos -= length
                tpos -= length
            else:
                qpos += length
                tpos += length
        elif op in ['I', 'D']:
            if cutoff <= length:
                if strand == '-':
                    # For minus strand, the indel spans from (current position - length) to current position.
                    event = (op, query_name, qpos - length, qpos,
                             target_name, tpos - length, tpos, length)
                else:
                    event = (op, query_name, qpos, qpos + length,
                             target_name, tpos, tpos + length, length)
                events.append(event)
            if strand == '-':
                qpos -= length
                tpos -= length
            else:
                qpos += length
                tpos += length
        else:
            if strand == '-':
                qpos -= length
                tpos -= length
            else:
                qpos += length
                tpos += length

    return events

def main():
    parser = argparse.ArgumentParser(
        description="Extract exact indel events from a PAF file based on cg:Z: field."
    )
    parser.add_argument("-in", dest="infile", required=True, help="Input PAF file")
    parser.add_argument("-out", dest="outfile", required=True, help="Output TSV file")
    parser.add_argument("-cutoff", dest="cutoff", type=int, default=0,
                        help="Minimum indel length to output (default: 0)")
    
    # If no arguments are provided, print help and exit.
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()

    all_events = []
    try:
        with open(args.infile, 'r') as infile:
            for line in infile:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    continue
                fields = line.split("\t")
                events = process_paf_line(fields, args.cutoff)
                all_events.extend(events)
    except Exception as e:
        sys.exit(f"Error reading input file: {e}")

    try:
        with open(args.outfile, 'w') as outfile:
            for evt in all_events:
                # evt is a tuple: (op, qname, qstart, qend, tname, tstart, tend, length)
                outfile.write("\t".join(str(x) for x in evt) + "\n")
    except Exception as e:
        sys.exit(f"Error writing output file: {e}")

if __name__ == "__main__":
    main()
