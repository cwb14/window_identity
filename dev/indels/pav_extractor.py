#!/usr/bin/env python3
import argparse
import re

def parse_cs(cs_str, qstart, tstart):
    """
    Parse a cs:Z: string and yield PAV events.
    
    For each operation, update query and target offsets.
    For an insertion (prefixed by '+'), record an event using the current query offset.
    For a deletion (prefixed by '-'), record an event using the current target offset.
    
    Yields dictionaries:
      {
         'type': 'ins' or 'del',
         'seq': the inserted/deleted sequence,
         'pav_start': computed coordinate (1-based),
         'pav_end': computed coordinate (1-based)
      }
      
    Note: This assumes qstart and tstart are the starting positions from the PAF file.
    """
    # Patterns for operations:
    # Match: :<number>
    # Insertion: +<letters>
    # Deletion: -<letters>
    # Substitution: *<letter><letter>
    pattern = re.compile(r'(?P<match>:(\d+))|(?P<ins>\+([a-zA-Z]+))|(?P<del>-(?:[a-zA-Z]+))|(?P<sub>\*[a-zA-Z]{2})')
    
    query_offset = int(qstart)  # using integer positions from input
    target_offset = int(tstart)
    
    for m in pattern.finditer(cs_str):
        if m.group('match'):
            # match length
            length = int(m.group(2))
            query_offset += length
            target_offset += length
        elif m.group('ins'):
            seq = m.group(4)
            event = {
                'type': 'ins',
                'seq': seq,
                'pav_start': query_offset,
                'pav_end': query_offset + len(seq) - 1
            }
            # In an insertion, only the query moves forward
            query_offset += len(seq)
            yield event
        elif m.group('del'):
            # deletion: remove the '-' and get the sequence
            seq = m.group('del')[1:]
            event = {
                'type': 'del',
                'seq': seq,
                'pav_start': target_offset,
                'pav_end': target_offset + len(seq) - 1
            }
            # In a deletion, only the target moves forward
            target_offset += len(seq)
            yield event
        elif m.group('sub'):
            # substitution: update both offsets by 1 and ignore for PAV output
            query_offset += 1
            target_offset += 1

def process_paf(paf_file, out_file, cutoff):
    """
    Process the input PAF file and write out a FASTA file of PAV events.
    
    For each line, extract:
      Query sequence, Query start, Query end,
      Target sequence, Target start, Target end,
      and the cs:Z: field.
      
    Then parse the cs:Z: string to extract insertion (using query coords)
    and deletion (using target coords) events.
    
    Only events with length >= cutoff are output.
    
    The FASTA header format is:
      >[sequenceID]:[PAV_start]..[PAV_end]
    and the sequence is the inserted/deleted sequence.
    For an insertion, sequenceID is Query sequence; for a deletion, it is Target sequence.
    """
    with open(paf_file, 'r') as fin, open(out_file, 'w') as fout:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            # Based on the example, the columns are:
            # 0: Query sequence ID
            # 2: Query start
            # 3: Query end
            # 5: Target sequence ID
            # 7: Target start
            # 8: Target end
            # cs:Z: field should be present in one of the optional fields (we look for a field starting with "cs:Z:")
            query_id = fields[0]
            qstart = fields[2]
            qend = fields[3]
            target_id = fields[5]
            tstart = fields[7]
            tend = fields[8]
            cs_field = None
            for field in fields[12:]:
                if field.startswith("cs:Z:"):
                    cs_field = field[5:]  # remove the "cs:Z:" prefix
                    break
            if cs_field is None:
                continue  # no cs field, skip line

            # Process the cs string to get events
            for event in parse_cs(cs_field, qstart, tstart):
                if len(event['seq']) < cutoff:
                    continue
                if event['type'] == 'ins':
                    header = f">{query_id}:{event['pav_start']}..{event['pav_end']}"
                elif event['type'] == 'del':
                    header = f">{target_id}:{event['pav_start']}..{event['pav_end']}"
                fout.write(header + "\n")
                fout.write(event['seq'] + "\n")

def main():
    parser = argparse.ArgumentParser(description="Convert PAF alignment with cs:Z: differences to a PAV FASTA file.")
    parser.add_argument("-paf", required=True, help="Input PAF file")
    parser.add_argument("-out", required=True, help="Output FASTA file")
    parser.add_argument("-c", type=int, default=25, help="Size cutoff for PAVs (default: 25 bp)")
    args = parser.parse_args()
    
    process_paf(args.paf, args.out, args.c)

if __name__ == "__main__":
    main()
