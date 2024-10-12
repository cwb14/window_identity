#!/usr/bin/env python

import sys
import re
import argparse

def parse_paf_line(line):
    fields = line.strip().split('\t')
    if len(fields) < 12:
        return None  # Not a valid PAF line
    query_name = fields[0]
    query_length = int(fields[1])
    query_start = int(fields[2])
    query_end = int(fields[3])
    strand = fields[4]
    target_name = fields[5]
    target_length = int(fields[6])
    target_start = int(fields[7])
    target_end = int(fields[8])
    num_match = int(fields[9])
    block_length = int(fields[10])
    mapping_quality = int(fields[11])
    tags = {}
    for field in fields[12:]:
        parts = field.split(':', 2)
        if len(parts) == 3:
            tag, type_, value = parts
            tags[f'{tag}:{type_}'] = value
    return {
        'query_name': query_name,
        'query_length': query_length,
        'query_start': query_start,
        'query_end': query_end,
        'strand': strand,
        'target_name': target_name,
        'target_length': target_length,
        'target_start': target_start,
        'target_end': target_end,
        'num_match': num_match,
        'block_length': block_length,
        'mapping_quality': mapping_quality,
        'tags': tags
    }

def parse_cs_tag(cs_string, query_start, query_end, target_start, target_end, strand):
    # Initialize positions
    if strand == '+':
        query_pos = query_start
        target_pos = target_start
        query_step = 1
        target_step = 1
    else:
        query_pos = query_end - 1  # Since we will be moving backwards
        target_pos = target_start
        query_step = -1
        target_step = 1
    # Regular expression to parse the cs:Z: tag
    pattern = re.compile(r'([:=*+-])(\d+|[a-z]+)', re.I)
    tokens = pattern.findall(cs_string)
    mismatches = []
    for op, val in tokens:
        if op == ':' or op == '=':
            # Match, advance positions
            length = int(val) if op == ':' else len(val)
            query_pos += length * query_step
            target_pos += length * target_step
        elif op == '*':
            # Substitution
            ref_base = val[0]
            query_base = val[1]
            mismatches.append({
                'ref_pos': target_pos,
                'ref_base': ref_base,
                'query_pos': query_pos,
                'query_base': query_base
            })
            query_pos += query_step
            target_pos += target_step
        elif op == '+':
            # Insertion in query relative to reference
            insertion = val
            query_pos += len(insertion) * query_step
            # target_pos remains the same
        elif op == '-':
            # Deletion from reference
            deletion = val
            target_pos += len(deletion) * target_step
            # query_pos remains the same
        else:
            # Unknown operation
            pass
    return mismatches

def main():
    parser = argparse.ArgumentParser(description='Process PAF file and output mismatches.')
    parser.add_argument('-paf', '--paf_file', required=True, help='Input PAF file')
    args = parser.parse_args()
    paf_file = args.paf_file
    with open(paf_file, 'r') as f:
        for line in f:
            result = parse_paf_line(line)
            if not result:
                continue
            if 'cs:Z' not in result['tags']:
                continue
            cs_string = result['tags']['cs:Z']
            mismatches = parse_cs_tag(cs_string, result['query_start'], result['query_end'],
                                      result['target_start'], result['target_end'], result['strand'])
            # For each mismatch, output the required info
            for m in mismatches:
                ref_contig = result['target_name']
                ref_start = m['ref_pos']
                ref_end = m['ref_pos'] + 1
                ref_allele = m['ref_base']
                query_allele = m['query_base']
                query_contig = result['query_name']
                query_start = m['query_pos']
                query_end = m['query_pos'] + 1
                print('\t'.join(map(str, [ref_contig, ref_start, ref_end, ref_allele, query_allele,
                                          query_contig, query_start, query_end])))
if __name__ == '__main__':
    main()
