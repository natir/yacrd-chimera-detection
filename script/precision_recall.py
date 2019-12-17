#!/usr/bin/env python3

import os
import re
import csv
import sys
import argparse

from collections import defaultdict

class Mapping:
    def __init__(self, reference, ref_beg, ref_end, que_beg, que_end):
        self.reference = reference
        self.ref_beg = ref_beg
        self.ref_end = ref_end
        self.que_beg = que_beg
        self.que_end = que_end

    def ovl_on_ref(self, other):
        return Mapping.ovl((self.ref_beg, self.ref_end), (other.ref_beg, other.ref_end))

    def ovl_on_que(self, other):
        return Mapping.ovl((self.que_beg, self.que_end), (other.que_beg, other.que_end))

    @staticmethod
    def ovl(query, target):
        if query[0] <= target[0] and target[0] <= query[1]:
            return True
        if query[0] <= target[1] and target[1] <= query[1]:
            return True
        if query[0] <= target[0] and target[1] <= query[1]:
            return True
        if target[0] <= query[0] and target[1] <= query[1]:
            return True
        
        return False
        
def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()

    parser.add_argument("mapping")
    parser.add_argument("yacrd")
    parser.add_argument("-c", "--circular", type=int, help="if a read map after this position it probably not be a chimera just a read map at the end and begin of genome", default=sys.maxsize)
    
    args = parser.parse_args(args)

    real_chimera = set(get_chimera(args.mapping, args.circular))
    yacrd_chimera = set(parse_yacrd_report(args.yacrd))

    precision = len(real_chimera & yacrd_chimera) / len(yacrd_chimera)
    recall = len(real_chimera & yacrd_chimera) / len(real_chimera)

    print(f"{args.yacrd},{precision},{recall},{2 * ((precision * recall) / (precision + recall))}");

def parse_yacrd_report(filename):
    reader = csv.reader(open(filename, "r"), delimiter='\t')

    for row in reader:
        if row[0] == "Chimeric":
            yield row[1]
    
def parse_mapping(filename):
    read2length = defaultdict(int)
    read2mapping = defaultdict(list)
    
    mapping = csv.reader(open(filename, "r"), delimiter='\t')
    for m in mapping:
        que_name = m[0]
        read2length[que_name] = int(m[1])
        
        ref_nam = m[5]
        ref_sta = int(m[7])
        ref_end = int(m[8])

        que_sta = int(m[2])
        que_end = int(m[3])

        read2mapping[que_name].append(Mapping(ref_nam, ref_sta, ref_end, que_sta, que_end))

    return read2length, read2mapping

def pos_on(pos):
    pos.sort(key=lambda x: x[0])

    end = list()
    while len(pos) >= 1:
        index = 0
        p = pos.pop(0)

        while len(pos) >= 1 and Mapping.ovl(p, pos[0]):
            p = merge_pos(p, pos.pop(0))
            index += 1

        end.append(p)
        
    return end

def pos_on_query(mapping):
    pos = [[m.que_beg, m.que_end] for m in mapping]
    
    return pos_on(pos)

def pos_on_ref(mapping):
    pos = [[m.ref_beg, m.ref_end] for m in mapping]

    return pos_on(pos)

def merge_pos(first, second):
    if second[0] < first[0]:
        first[0] = second[0]
    if second[1] > first[1]:
        first[1] = second[1]
        
    return first

def not_to_fare(len1, len2):
    return max(len1, len2) / max(len1, len2) > 1.1


def get_chimera(filename, circular_limit):
    read2length, read2mapping = parse_mapping(filename)
    
    for read, mapping in read2mapping.items():
        
        if len(read2mapping[read]) < 2 and not_to_fare(mapping[0].que_end - mapping[0].que_end, read2length[read]):
            continue
        
        pos = pos_on_query(mapping)
        if len(pos) < 2:
            continue

        if any([elt > circular_limit for sublist in pos for elt in sublist]):
            continue

        yield read

    
if __name__ == "__main__":
    main(sys.argv[1:])
