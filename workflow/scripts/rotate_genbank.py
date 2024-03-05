#!/usr/bin/env python3
import os
import sys

from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio.Seq import Seq
from Bio import SeqIO


def split_location_on_cut(location, cut_point):
    new_parts = []
    for part in location.parts:
        if cut_point not in part:
            new_parts.append(part)
            continue
        new_parts.append(FeatureLocation(part.start, cut_point, part.strand))
        new_parts.append(FeatureLocation(cut_point, part.end, part.strand))
    return CompoundLocation(new_parts)

assert split_location_on_cut(FeatureLocation(5, 12, -1), 7).parts == [FeatureLocation(5, 7, -1), FeatureLocation(7, 12, -1)]


def shift_location(location, cut_point, padding, rec_len):
    new_parts = []
    for part in location.parts:
        if part.start < cut_point:
            new_parts.append(FeatureLocation(part.start + (rec_len - cut_point) + padding,
                                             part.end + (rec_len - cut_point) + padding,
                                             part.strand))
        else:
            new_parts.append(FeatureLocation(part.start - cut_point,
                                             part.end - cut_point,
                                             part.strand))
    if len(new_parts) == 1:
        return new_parts[0]
    return CompoundLocation(new_parts)


assert shift_location(FeatureLocation(7, 8, -1), 15, 10, 20) == FeatureLocation(22, 23, -1), shift_location(FeatureLocation(7, 8, -1), 15, 10, 20)
assert shift_location(FeatureLocation(16, 17, 1), 15, 10, 20) == FeatureLocation(1, 2, 1), shift_location(FeatureLocation(16, 17, 1), 15, 10, 20)


def rotate(filename, cut_point, padding=0):
    recs = list(SeqIO.parse(filename, 'genbank'))
    if len(recs) > 1:
        raise ValueError("Can only rotate a single record, not multiple")
    rec = recs[0]
    if len(rec) < cut_point:
        raise ValueError("Invalid cut location: outside record")

    new_seq = rec.seq[cut_point:] + Seq("A" * padding) + rec.seq[:cut_point]

    for feature in rec.features:
        if feature.type == "source":
            feature.location = FeatureLocation(feature.location.start, feature.location.end+padding, feature.location.strand)
            continue
        location = feature.location
        if cut_point in location:
            location = split_location_on_cut(location, cut_point)
        location = shift_location(location, cut_point, padding, len(rec))
        assert not (location.start == 0 and location.end == len(rec))
        feature.location = location
    rec.annotations["topology"] = "circular"
    rec.seq = new_seq
    SeqIO.write([rec], filename, "genbank")


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print(f"Usage: {os.path.basename(sys.argv[0])} input_file cut_point [padding_size]", file=sys.stderr)
        exit(1)
    rotate(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 0)
    exit(0)