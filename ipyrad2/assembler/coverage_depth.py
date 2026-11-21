"""Stream single-sample pileup depth into run-length encoded observed-depth intervals."""
import sys


def encode_depth(source, target, minimum):
    pending = None
    for line in source:
        chrom, pos, depth = line.split()[:3]
        start, depth = int(pos) - 1, int(depth)
        if pending and (chrom != pending[0] or start != pending[2] or depth != pending[3]):
            target.write("\t".join(map(str, pending)) + "\n")
            pending = None
        if depth < minimum:
            continue
        if pending:
            pending[2] = start + 1
        else:
            pending = [chrom, start, start + 1, depth]
    if pending:
        target.write("\t".join(map(str, pending)) + "\n")


if __name__ == "__main__":
    encode_depth(sys.stdin, sys.stdout, int(sys.argv[1]))
