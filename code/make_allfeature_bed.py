#!/usr/bin/env python3

#----- IMPORT LIBRARIES
import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *


def makefeaturebed(maturetRNAs, lociFile, ensGTF):

    #----- Open the output file for writing
    allfeatfile = open("allfeats.bed", "w")

    #----- Process the mature tRNAs, tRNA loci, and smRNA gene features
    for currfeature in readbed(maturetRNAs):
        print(currfeature.bedstring(), file=allfeatfile)
    for currfeature in readbed(lociFile):
        print(currfeature.bedstring(), file=allfeatfile)
    for currfeature in readgtf(ensGTF):
        print(currfeature.bedstring(name = currfeature.data["genename"]), file=allfeatfile)
    allfeatfile.close()

def main():
    #----- Argument parser
    parser = argparse.ArgumentParser(description="Generate a BED file with mature tRNAs, loci, and gene features.")
    parser.add_argument("maturetRNAs", help="Path to the mature tRNA BED file")
    parser.add_argument("lociFile", help="Path to the tRNA loci BED file")
    parser.add_argument("ensGTF", help="Path to the Ensembl GTF file")

    args = parser.parse_args()

    # Call the function with provided arguments
    makefeaturebed(args.maturetRNAs, args.lociFile, args.ensGTF)

if __name__ == "__main__":
    main()