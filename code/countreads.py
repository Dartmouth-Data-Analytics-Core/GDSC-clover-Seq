#!/usr/bin/env python3
"""
countreads.py — Count tRNA and other RNA features from sample BAM files.

For each sample BAM, reads are tallied per tRNA gene (mature tRNA and pre-tRNA
loci), Ensembl biotypes, and arbitrary BED features. Outputs gene-level count
tables, fragment-type breakdowns (whole / 5' / 3'), CCA-end distributions, and
amino/anticodon-level isotype counts. Parallelises across samples via
multiprocessing.Pool when --cores > 1.
"""

import pysam
import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *
import itertools
import threading
import time
from multiprocessing import Process, Queue, Pool



def getdupes(namelist):
    allset = set()
    for currname in namelist:
        if currname in allset:
            yield currname
        else:
            allset.add(currname)

def enddict():
    return defaultdict(int)
class featurecount:
    """Per-sample tRNA and feature read count accumulator.

    Tracks counts at multiple resolution levels: total tRNA reads, per-gene
    counts, fragment types (whole/5'/3'), CCA-end types, amino-acid and
    anticodon isotype counts, and average read lengths per feature.
    """
    def __init__(self, samplename, bamfile, trnas = list(), trnaloci = list(), emblgenes = list(), otherfeats = list()):
        self.samplename = samplename
        self.bamfile = bamfile
        self.trnas = trnas
        self.trnaloci = trnaloci
        self.emblgenes = emblgenes
        self.otherfeats = otherfeats

        self.counts = defaultdict(int)
        self.trnacounts = defaultdict(int)
        self.antitrnacount = defaultdict(int)
        self.trnawholecounts = defaultdict(int)
        self.trnafivecounts = defaultdict(int)
        self.trnathreecounts = defaultdict(int)
        self.trnalocuscounts = defaultdict(int)
        self.trnalocustrailercounts = defaultdict(int)
        self.partialtrnalocuscounts = defaultdict(int)
        self.fulltrnalocuscounts  = defaultdict(int)
        self.trnauniquecounts = defaultdict(int)
        self.aminocounts  = defaultdict(int)
        self.anticodoncounts =  defaultdict(int)
        self.trnaendtypecounts = defaultdict(enddict)
        self.lengthsum = defaultdict(int)
        self.lengthtotal = defaultdict(int)

        self.gcpercent = defaultdict(int)
        self.gctotal = defaultdict(int)

        self.genetypes = dict()

    def setgenetype(self, genename, genetype):
        self.genetypes[genename] = genetype
    def addcount(self, genename):
       self.counts[genename] += 1
    def addantitrnacount(self, genename):
       self.antitrnacount[genename] += 1

    def addlocuscount(self, genename):
       self.trnalocuscounts[genename] += 1
    def addpartiallocuscount(self, genename):
       self.partialtrnalocuscounts[genename] += 1
    def addfulllocuscount(self, genename):
       self.fulltrnalocuscounts[genename] += 1
    def addlocustrailercount(self, genename):
       self.trnalocustrailercounts[genename] += 1
    def addtrnacount(self, genename):
       self.trnacounts[genename] += 1
    def adduniquecount(self, genename):
       self.trnauniquecounts[genename] += 1
    def addaminocount(self, amino):
       self.aminocounts[amino] += 1
    def addanticodoncount(self, anticodon):
       self.anticodoncounts[anticodon] += 1
    def addfragcount(self, featname, fragtype):
        if fragtype == "Whole":
            self.trnawholecounts[featname] += 1
        elif fragtype == "Fiveprime":
            self.trnafivecounts[featname] += 1
        elif fragtype == "Threeprime":
            self.trnathreecounts[featname] += 1
    def addendcount(self, featname, endtype):
        if endtype is not None:
            self.trnaendtypecounts[featname][endtype] += 1

    def addreadlength(self, genename, length):
       self.lengthsum[genename] += length
       self.lengthtotal[genename] += 1

    def addgc(self, genename, gc, length):
       self.gcpercent[genename] += gc
       self.gctotal[genename] += length


    def getgenecount(self, genename):
       return self.counts[genename]
    def getantitrnacount(self, genename):
       return self.antitrnacount[genename]
    def getlocuscount(self, genename):
       return self.trnalocuscounts[genename]
    def getpartiallocuscount(self, genename):
       return self.partialtrnalocuscounts[genename]
    def getfulllocuscount(self, genename):
       return self.fulltrnalocuscounts[genename]
    def getlocustrailercount(self, genename):
       return self.trnalocustrailercounts[genename]
    def gettrnacount(self, genename):
       return self.trnacounts[genename]
    def getuniquecount(self, genename):
       return self.trnauniquecounts[genename]
    def getaminocount(self, amino):
       return self.aminocounts[amino]
    def getanticodoncount(self, anticodon):
       return self.anticodoncounts[anticodon]

    def getfivecount(self, genename):
       return self.trnafivecounts[genename]
    def getthreecount(self, genename):
       return self.trnathreecounts[genename]
    def getwholecount(self, genename):
       return self.trnawholecounts[genename]
    def getendtypecount(self, genename):
       return self.trnaendtypecounts[genename]
    def getreadavglength(self, genename):
        if self.lengthtotal[genename] == 0:
            return 0
        else:
            return self.lengthsum[genename] / self.lengthtotal[genename]
    def getreadavggc(self, genename):
        if self.gctotal[genename] == 0:
            return 0
        else:
            return self.gcpercent[genename] / self.gctotal[genename]

def getbamcounts(bamfile, samplename, trnainfo, trnaloci, trnalist, featurelist=dict(), otherseqdict=dict(), embllist=list(), nomultimap=False, allowindels=True, maxmismatches=None):
    """Count reads per tRNA and feature for one sample BAM.

    Iterates mature tRNA features (trnalist), pre-tRNA loci (trnaloci), Ensembl
    genes (embllist), and arbitrary BED features (featurelist), accumulating
    read counts and fragment types in a featurecount object. Returns the
    populated featurecount for the sample.
    """
    samplecounts = featurecount(samplename, bamfile, trnas=trnalist, trnaloci=trnaloci, emblgenes=embllist, otherfeats=featurelist)
    fullpretrnathreshold = 2
    minpretrnaextend = 5
    minmapq = 2 if nomultimap else 0
    minreads = 5
    genetypes = dict()
    currbam = bamfile

    try:
        if not os.path.isfile(currbam+".bai") or os.path.getmtime(currbam+".bai") < os.path.getmtime(currbam):
            pysam.index(""+currbam)
        bamfile = pysam.Samfile(""+currbam, "rb")
    except IOError as xxx_todo_changeme1:
        (strerror) = xxx_todo_changeme1
        print(strerror, file=sys.stderr)
        sys.exit(1)
    except pysam.utils.SamtoolsError:
        print("Can not index "+currbam, file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(1)

    for currfile in featurelist.keys():
        for currfeat in featurelist[currfile]:
            try:
                for currread in getbamrange(bamfile, currfeat, singleonly=nomultimap, maxmismatches=maxmismatches, allowindels=allowindels):
                    if currfeat.coverage(currread) > 10:
                        samplecounts.addcount(currfeat.name)
                        samplecounts.addreadlength(currfeat.name, currread.length())
                        samplecounts.setgenetype(currfeat.name, os.path.basename(currfile))
            except ValueError:
                pass

    for currtype in otherseqdict.keys():
        for currfeat in otherseqdict[currtype]:
            for currread in getbamrange(bamfile, currfeat, singleonly=nomultimap, maxmismatches=maxmismatches, allowindels=allowindels):
                samplecounts.addcount(currfeat.name)
                samplecounts.addreadlength(currfeat.name, currread.length())
                samplecounts.setgenetype(currfeat.name, currtype)

    for genename, featset in itertools.groupby(embllist, lambda x: x.data["genename"]):
        try:
            for currfeat in list(featset):
                for currread in getbamrangeshort(bamfile, currfeat, singleonly=nomultimap, maxmismatches=maxmismatches, allowindels=allowindels, skiptags=True):
                    if currfeat.coverage(currread) > 10:
                        samplecounts.addcount(genename)
                        samplecounts.addreadlength(currfeat.name, currread.length())
                        samplecounts.setgenetype(genename, currfeat.data["biotype"])
        except ValueError:
            pass

    for currfeat in trnaloci:
        for currread in getbamrangeshort(bamfile, currfeat.addmargin(30), singleonly=nomultimap, maxmismatches=maxmismatches, allowindels=allowindels, skiptags=True):
            if currfeat.coverage(currread) > 10 and (currread.start + minpretrnaextend <= currfeat.start or currread.end - minpretrnaextend >= currfeat.end):
                samplecounts.addlocuscount(currfeat.name)
                samplecounts.addreadlength(currfeat.name, currread.length())
                if currread.start + fullpretrnathreshold < currfeat.start and currread.end - fullpretrnathreshold + 3 > currfeat.end:
                    samplecounts.addfulllocuscount(currfeat.name)
                else:
                    samplecounts.addpartiallocuscount(currfeat.name)
            elif currfeat.getdownstream(30).coverage(currread) > 10:
                samplecounts.addlocuscount(currfeat.name)
                samplecounts.addlocustrailercount(currfeat.name)

    for currfeat in trnalist:
        featreads = 0
        for currread in getbam(bamfile, currfeat, singleonly=nomultimap, allowindels=allowindels):
            if maxmismatches is not None and currread.getmismatches() > maxmismatches:
                continue
            samplecounts.addreadlength(currfeat.name, currread.length())
            featreads += 1
            if not currfeat.strand == currread.strand:
                samplecounts.addantitrnacount(currfeat.name)
                continue
            if not currfeat.coverage(currread) > 10:
                continue
            curramino = trnainfo.getamino(currfeat.name)
            curranticodon = trnainfo.getanticodon(currfeat.name)
            samplecounts.addtrnacount(currfeat.name)
            fragtype = getfragtype(currfeat, currread)
            samplecounts.addfragcount(currfeat.name, fragtype)
            endtype = getendtype(currfeat, currread)
            samplecounts.addendcount(currfeat.name, endtype)
            if currread.isuniquetrnamapping():
                samplecounts.adduniquecount(currfeat.name)
            if currread.isuniqueaminomapping():
                samplecounts.addaminocount(curramino)
                if currread.isuniqueacmapping():
                    samplecounts.addanticodoncount(curranticodon)

    return samplecounts

def printcountfile(countfile, samples, samplecounts, trnalist, trnaloci, featurelist, embllist, otherseqdict=dict(), minreads=5, includebase=False):
    """Write gene-level count table to countfile, filtering features below minreads.

    When includebase is True, writes simple total counts per gene. When False,
    writes fragment-type breakdown (whole/5'/3'/other/antisense) per gene.
    """
    print("\t".join(samples), file=countfile)
    trnanames = set()
    for currfeat in trnalist:
        if max(itertools.chain((samplecounts[currsample].gettrnacount(currfeat.name) for currsample in samples), [0])) < minreads:
            continue
        if includebase:
            print(currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].gettrnacount(currfeat.name)) for currsample in samples), file=countfile)
            print(currfeat.name+"_antisense\t"+"\t".join(str(samplecounts[currsample].getantitrnacount(currfeat.name)) for currsample in samples), file=countfile)
        else:
            print(currfeat.name+"_wholecounts\t"+"\t".join(str(samplecounts[currsample].getwholecount(currfeat.name)) for currsample in samples), file=countfile)
            print(currfeat.name+"_fiveprime\t"+"\t".join(str(samplecounts[currsample].getfivecount(currfeat.name)) for currsample in samples), file=countfile)
            print(currfeat.name+"_threeprime\t"+"\t".join(str(samplecounts[currsample].getthreecount(currfeat.name)) for currsample in samples), file=countfile)
            print(currfeat.name+"_other\t"+"\t".join(str(samplecounts[currsample].gettrnacount(currfeat.name) - (samplecounts[currsample].getwholecount(currfeat.name) + samplecounts[currsample].getfivecount(currfeat.name) + samplecounts[currsample].getthreecount(currfeat.name))) for currsample in samples), file=countfile)
            print(currfeat.name+"_antisense\t"+"\t".join(str(samplecounts[currsample].getantitrnacount(currfeat.name)) for currsample in samples), file=countfile)

    for currfeat in trnaloci:
        if max(itertools.chain((samplecounts[currsample].getlocuscount(currfeat.name) for currsample in samples), [0])) < minreads:
            continue
        if includebase:
            print(currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].getlocuscount(currfeat.name)) for currsample in samples), file=countfile)
        else:
            print(currfeat.name+"_wholeprecounts\t"+"\t".join(str(samplecounts[currsample].getfulllocuscount(currfeat.name)) for currsample in samples), file=countfile)
            print(currfeat.name+"_partialprecounts\t"+"\t".join(str(samplecounts[currsample].getpartiallocuscount(currfeat.name)) for currsample in samples), file=countfile)
            print(currfeat.name+"_trailercounts\t"+"\t".join(str(samplecounts[currsample].getlocustrailercount(currfeat.name)) for currsample in samples), file=countfile)

    for currbed in featurelist.keys():
        for currfeat in featurelist[currbed]:
            if currfeat.name in trnanames:
                continue
            trnanames.add(currfeat.name)
            if max(samplecounts[currsample].getgenecount(currfeat.name) for currsample in samples) > minreads:
                print(currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].getgenecount(currfeat.name)) for currsample in samples), file=countfile)
    for currtype in otherseqdict.keys():
        for currfeat in otherseqdict[currtype]:
            trnanames.add(currfeat.name)
            if max(samplecounts[currsample].getgenecount(currfeat.name) for currsample in samples) > minreads:
                print(currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].getgenecount(currfeat.name)) for currsample in samples), file=countfile)
    for currfeat in embllist:
        genename = currfeat.data['genename']
        if genename in trnanames:
            continue
        trnanames.add(genename)
        if genename is None:
            print(currfeat.name, file=sys.stderr)
            sys.exit(1)
        if max(samplecounts[currsample].getgenecount(genename) for currsample in samples) > minreads:
            print(genename+"\t"+"\t".join(str(samplecounts[currsample].getgenecount(genename)) for currsample in samples), file=countfile)


def averagesamples(allcounts, genename, samples):
    return str(sum(allcounts[currsample].lengthsum[genename] for currsample in samples)/(.01+1.*sum(allcounts[currsample].lengthtotal[genename] for currsample in samples)))

def gcsamples(allcounts, genename, samples):
    return str(sum(allcounts[currsample].gcpercent[genename] for currsample in samples)/(.01+1.*sum(allcounts[currsample].gctotal[genename] for currsample in samples)))

def printtypefile(genetypeout, samples, allcounts, trnalist, trnaloci, featurelist, embllist, otherseqdict=dict(), minreads=5):
    """Write gene-to-biotype mapping table with chromosome and average read length.

    Each row: gene_name, biotype, chromosome, mean_read_length. Used downstream
    for DESeq2 grouping and size-factor normalisation.
    """
    trnanames = set()
    genetypes = dict()
    genelengths = dict()
    for currsample in samples:
        genetypes.update(allcounts[currsample].genetypes)
    for currbed in featurelist.keys():
        for currfeat in featurelist[currbed]:
            if currfeat.name in trnanames:
                continue
            trnanames.add(currfeat.name)
            if max(allcounts[currsample].counts[currfeat.name] for currsample in samples) > minreads:
                print(currfeat.name+"\t"+genetypes[currfeat.name]+"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)

    for currfeat in trnaloci:
        print(currfeat.name+"_wholeprecounts"+"\t"+"trna_wholeprecounts"+"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+"_partialprecounts"+"\t"+"trna_partialprecounts"+"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+"_trailercounts"+"\t"+"trna_trailercounts"+"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+"\t"+"tRNA_locus"+"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
    for currfeat in trnalist:
        print(currfeat.name+"_wholecounts"+"\t"+"trna_wholecounts"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+"_fiveprime"+"\t"+"trna_fiveprime"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+"_threeprime"+"\t"+"trna_threeprime"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+"_other"+"\t"+"trna_other"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+"_antisense"+"\t"+"trna_antisense"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+"\t"+"tRNA"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)

    for currfeat in embllist:
        genename = currfeat.data['genename']
        if genename in trnanames:
            continue
        trnanames.add(genename)
        if genename is None:
            continue
        if max(allcounts[currsample].counts[genename] for currsample in samples) > minreads:
            print(genename+"\t"+genetypes[genename]+"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
    for currtype in otherseqdict.keys():
        for currfeat in otherseqdict[currtype]:
            genename = currfeat.name
            if genename in trnanames:
                continue
            trnanames.add(genename)
            if genename is None:
                continue
            if max(allcounts[currsample].counts[genename] for currsample in samples) > minreads:
                print(genename+"\t"+genetypes[genename]+"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)


def printtrnauniquecountcountfile(trnauniquefile, samples, samplecounts, trnalist, trnaloci, minreads=5):
    """Write uniquely-mapped tRNA read counts per sample (one tRNA per row)."""
    trnauniquefile = open(trnauniquefile, "w")
    print("\t".join(currsample for currsample in samples), file=trnauniquefile)
    for currfeat in trnalist:
        if max(samplecounts[currsample].getuniquecount(currfeat.name) for currsample in samples) < minreads:
            continue
        print(currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].getuniquecount(currfeat.name)) for currsample in samples), file=trnauniquefile)
    trnauniquefile.close()

def printtrnacountfile(trnacountfilename, samples, samplecounts, trnalist, trnaloci, minreads=5):
    """Write tRNA locus and mature tRNA read counts per sample (one tRNA per row)."""
    trnacountfile = open(trnacountfilename, "w")
    print("\t".join(currsample for currsample in samples), file=trnacountfile)
    for currfeat in trnaloci:
        if max(samplecounts[currsample].getlocuscount(currfeat.name) for currsample in samples) < minreads:
            continue
        print(currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].getlocuscount(currfeat.name)) for currsample in samples), file=trnacountfile)
    for currfeat in trnalist:
        if max(samplecounts[currsample].gettrnacount(currfeat.name) for currsample in samples) < minreads:
            continue
        print(currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].gettrnacount(currfeat.name)) for currsample in samples), file=trnacountfile)
    trnacountfile.close()

trnaends = list(["CCA","CC","C",""])
def printtrnaendfile(trnaendfilename, samples, samplecounts, trnalist, trnaloci, minreads=5):
    """Write CCA-end type counts per tRNA and sample.

    End types are 'CCA' (full), 'CC', 'C', and 'Trimmed' (no CCA). Only tRNAs
    with at least minreads total counts across samples are included.
    """
    if trnaendfilename is None:
        return
    trnaendfile = open(trnaendfilename, "w")
    print("end\t"+"\t".join(currsample for currsample in samples), file=trnaendfile)
    for currfeat in trnalist:
        if max(samplecounts[currsample].gettrnacount(currfeat.name) for currsample in samples) < minreads:
            continue
        for currend in trnaends:
            endstring = currend if currend != "" else "Trimmed"
            print(currfeat.name+"\t"+endstring+"\t"+"\t".join(str(samplecounts[currsample].getendtypecount(currfeat.name)[currend]) for currsample in samples), file=trnaendfile)
    trnaendfile.close()

def getbamcountsthr(results, currsample, *args, **kwargs):
    results[currsample] = getbamcounts(*args, **kwargs)
def getbamcountsqueue(countqueue, currsample, *args, **kwargs):
    countqueue.put([currsample, getbamcounts(*args, **kwargs)])

def countreadspool(args):
    return getbamcounts(*args[0], **args[1])

def compressargs(*args, **kwargs):
    return tuple([args, kwargs])

def testmain(**argdict):
    """Main entry point: load reference files, count reads per sample, write output tables."""
    trnauniquefilename = None
    argdict = defaultdict(lambda: None, argdict)
    includebase = argdict["nofrag"]
    fullpretrnasonly = argdict["onlyfullpretrnas"]
    trnatable = argdict["trnatable"]
    removepseudo = argdict["removepseudo"]
    ensemblgtf = argdict["ensemblgtf"]
    nomultimap = argdict["nomultimap"]
    if argdict["maxmismatches"] is not None:
        maxmismatches = int(argdict["maxmismatches"])
    else:
        maxmismatches = None
    cores = argdict["cores"] if argdict["cores"] is not None else 1
    trnaendfilename = argdict["trnaends"]
    threadmode = cores > 1

    if argdict["bamdir"] is not None:
        bamdir = argdict["bamdir"]
    else:
        bamdir = "./"
    sampledata = samplefile(argdict["samplefile"], bamdir=bamdir)

    bedfiles = list()
    if argdict["trnauniquecounts"] is not None:
        trnauniquefilename = argdict["trnauniquecounts"]
    if argdict["bedfile"] is not None:
        bedfiles = argdict["bedfile"]
    trnalocifiles = list()
    if argdict["trnaloci"] is not None:
        trnalocifiles = argdict["trnaloci"]
    maturetrnas = list()
    if argdict["maturetrnas"] is not None:
        maturetrnas = argdict["maturetrnas"]

    genetypefile = argdict["genetypefile"]
    trnacountfilename = argdict["trnacounts"]

    trnainfo = transcriptfile(trnatable)

    samples = sampledata.getsamples()
    otherseqdict = dict()
    try:
        featurelist = dict()
        trnaloci = list()
        for currfile in bedfiles:
            bedfeatures = list(readfeatures(currfile, removepseudo=removepseudo))
            featurelist[currfile] = bedfeatures
        trnalist = list()
        for currfile in trnalocifiles:
            trnaloci.extend(list(readbed(currfile)))
        for currfile in maturetrnas:
            trnalist.extend(list(readbed(currfile)))
        if ensemblgtf is not None:
            embllist = list(readgtf(ensemblgtf, filterpsuedo=removepseudo))
        else:
            embllist = list()
    except IOError as e:
        print(e, file=sys.stderr)
        sys.exit()

    allfeats = trnaloci + trnalist
    if len(set(curr.name for curr in allfeats)) < len(list(curr.name for curr in allfeats)):
        print("Duplicate names in feature list", file=sys.stderr)

    allcounts = dict()
    starttime = time.time()
    print(maxmismatches, file=sys.stderr)

    if threadmode:
        countpool = Pool(processes=cores)
        arglist = list()
        for currsample in samples:
            currbam = sampledata.getbam(currsample)
            arglist.append(compressargs(currbam, currsample, trnainfo, trnaloci, trnalist, otherseqdict=otherseqdict, embllist=embllist, featurelist=featurelist, maxmismatches=maxmismatches))
        results = countpool.map(countreadspool, arglist)
        for i, curr in enumerate(samples):
            allcounts[curr] = results[i]
    else:
        for currsample in samples:
            currbam = sampledata.getbam(currsample)
            allcounts[currsample] = getbamcounts(currbam, currsample, trnainfo, trnaloci, trnalist, otherseqdict=otherseqdict, embllist=embllist, featurelist=featurelist, maxmismatches=maxmismatches)

    if argdict["countfile"] is None or argdict["countfile"] == "stdout":
        countfile = sys.stdout
    else:
        countfile = open(argdict["countfile"], "w")
    printcountfile(countfile, samples, allcounts, trnalist, trnaloci, featurelist, embllist, otherseqdict=otherseqdict, includebase=includebase)

    if genetypefile is not None:
        genetypeout = open(genetypefile, "w")
        printtypefile(genetypeout, samples, allcounts, trnalist, trnaloci, featurelist, embllist, otherseqdict=otherseqdict)

    if trnacountfilename is not None:
        printtrnacountfile(trnacountfilename, samples, allcounts, trnalist, trnaloci)
        printtrnaendfile(trnaendfilename, samples, allcounts, trnalist, trnaloci)

    if trnauniquefilename is not None:
        printtrnauniquecountcountfile(trnauniquefilename, samples, allcounts, trnalist, trnaloci)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Count tRNA and other RNA features from BAM files.')
    parser.add_argument('--samplefile',
                       help='Sample file')
    parser.add_argument('--bedfile', nargs='+', default=None,
                       help='BED file(s) with non-tRNA features')
    parser.add_argument('--gtffile', nargs='+', default=None,
                       help='GTF file(s) with non-tRNA features')
    parser.add_argument('--ensemblgtf',
                       help='Ensembl GTF file')
    parser.add_argument('--trnaloci', nargs='+', default=None,
                       help='BED file(s) with tRNA loci')
    parser.add_argument('--maturetrnas', nargs='+', default=None,
                       help='BED file(s) with mature tRNA features')
    parser.add_argument('--onlyfullpretrnas', action="store_true", default=False,
                       help='Only include full pre-tRNAs')
    parser.add_argument('--trnatable',
                       help='Table of tRNA features')
    parser.add_argument('--removepseudo', action="store_true", default=False,
                       help='Remove pseudogenes from Ensembl GTFs')
    parser.add_argument('--genetypefile',
                       help='Output file with gene types')
    parser.add_argument('--trnacounts',
                       help='Output file with tRNA gene counts')
    parser.add_argument('--nofrag', action="store_true", default=False,
                       help='Disable fragment determination')
    parser.add_argument('--nomultimap', action="store_true", default=False,
                       help='Do not count multiply mapped reads')
    parser.add_argument('--maxmismatches', default=None,
                       help='Maximum number of allowable mismatches')
    parser.add_argument('--trnaends', default=None,
                       help='Output file for tRNA end counts')
    parser.add_argument('--trnauniquecounts', default=None,
                       help='Output file for unique tRNA counts')
    parser.add_argument('--countfile', default=None,
                       help='Output file for gene-level read counts (default: stdout)')
    parser.add_argument('--cores', type=int, default=1,
                       help='Number of parallel worker processes (default: 1)')

    args = parser.parse_args()
    argvars = vars(args)
    testmain(**argvars)
