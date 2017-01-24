#!/gpfs/gpfs0/tools/biobuilds/BioBuilds-2015.11/bin/python
import sys
import argparse
import re

# default shebang bath is 
# #!/usr/bin/python
def quickPrint(d):
    s = ""
    for k,v in d.items():
        s += str(k) + ","
        s += str(v) + ","
    return s
def quickPrintGeneName(d):
    s = ""
    for k,v in d.items():
        if k == 'gene_name':
            s += str(v)
    return s
def parseDataAsList(l):
    r = []
    with open(l) as f:
        for i in f:
            i.rstrip('\r\n')
            iLine = i.split(',')
            r.append(iLine)
        return r
def parseFormattedDataAsList(l):
    r = []
    with open(l) as f:
        for i in f:
            i = i.rstrip('\n\r')
            r.append(i)
        return r
def genesFromhg19(d, genes):
    emptyD = dict()
    tDict = dict()
    for gItm in genes:
        genesItem = gItm
        for g in genesItem:
            gRegex = re.compile(r"\b" + g + r"\b")
            m = ""
            try:
                m = d['gene_name']
            except KeyError:
                return emptyD
            regex = re.search(gRegex, m)
            if regex:
                tDict[m] = 1
                return (d, tDict)
    return (emptyD, tDict)
def parseGeneAndPathwayData(f, genesList):
    # format is
    # (chrNm, (exonStart, (exon, exonGeneData)))
    tupleList = []
    tupleDict = dict()
    with open(f) as l:
        for i in l:
            i = i.rstrip('\n')
            t1 = i.split('\t') # NEEDTODO: add catch for nonstandard lines
            if t1[2] != "exon":
                continue # only get exon GTF entries at this point 
            tSplit = t1[8]
            tSplit = tSplit.replace('"', '')
            t2 = tSplit.split(';') # NEEDTODO: check for nonstandard lines and debug
            exon = dict()
            exonStartInt = int(t1[3])
            exonStopInt = int(t1[4])
            exon["Start"] = exonStartInt
            exon["Stop"] = exonStopInt
            d = dict()
            for itm in t2:
                itmSplit = itm.split(' ')
                if len(itmSplit) == 2:
                    d[itmSplit[0]] = itmSplit[1]
                elif len(itmSplit) == 3:
                    itmSplit.pop(0)
                    d[itmSplit[0]] = itmSplit[1]
                else:
                    continue
            returnTuple = genesFromhg19(d, genesList)
            if returnTuple[0]:
                tupleList.append((t1[0], (t1[3], (exon, d))))
                t = returnTuple[1]
                tupleDict.update(t)
        return (tupleList, tupleDict)

def parseGeneData(f):
    # format is
    # (chrNm, (exonStart, (exon, exonGeneData)))
    tupleList = []
    with open(f) as l:
        for i in l:
            i = i.rstrip('\n')
            t1 = i.split('\t') # NEEDTODO: add catch for nonstandard lines
            if t1[2] != "exon":
                continue # only get exon GTF entries at this point 
            tSplit = t1[8]
            tSplit = tSplit.replace('"', '')
            t2 = tSplit.split(';') # NEEDTODO: check for nonstandard lines and debug
            exon = dict()
            exonStartInt = int(t1[3])
            exonStopInt = int(t1[4])
            exon["Start"] = exonStartInt
            exon["Stop"] = exonStopInt
            d = dict()
            for itm in t2:
                itmSplit = itm.split(' ')
                if len(itmSplit) == 2:
                    d[itmSplit[0]] = itmSplit[1]
                elif len(itmSplit) == 3:
                    itmSplit.pop(0)
                    d[itmSplit[0]] = itmSplit[1]
                else:
                    continue
            tupleList.append((t1[0], (t1[3], (exon, d))))
        return tupleList
def parseMethylatedData(l):
    # format is 
    # [(Chromosome, ({'Start': 0, 'Stop': 1}, (int, int)))]
    mSiteList = []
    for i in l:
        iLine = i
        for j in iLine:
            j = j.strip('\n')
            iList = j.split('\t')
            cName = iList[0]
            mStart = int(iList[1])
            mStop = int(iList[2])
            mStartAndStop = dict()
            mStartAndStop["Start"] = mStart
            mStartAndStop["Stop"] = mStop
            mPercent = float(iList[3])
            numMethylated = int(iList[4])
            numNotMethylated = int(iList[5])
            if mPercent > 50.0:
                covCheck = numMethylated + numNotMethylated
                if covCheck > 9:
                    mSiteList.append((cName, (mStartAndStop, (mPercent, (numMethylated, numNotMethylated)))))
                else:
                    continue
            else:
                continue
    return mSiteList
def formatMethylatedData(l):
    # format is
    # [{ 'geneName': 'GENE', 'geneCDS': [1, 100], 
    #    'Methylation': 
    #         [{ 'chromosome': '1', 'site': 0, 'percent': 100}, 
    #          {'chromosome': '1', 'site': 1, 'percent':99}]}]
    annotationFlag = 0
    methylationFlag = 0
    d = dict()
    m = dict()
    m_list = []
    r = []
    geneName = ""
    firstIteration = 1
    for i in l:
        if i == 'MethylationSite':
            methylationFlag = 1
            continue
        elif i == 'SiteAnnotation':
            annotationFlag = 1
            continue
        else:
            if i[0:6] == 'intron':
                annotationFlag = 1
                methylationFlag = 1
                continue
            if annotationFlag:
                iLine = i.split(',')
                annotationFlag = 0
                if iLine[4] == geneName:
                    continue
                else:
                    if firstIteration:
                        firstIteration = 0
                        geneName = iLine[4]
                        d["geneName"] = geneName
                        d["geneCDS"] = [iLine[1],iLine[3]]
                    else:
                        d["Methylation"] = m_list
                        r.append(d)
                        d = dict()
                        m_list = []
                        geneName = iLine[4]
                        d["geneName"] = geneName
                        d["geneCDS"] = [iLine[1],iLine[3]]
            elif methylationFlag:
                m = dict()
                iLine = i.split(',')
                methylationFlag = 0
                m["chromosome"] = iLine[1]
                m["site"] = iLine[3]
                m["percent"] = iLine[7]
                m_list.append(m)
            else:
                continue
    # required for last iteration
    if m_list:
        d["Methylation"] = m_list
    r.append(d)
    return r
def main():
    gtfGeneList = parseGeneData('genes.gtf')
    parser = argparse.ArgumentParser(description='Parse Methylation Data')
    parser.add_argument('-i', '--input', dest='input', help='Input bisulfite data file')
    # action taken
    parser.add_argument('-c', '--cmd', dest='cmd', choices=['sortcov', 'sortmethyl', 'sortmethylP'], help='command can be to sort a coverage file from bismark, parse a sorted methylation file, or parse a methylation file using a pathway file')
    parser.add_argument('-f', '--file', dest='fileList', help='if using sortmethylP, add file name of pathway list, otherwise leave blank')
    args = parser.parse_args()
    # ALGORITHM:
    # if ARGS == file:    coverage file (must be tab separated), 
    #                     sorted methylationFile (must be in format MethylationSite\ndata\nSiteAnnotation\ndata),
    #                     sorted MethylationFileWithPathwayFile (must be in format MethylationSite\ndata\nSiteAnnotation\ndata with file named '' in cwd)
    mFileIn = args.input    
    cmd = args.cmd
    fileList = args.fileList
    # import pathway List, if it exists
    pathwayList = []
    if fileList:
        pathwayList = parseDataAsList(fileList)  
    # import methylation data file, if present
    mFileInAsList = []
    if mFileIn:
        mFileInAsList = parseDataAsList(mFileIn)
    else:
        print "Error, a Methylation data file is required but could not be found.  Please check input parameters"
        print mFileIn
        print cmd
        print fileList
        sys.exit()
    if cmd == 'sortcov':
        parsedMethylatedDataAsList = parseMethylatedData(mFileInAsList)
        pLen = len(parsedMethylatedDataAsList)
        methylatedGenes = []
        fDataList = []
        for v in xrange(0, pLen):
            mVal = parsedMethylatedDataAsList[v]
            intronOrExtragenic = 0
            for x in gtfGeneList:
                if x[0] == mVal[0]:
                    checkPos = int(x[1][0])
                    if mVal[1][0]["Start"] > checkPos:
                        if mVal[1][0]["Stop"] < x[1][1][0]["Stop"]:
                            methylatedGenes.append((mVal, (x[1][1][1], x[1][1][0])))
                            intronOrExtragenic = 1
                            break
            if intronOrExtragenic:
                continue
            else:
                methylatedGenes.append((mVal, (0, {'MethylationSite': 'intron or extragenic'})))
        for itm in methylatedGenes:
            mName = itm[0]
            itmParse = itm[1]
            if itmParse[0] == 0:
                s1 = "SiteAnnotation"
                print s1
                fDataList.append(s1)
                s1Name = itmParse[1]['MethylationSite']
                print s1Name
                fDataList.append(s1Name)
                sName = "MethylationSite"
                print sName
                fDataList.append(sName)
                pString =  "Chromosome" + "," + mName[0] + "," + "Start" + "," + str(mName[1][0]['Start']) + "," + "Stop" + "," + str(mName[1][0]['Stop']) + "," + "MethylationPercent" + "," + str(mName[1][1][0])
                print pString
                fDataList.append(pString)
            else:
                sName = "SiteAnnotation"
                fDataList.append(sName)
                s1 = quickPrintGeneName(itmParse[0])
                s = quickPrint(itmParse[1])
                s += s1
                print s
                fDataList.append(s)
                pString =  "Chromosome" + "," + mName[0] + "," + "Start" + "," + str(mName[1][0]['Start']) + "," + "Stop" + "," + str(mName[1][0]['Stop']) + "," + "MethylationPercent" + "," + str(mName[1][1][0])
                pName = "MethylationSite"
                print pName
                print pString
                print "##"
                fDataList.append(pName)
                fDataList.append(pString)
        print ""
        print "########"
        print ""
        formattedMethylatedData = formatMethylatedData(fDataList)
        # write better print routine later
        for k in formattedMethylatedData:
            kItm = k
            for k1,v1 in kItm.items():
                print k1
                print v1
            print "##"
    elif cmd == 'sortmethyl':
        mFileInAsList = parseFormattedDataAsList(mFileIn)
        formattedMethylatedData = formatMethylatedData(mFileInAsList)
        for k in formattedMethylatedData:
            print k 
    elif cmd == 'sortMethylP':
        parsedMethylatedDataAsList = parseMethylatedData(mFileInAsList)
        pLen = len(parsedMethylatedDataAsList)
        tupleOfLists = parseGeneAndPathwayData('genes.gtf', pathwayList)
        itmList = tupleOfLists[0]
        itmDict = tupleOfLists[1]
        mDataLen = len(parsedMethylatedDataAsList)
        methylatedGenes = []
        for v in xrange(0, mDataLen):
            mVal = parsedMethylatedDataAsList[v]
            intronOrExtragenic = 0
            for x in itmList:
                if x[0] == mVal[0]: # same chromosome
                    checkPos = int(x[1][0])
                    if mVal[1][0]["Start"] > checkPos:
                        if mVal[1][0]["Stop"] < x[1][1][0]["Stop"]:
                        # both coord-sorted, so the first hit is a likely match for a region
                        # methylatedGenes.append((mVal, x[1][1][0]))
                        # print x[1][1][1]
                        # print x[1][1][0]
                            print mVal
                            print x[1][1][1]
                            print x[1][1][0]
                            methylatedGenes.append((mVal, (x[1][1][1], x[1][1][0])))
                            intronOrExtragenic = 1
                            break
            if intronOrExtragenic:
                continue
            else:
                methylatedGenes.append((mVal, (0, {'MethylationSite': 'intron or extragenic'})))
        for itm in methylatedGenes:
            mName = itm[0]
            itmParse = itm[1]
            if itmParse[0] == 0:
                continue # no annotation entry
            checkGeneName = ""
            try:
                checkGeneName = itmParse[0]['gene_name']
            except TypeError:
                print "Need to fix this!"
                print itmParse
                sys.exit()
            if checkGeneName not in itmDict:
                print "WARNING! " + checkGeneName + " not in gene List!"
        for itm in methylatedGenes:
            mName = itm[0]
#        print "MethylationSite" 
#        pString =  "Chromosome" + "," + mName[0] + "," + "Start" + "," + str(mName[1][0]['Start']) + "," + "Stop" + "," + str(mName[1][0]['Stop']) + "," + "MethylationPercent" + "," + str(mName[1][1][0])
#        print pString
#        print "SiteAnnotation"
            itmParse = itm[1]
            if itmParse[0] == 0:
            # s = itmParse[1]['MethylationSite']
            # print s
                continue
            else:
                print "SiteAnnotation"
                s_tmp = quickPrintGeneName(itmParse[0])
                s = quickPrint(itmParse[1])
                s += s_tmp
                print s
                pString =  "Chromosome" + "," + mName[0] + "," + "Start" + "," + str(mName[1][0]['Start']) + "," + "Stop" + "," + str(mName[1][0]['Stop']) + "," + "MethylationPercent" + "," + str(mName[1][1][0])
                print "MethylationSite"
                print pString
                print "##"
    else:
        print 'sorry unrecognized command'
if __name__ == "__main__":
    main()
