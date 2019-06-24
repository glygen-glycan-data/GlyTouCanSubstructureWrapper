#!/bin/env python27

import re
import os
import sys
import math
import time
import urllib
import urllib2


class GlyTouCan():
    endpt = 'http://ts.glytoucan.org/sparql'
    substr_endpt = 'http://test.ts.glytoucan.org/sparql'
    api = 'https://api.glytoucan.org/'


    def __init__(self):
        self.delaytime = .2
        self.delaybatch = 1
        self.maxattempts = 3
        self._lastrequesttime = 0
        self._lastrequestcount = 0

        self.glygen = open(os.path.dirname(os.path.abspath(__file__)) + "/GlyGen_accessions").read().split()


    def _wait(self, delaytime=None):
        if delaytime != None:
            time.sleep(delaytime)
            return

        if (self._lastrequestcount % self.delaybatch) == 0 and self._lastrequestcount > 0:
            time.sleep(self.delaytime)
        self._lastrequesttime = time.time()
        self._lastrequestcount += 1



    getseq_sparql = """
PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
PREFIX glytoucan: <http://www.glytoucan.org/glyco/owl/glytoucan#>

SELECT DISTINCT ?WURCS_label
WHERE{
  VALUES ?accNum {"%s"}
  ?saccharide glytoucan:has_primary_id ?accNum .

  ?saccharide glycan:has_glycosequence ?wcsSeq .
  ?wcsSeq rdfs:label ?wcsLabel .
  BIND(STR(?wcsLabel) AS ?WURCS_label)
  ?wcsSeq glycan:in_carbohydrate_format glycan:carbohydrate_format_wurcs .
}"""

    def getseq(self, accession, format="wurcs"):
        # assert (format in ("wurcs", "glycoct", "iupac_extended", "iupac_condensed"))
        url = "https://ts.glytoucan.org/sparql?query=" + urllib.quote(self.getseq_sparql%accession)
        data = urllib2.urlopen(url).read()
        wurcs = re.compile(r"<literal>WURCS.*</literal>").findall(data)[0].lstrip("<literal>").rstrip("</literal>")
        return wurcs

    def getsubstructure(self, acc):
        wurcs = self.getseq(acc)
        firstbatch, totalcount = self.getsubstructurebywurcs20each(wurcs)

        if totalcount <= 20 and len(firstbatch) == totalcount:
            return firstbatch
        else:
            if totalcount <= 20 and len(firstbatch) != totalcount:
                raise RuntimeError("Total count of super structure doesn;t match with the result accession count")

            res = firstbatch
            requestcount = int(math.ceil(totalcount/20.0))
            for i in range(1, requestcount):
                print "%s/%s glycan to go" % (totalcount-i*20,totalcount)
                thisbatch, thiscount = self.getsubstructurebywurcs20each(wurcs, offset=i*20)
                if thiscount != totalcount:
                    raise RuntimeError("Total count not match")
                res += thisbatch

            if len(res) != totalcount:
                raise RuntimeError("Total count not match")
            return res

    def getsubstructurebywurcs20each(self, wurcs, offset=0):
        url = "https://glytoucan.org/connect?url="
        purl = "http://stanza.glytoucan.org/stanza/substructure_search?order=ASC&orderkey=AccessionNumber&lang=1&"
        para = {}
        para["wurcs"] = wurcs
        para["offset"] = offset
        purl += "wurcs=" + para["wurcs"]
        purl += "&offset=" + str(para["offset"])

        url += urllib.quote(purl)
        # print url
        attempt = 0
        html = None
        while html == None and attempt < self.maxattempts:
            try:
                attempt += 1
                data = urllib2.urlopen(url)
                html = data.read()
            except:
                # traceback.print_exc()
                self._wait(self.delaytime ** attempt)

        if html == None:
            raise IOError("Cannot get result from glytoucan stanza")

        accs = list(set(re.compile(r"G\d{5}\w\w").findall(html)))

        regex_res = re.compile(r"<hr />\s{1,5}\d{1,10}").findall(html)
        if len(regex_res) != 1:
            raise RuntimeError("Cannot found total count for super structure")
        temp = re.compile(r"\d{1,10}").findall(regex_res[0])
        if len(temp) != 1:
            raise RuntimeError("Cannot found total count for super structure")
        totalSuperStructureCount = int(temp[0])
        return accs, totalSuperStructureCount

    def getSubstructureWithinGlyGen(self, acc):
        res = self.getsubstructure(acc)
        return filter(lambda x:x in self.glygen, res)



if __name__ == "__main__":

    gtc = GlyTouCan()
    if len(sys.argv)>1:
        acc = sys.argv[1]
    else:
        print >> sys.stderr, "No GlyTouCan accession was given"
        sys.exit(1)

    if not re.compile(r"G\d{5}\w{2}").findall(acc):
        raise RuntimeError("GlyTouCan accession format error")

    for acc in gtc.getSubstructureWithinGlyGen(acc):
        print acc