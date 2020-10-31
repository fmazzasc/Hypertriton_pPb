#!/usr/bin/env python3

import os
from ROOT import gROOT

gROOT.SetBatch(True)

gROOT.LoadMacro("GenerateTableFromMC.cc")
gROOT.LoadMacro("GenerateTableFromData.cc")
from ROOT import GenerateTableFromMC, GenerateTableFromData


print("++++++++++++++++++++++++++++++++++++++++++")
print("Generate Signal Table")
print("++++++++++++++++++++++++++++++++++++++++++")
GenerateTableFromMC(True, "mtexp")
print("++++++++++++++++++++++++++++++++++++++++++")
print("Generate Data Table")
print("++++++++++++++++++++++++++++++++++++++++++")
GenerateTableFromData(False)
print("++++++++++++++++++++++++++++++++++++++++++")
print("Generate Like-Sign Backgoundd Table")
print("++++++++++++++++++++++++++++++++++++++++++")
GenerateTableFromData(True)