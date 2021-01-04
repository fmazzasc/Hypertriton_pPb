#!/usr/bin/env python3

import os
from ROOT import gROOT

gROOT.SetBatch(True)

gROOT.LoadMacro("GenerateTableFromMC.cc")
gROOT.LoadMacro("GenerateTableFromData.cc")
from ROOT import GenerateTableFromMC, GenerateTableFromData


print("++++++++++++++++++++++++++++++++++++++++++")
print("Generate Mt exp")
print("++++++++++++++++++++++++++++++++++++++++++")
GenerateTableFromMC(True, "mtexp")
print("++++++++++++++++++++++++++++++++++++++++++")
print("Generate pt exp")
print("++++++++++++++++++++++++++++++++++++++++++")
GenerateTableFromMC(True, "ptexp")
print("++++++++++++++++++++++++++++++++++++++++++")
print("Generate bol")
print("++++++++++++++++++++++++++++++++++++++++++")
GenerateTableFromMC(True, "bol")
print("++++++++++++++++++++++++++++++++++++++++++")
print("Generate bw")
print("++++++++++++++++++++++++++++++++++++++++++")
GenerateTableFromMC(True, "bw")
print("++++++++++++++++++++++++++++++++++++++++++")