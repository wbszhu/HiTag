import sys

ValidPairs=sys.argv[1]
PeakFile=sys.argv[2]
OutDir=sys.argv[3]
LowDistThr=sys.argv[4]
ChrSizeFile=sys.argv[5]
UseP2PBackgrnd=sys.argv[6]

content_to_addvalue = {
    "ValidPairs":ValidPairs,
    "PeakFile":PeakFile,
    "OutDir":OutDir,
    "LowDistThr":LowDistThr,
    "ChrSizeFile":ChrSizeFile,
    "UseP2PBackgrnd":UseP2PBackgrnd,
}


temp = """
#==================================== 
# Sample configuration file for running FitHiChIP
#====================================  

#*****************************
# important parameters
#*****************************

# File containing the valid pairs from HiCPro pipeline 
# Can be either a text file, or a gzipped text file 
ValidPairs={ValidPairs}

# File containing the bin intervals (according to a specified bin size)
# which is an output of HiC-pro pipeline
# If not provided, this is computed from the parameter 1
Interval=

# File storing the contact matrix (output of HiC-pro pipeline)
# should be accompanied with the parameter 2
# if not specified, computed from the parameter 1
Matrix=

# Pre-computed locus pair file
# of the format: 
# chr1 	start1 	end1 	chr2 	start2 	end2 	contactcounts
Bed=


# File containing reference ChIP-seq / HiChIP peaks (in .bed format)
# mandatory parameter
PeakFile={PeakFile}

# Output base directory under which all results will be stored
OutDir={OutDir}

#Interaction type - 1: peak to peak 2: peak to non peak 3: peak to all (default) 4: all to all 5: everything from 1 to 4.
IntType=3

# Size of the bins [default = 5000], in bases, for detecting the interactions.
BINSIZE=5000

# Lower distance threshold of interaction between two segments
# (default = 20000 or 20 Kb)
LowDistThr=20000

# Upper distance threshold of interaction between two segments
# (default = 2000000 or 2 Mb)
UppDistThr=2000000

# Applicable only for peak to all output interactions - values: 0 / 1
# if 1, uses only peak to peak loops for background modeling - corresponds to FitHiChIP(S)
# if 0, uses both peak to peak and peak to nonpeak loops for background modeling - corresponds to FitHiChIP(L)
UseP2PBackgrnd={UseP2PBackgrnd}

# parameter signifying the type of bias vector - values: 1 / 2
# 1: coverage bias regression	2: ICE bias regression
BiasType=1

# following parameter, if 1, means that merge filtering (corresponding to either FitHiChIP(L+M) or FitHiChIP(S+M))
# depending on the background model, would be employed. Otherwise (if 0), no merge filtering is employed. Default: 1
MergeInt=1

# FDR (q-value) threshold for loop significance
QVALUE=0.01

# File containing chromomosome size values corresponding to the reference genome.
ChrSizeFile={ChrSizeFile}

# prefix string of all the output files (Default = 'FitHiChIP').
PREFIX=FitHiChIP

# Binary variable 1/0: if 1, overwrites any existing output file. otherwise (0), does not overwrite any output file.
OverWrite=1
"""


new_file_content = temp.format(**content_to_addvalue)

print(new_file_content)
