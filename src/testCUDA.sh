#!/bin/bash

thisDir="$(cd "$(dirname "$0")" && pwd)"
topDir=`dirname ${thisDir}`
binDir=${topDir}/bin
dataDir=${topDir}/Catalogs
outputDir=${topDir}/output

inputFile="Aardvark/Catalog_v1.0/truth_oscillationcorrected_unrotated/Aardvark_v1.0c_truth.190.root"
outputFile="reference_CUDA.root"
${binDir}/cudaref  ${dataDir}/${inputFile}  ${outputDir}/${outputFile} 50000 0 2
