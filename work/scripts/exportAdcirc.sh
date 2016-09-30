#!/bin/bash

#...Get the adcirc root directory
parentdir="$(dirname "$(pwd)")"
parentdir="$(dirname "$parentdir")"

#...Generate version_autogen.F
version=$(./adcircVersion.sh $parentdir)
if [ $? != 0 ] ; then
    echo "ERROR: Export Failed @ version stage."
    exit 1
fi

#...Get the name of the current branch
branch=$(git branch | grep \* | cut -d\* -f2 | cut -c2-)
if [ $? != 0 ] ; then
    echo "ERROR: Export Failed @ branch stage."
    exit 1
fi

#...Move to top level directory
cd ../../

#...Generate the changelog
./work/scripts/generateChangelog.sh -a 
if [ $? != 0 ] ; then
    echo "ERROR: Export Failed @ changelog stage."
    exit 1
fi

#...Generate the names of the output files
outputDirectory=adcirc_$branch"_"$version
outputFile=adcirc_$branch"_"$version.tar.gz

#...Check if the folder/file already exists
if [ -d work/scripts/$outputDirectory ] ; then
    echo "ERROR: Output directory already exists."
    exit 1
fi

if [ -s work/scripts/$outputFile ] ; then
    echo "ERROR: Output archive already exists."
    exit 1
fi

#...Make a directory to hold the output data
mkdir work/scripts/$outputDirectory

#...Export the repository snapshot to folder
git archive $branch | tar -x -C work/scripts/$outputDirectory
if [ $? != 0 ] ; then
    echo "ERROR: Export Failed @ git archive stage."
    exit 1
fi

#...Grab the version.F file to include
mv version.F work/scripts/$outputDirectory
if [ $? != 0 ] ; then
    echo "ERROR: Export Failed @ tar stage."
    exit 1
fi

#...Grab the changelog file to include
mv changelog work/scripts/$outputDirectory/.
if [ $? != 0 ] ; then
    echo "ERROR: Export Failed @ tar stage."
    exit 1
fi

#...Create the .tar.gz package
cd work/scripts
tar -czf $outputFile $outputDirectory 
if [ $? != 0 ] ; then
    echo "ERROR: Export Failed @ gzip stage."
    exit 1
fi
rm -r $outputDirectory
