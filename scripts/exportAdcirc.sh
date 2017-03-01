#!/bin/bash

#...Compression Format
# 0 = No Compression (tar)
# 1 = gzip (.tar.gz)
# 2 = bzip (.tar.bz)
# 3 = xz (.tar.xz)
# 4 = 7zip LZMA (.7z)
compression=1

#...Get the adcirc root directory
parentdir="$(dirname "$(pwd)")"

#...Generate version_autogen.F
version=$(./adcircVersion.sh $parentdir)
if [ $? != 0 ] ; then
    echo "ERROR: Export Failed @ version stage."
    exit 1
fi

#...Check if version is modified
mod=$(git diff-index --quiet HEAD)
if [ "x$?" != "x0" ] ; then
    echo "WARNING: The version of ADCIRC being exported has been modified from its repository state."
fi

#...Get the name of the current branch
branch=$(git branch | grep \* | cut -d\* -f2 | cut -c2-)
if [ $? != 0 ] ; then
    echo "ERROR: Export Failed @ branch stage."
    exit 1
fi

#...Generate the names of the output files
outputDirectory=adcirc_$branch"_"$version
if [ $compression == 0 ] ; then
    outputFile=adcirc_$branch"_"$version.tar
elif [ $compression == 1 ] ; then
    outputFile=adcirc_$branch"_"$version.tar.gz
elif [ $compression == 2 ] ; then
    outputFile=adcirc_$branch"_"$version.tar.bz2
elif [ $compression == 3 ] ; then
    outputFile=adcirc_$branch"_"$version.tar.xz
elif [ $compression == 4 ] ; then
    outputFile=adcirc_$branch"_"$version.7z
else
    echo "ERROR: Invalid compression method."
    exit 1
fi

#...Move to top level directory
cd ../

#...Generate the changelog
./scripts/generateChangelog.sh -a 
if [ $? != 0 ] ; then
    echo "ERROR: Export Failed @ changelog stage."
    exit 1
fi


#...Check if the folder/file already exists
if [ -d ./scripts/$outputDirectory ] ; then
    echo "ERROR: Output directory already exists."
    exit 1
fi

if [ -s ./scripts/$outputFile ] ; then
    echo "ERROR: Output archive already exists."
    exit 1
fi

#...Make a directory to hold the output data
mkdir ./scripts/$outputDirectory

#...Export the repository snapshot to folder
git archive $branch | tar -x -C ./scripts/$outputDirectory
if [ $? != 0 ] ; then
    echo "ERROR: Export Failed @ git archive stage."
    exit 1
fi

#...Remove default version
rm ./scripts/$outputDirectory/version_default.F

#...Grab the version.F file to include
mv version.F ./scripts/$outputDirectory
if [ $? != 0 ] ; then
    echo "ERROR: Export Failed @ tar stage."
    exit 1
fi

#...Grab the changelog file to include
mv changelog ./scripts/$outputDirectory/.
if [ $? != 0 ] ; then
    echo "ERROR: Export Failed @ tar stage."
    exit 1
fi

#...Create the .tar.gz package
cd ./scripts

if [ $compression == 0 ] ; then
    tar -cf $outputFile $outputDirectory 
    if [ $? != 0 ] ; then
        echo "ERROR: Export Failed @ tar stage."
        exit 1
    fi
elif [ $compression == 1 ] ; then
    tar -czf $outputFile $outputDirectory 
    if [ $? != 0 ] ; then
        echo "ERROR: Export Failed @ tar/gzip stage."
        exit 1
    fi
elif [ $compression == 2 ] ; then
    tar -cjf $outputFile $outputDirectory 
    if [ $? != 0 ] ; then
        echo "ERROR: Export Failed @ tar/bzip2 stage."
        exit 1
    fi
elif [ $compression == 3 ] ; then
    tar -cJf $outputFile $outputDirectory 
    if [ $? != 0 ] ; then
        echo "ERROR: Export Failed @ tar/xz stage."
        exit 1
    fi
elif [ $compression == 4 ] ; then
    7za a -m0=lzma2 $outputFile $outputDirectory >/dev/null 
    if [ $? != 0 ] ; then
        echo "ERROR: Export Failed @ 7Zip stage."
        exit 1
    fi
else
    echo "ERROR: Invalid compression selection."
fi    
rm -r $outputDirectory
