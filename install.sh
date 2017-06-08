#!/bin/bash

message () {
echo; printf '%.0s=' {1..80}; echo
echo "$1"
printf '%.0s=' {1..80}; echo; echo
}

ver="v1.2.2"

# Retrieving CAVA directory and Python version
dir=$(pwd)
pver=$(python -c "import sys; print sys.version[:3]")

if [ $pver != '2.7' ]; then
	echo; echo "CAVA v1.2.2 requires Python 2.7.x"; echo
	exit
fi

# Index genome by Stampy
message "Installing CAVA (Clinical Annotation of VAriants) "$ver

# Creating temporary and final directories for pysam
mkdir -p $dir/tmp/lib/python$pver/site-packages/
mkdir -p $dir/tmp/lib64/python$pver/site-packages/
mkdir $dir/pysamdir

# Unpacking pysam 0.7.7
tar -zxvf pysam-0.7.7.tar.gz

# Building pysam 0.7.7
cd pysam-0.7.7
python setup.py build

# Installing pysam 0.7.7 to temporary directory
PYTHONPATH="${PYTHONPATH}:$dir/tmp/lib/python$pver/site-packages:$dir/tmp/lib64/python$pver/site-packages"
export PYTHONPATH 
python setup.py install --prefix=$dir/tmp
cd ..

# Copying pysam from temporary to final directory
if [ -e "tmp/lib/python2.7/site-packages/setuptools.pth" ]; then
		dest=$(ls -d $dir/tmp/lib/python2.7/site-packages/*/)
		cp -r $dest/* $dir/pysamdir/
fi
if [ -e "tmp/lib64/python2.7/site-packages/setuptools.pth" ]; then
		dest=$(ls -d $dir/tmp/lib64/python2.7/site-packages/*/)
		cp -r $dest/* $dir/pysamdir/
fi

# Removing temporary directories of pysam
rm -r tmp
rm -r pysam-0.7.7

# Checking if installation was successful
if [ -e "pysamdir/pysam/Pileup.py" ]
then
  message "CAVA $ver successfully installed. (Usage help: python cava.py -h)"
else
  rm -r pysamdir
  message "CAVA $ver installation failed."
fi

