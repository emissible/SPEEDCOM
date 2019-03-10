#!/bin/bash

# Download the PHOTOCHEMCAD3 database for the user, will unpack it into
# the DATA directory (not tracked by github)

URL="http://www.photochemcad.com/download/PhotochemCAD3.zip"
FLNM=$(basename $URL)
PREFIX=DATA

#make data directory if needed
ls $PREFIX 2>/dev/null
if [[ $? -ne 0 ]]; then
  echo "Making $PREFIX directory in $PWD"
  mkdir $PREFIX || \
  { echo "Unable to create directory, check permissions" && exit 1; }
else
  echo "DATA directory already exists, be careful of existing files"
fi

#Downloading the *.zip archive
which wget 2>/dev/null
if [[ $? -eq 0 ]]; then
  wget -P $PREFIX/ $URL ||\
  { echo "Unable to download, check internet and try again" && exit 1; }
else
  curl -o $PREFIX/$FLNM $URL ||\
  { echo "Unable to download, check curl installed and internet" && exit 1; }
fi

#Unzipping the *.zip archive
which unzip 2>/dev/null
if [[ $? -eq 0 ]]; then
  unzip $PREFIX/$FLNM -d DATA || \
  { echo "Unable to extract, check internet and is zip" && exit 1; }
else
  echo "You need to install 'unzip' for this script"
  exit 1
fi

#Let user know finished
echo "PHOTOCHEM archive sucessfully downloaded and unziped"
exit 0
