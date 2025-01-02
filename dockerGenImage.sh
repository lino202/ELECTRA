#!/bin/bash

# Save for checking if we are on ELECTRA root dir
current_directory=$(basename "$PWD")
desired_directory="ELECTRA"

# Read image version number from user
echo
echo \*\*\*\* ELECTRA Docker Image Generator \*\*\*\*
echo
echo "Hello $USER and welcome to to the ELECTRA docker image generator."

if [ "$current_directory" != "$desired_directory" ];then
    echo "WARNING: Please execute this command from the ELECTRA root folder and not from ${current_directory}"
fi

echo "Please provide the desired version number in the format major.minor.tweak \(eg., 1.0.0\)"
read -p 'ELECTRA version (x.x.x): ' version
echo
echo "Generating Docker image: electra-docker:${version} ..."

# Create docker image
mkdir images
docker image build -t electra-docker:$version .
docker save electra-docker:$version | gzip > ./images/electra-docker_$version.tar.gz

echo
echo Generated successfuly Docker image: electra-docker:$version
