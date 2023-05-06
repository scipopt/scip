#! /bin/bash




mkdir -p instances/miplib2017/

# download the archive benchmark.zip
if ! [ -f instances/benchmark.zip ]
then
    echo Downloading benchmark instances
    wget http://miplib.zib.de/downloads/benchmark.zip -P instances/
else
    echo Benchmark instance archive already exists
fi

unzip -u instances/benchmark.zip -d instances/miplib2017

echo Instances can be found in instances/miplib2017/