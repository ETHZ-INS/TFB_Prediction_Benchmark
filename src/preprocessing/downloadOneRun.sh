#!/bin/bash
# Pierre-Luc Germain, 28.11.2021, released under GPL-3

if [ $# -eq 0 ] || [ $# -lt 2 ]; then
  echo "Downloads a SRA run specified by a SRR* id,"
  echo "first trying via EBI, otherwise fasterq-dump."
  echo "Works on SE/PE data, outputs fastq.gz files."
  echo "(curl, pigz and fasterq-dump are required)"
  echo "Usage:"
  echo "`basename $0` SRRid [output_folder]"
  exit
fi

id=$1
outdir=$2
if [ -z "$2" ]; then
  outdir=$id
fi
ncpus=4
nretry=25

# build the european base link:
baseUrl="ftp://ftp.sra.ebi.ac.uk/vol1/fastq"
if [ ${#id} -ge 10 ]; then
  baseLink=$baseUrl/`echo $id | cut -c 1-6`/00${id: -1}/$id/$id
else
  baseLink=$baseUrl/`echo $id | cut -c 1-6`/$id/$id
fi

SECONDS=0
# first check if the file(s) exist(s) on the (faster) ENA server
# check for PE files:
if curl -s --retry 1 -o /dev/null --head --fail "$baseLink"_1.fastq.gz ; then
  printf "Downloading $id (PE) via curl-EBI: "
  curl -C - --retry $nretry -s -o $outdir/$id"_1.fastq.gz" "$baseLink"_1.fastq.gz & \
  curl -C - --retry $nretry -s -o $outdir/$id"_2.fastq.gz" "$baseLink"_2.fastq.gz & wait
  # check for an eventual third file:
  if curl -s --retry 1 -o /dev/null --head --fail "$baseLink"_3.fastq.gz ; then
    curl -C - --retry $nretry -s -o $outdir/$id"_3.fastq.gz" "$baseLink"_3.fastq.gz
  fi
  # if nothing found, check for SE file:
elif curl -s --retry 1 -o /dev/null --head --fail "$baseLink.fastq.gz" ; then
  printf "Downloading $id (SE) via curl-EBI: "
  curl -C - --retry $nretry -s -o $outdir/$id.fastq.gz $baseLink.fastq.gz
  # not on ENA, resort to SRA:
else
  echo "Downloading $id via fasterq-dump: "
  mkdir -p $outdir/tmp
  fasterq-dump $id -e $ncpus --temp $outdir/tmp --skip-technical --split-files -p -O $outdir && \
    pigz -p $ncpus $outdir/*.fastq
  if [ -d "$outdir/tmp" ]; then rm -r "$outdir/tmp"; fi
fi

if !(ls $outdir/$id*.fastq.gz 1> /dev/null 2>&1); then
  if [ -z "$3" ]; then
    echo "Something went wrong, retrying..."
    sleep 1
    /bin/sh $0 $1 "$2" noretry
  else
    echo $id >> $outdir/failed.txt
  fi
fi

if [ -z "$3" ]; then 
  echo "done in $(($SECONDS/60))min$(($SECONDS % 60))s"
fi


