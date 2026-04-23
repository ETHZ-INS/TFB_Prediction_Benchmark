#!/bin/bash
# Pierre-Luc Germain, 28.11.2021, released under GPL-3

runtable=${1:-SraRunTable.txt}
outdir="$2"

if [ ! -r "$runtable" ]; then
  echo "Could not read run table."
  echo "Usage: `basename $0` SraRunTable.txt [output_directory]"
  echo "This downloads the SRA runs specified in the SraRunTable.txt."
  echo "To get such a table, use the SRA run selector:"
  echo "https://www.ncbi.nlm.nih.gov/Traces/study/"
  echo "(See https://github.com/plger/SRAdownload )"
  exit
fi

if [ ! -w "$outdir" ]; then
  echo "Target output directory ($outdir) is not writable!"
  exit
fi

# we first group SRR runs by experiment:
tmpfile=$outdir/$(mktemp runsByExp.XXXXXX)
awk -vFPAT='([^,]*)|("[^"]+")' -vOFS=, '
NR==1 {
    expi = 1
    geo = ""
    for (i=1; i<=NF; i++){
        if($i=="Experiment") expi=i
        if($i=="GEO_Accession (exp)") geo=i
    }
    next
};
{
    if(geo!=""){
        print $expi"."$geo"\t"$1
    }else{
        print $expi"\t"$1
    }
}' $runtable | sort | awk '
$1==x{
    printf ",%s", $2
    next
}
{
    if(NR>1) printf "\n"
    printf $1"\t"$2
    x=$1
}
END {
    printf "\n"
}' > $tmpfile

# we next proceed to download them
while read -u 9 -r name runs
do
  # for each SRX/GSM ID
  # first check if a folder already exists, if so skip
  if [[ -d "$outdir/$name" ]] ; then
    echo "$name exists, skipping..."
  else
    echo $name
    mkdir "$outdir/$name"
    for SRR in `echo $runs | tr "," "\n"`; do
      /bin/bash `dirname $0`/downloadOneRun.sh "$SRR" "$outdir/$name"
    done
  fi
done 9< $tmpfile

rm $tmpfile

echo "Done!"
