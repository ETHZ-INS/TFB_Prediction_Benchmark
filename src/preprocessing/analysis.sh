#!/bin/bash
datadir="$1"
# H1: adapters=/reference/adapters/truseqPE.fa
adapters="$2"
outdir="$3"
ref="$4"
chrSizes="$5"
base="$6"

mkdir -p $datadir/trimmed
mkdir -p $datadir/aligned

bam=$datadir/aligned/"$base".bam
trimdir=$datadir/trimmed
alignedDir=$datadir/aligned


#col_num=$(head -n 1 "$file" | tr ',' '\n' | grep -n "^Name$" | cut -d: -f1)
#awk -F',' -v col="$col_num" 'NR>1 {print $col}' "$file" | while read -r entry; do

dir_pattern="SRX*"
find $datadir -name "$dir_pattern" -type d | while read -r entry; do
   for file in "$entry"/*_1.fastq.gz; do
      if [[ -f "$file" ]]; then
          sraRunName=$(basename "$file" "_1.fastq.gz")
          echo $sraRunName
          echo "processing ${sraRunName}"
          
          trimmomatic PE -threads 6 -summary $trimdir/"$base".stats -phred33 $entry/${sraRunName}_1.fastq.gz $entry/${sraRunName}_2.fastq.gz $trimdir/${sraRunName}_1.paired.fastq.gz $trimdir/${sraRunName}_1.unpaired.fastq.gz $trimdir/${sraRunName}_2.paired.fastq.gz $trimdir/${sraRunName}_2.unpaired.fastq.gz ILLUMINACLIP:$adapters:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25
          (bowtie2 -p 12 --dovetail --no-mixed --no-discordant -I 15 -X 2000 -x $ref -1 $trimdir/"${sraRunName}"_1.paired.fastq.gz -2 $trimdir/"${sraRunName}"_2.paired.fastq.gz) 2> $alignedDir/"${sraRunName}".bowtie2 | samtools view -bS - | samtools sort -@4 -m 2G - > $alignedDir/${sraRunName}.bam
          java -jar /common/picard.jar MarkDuplicates I=$alignedDir/${sraRunName}.bam O=$alignedDir/filtered_${sraRunName}.bam.2 M=$alignedDir/${sraRunName}.picard.dupMetrics.txt
          samtools sort $alignedDir/${sraRunName}.bam -o $alignedDir/sorted_${sraRunName}.bam #check if filtered file should be used!!!
          
          ##samtools index $bam && rm $trimdir/"$base"_*.fastq.gz
    fi
    done
    
done

echo "Merging bam-files"
bam_files=($alignedDir/sorted_*.bam)
samtools merge -o $datadir/merged_${base}.bam "${bam_files[@]}" -f
    
mv $datadir/merged_${base}.bam ${outdir}/rawMerged_${base}.bam
samtools sort -@ 4 ${outdir}/rawMerged_${base}.bam -o ${outdir}/sorted_rawMerged_${base}.bam
rm ${outdir}/rawMerged_${base}.bam

echo "Subsetting bam-files to chr1-chr22"   
samtools index ${outdir}/sorted_rawMerged_${base}.bam

# subset to autosomes
samtools view -b -o ${outdir}/merged_${base}.bam ${outdir}/sorted_rawMerged_${base}.bam chr{1..22}
   

# # 1. convert to bedpe
# bedtools bamtobed -bedpe -i ${outdir}/"rawMerged_${base}.bam" > ${outdir}/"merged_${base}.bedpe"
# rm ${outdir}/"rawMerged_${base}.bam"
#     
# # 3. subset to chromosomes 1-22
# awk '($1 ~ /^[0-9]+$/ && $1 >= 1 && $1 <= 22)' ${outdir}/"merged_${base}.bedpe" > ${outdir}/"merged_sub_${base}.bedpe"
# rm ${outdir}/"merged_${base}.bedpe"
#     
# echo "conversion to bedpe successful"
#     
# # 4. shift +4/-4
# awk -F'\t' -v OFS='\t' '{{$2 = $2 + 4; $3 = $3 - 4; print}}' ${outdir}/"merged_sub_${base}.bedpe" > ${outdir}/"merged_shifted_${base}.bedpe"
# rm ${outdir}/"merged_sub_${base}.bedpe"
# awk -F'\t' -v OFS='\t' '{{$5 = $5 + 4; $6 = $6 - 4; print}}' ${outdir}/"merged_shifted_${base}.bedpe" > ${outdir}/"merged_shifted2_${base}.bedpe"
# rm ${outdir}/"merged_shifted_${base}.bedpe"
# 
# # switch chromosome naming convention
# awk 'BEGIN{OFS="\t"} {$1="chr"$1; print}' ${outdir}/"merged_shifted2_${base}.bedpe" > ${outdir}/"merged_chr_${base}.bedpe"
# rm ${outdir}/"merged_shifted2_${base}.bedpe"    
# awk 'BEGIN{OFS="\t"} {$4="chr"$4; print}' ${outdir}/"merged_chr_${base}.bedpe" > ${outdir}/"merged_chr2_${base}.bedpe"
# rm ${outdir}/"merged_chr_${base}.bedpe"
#     
# # 5. convert back to bam
# bedtools bedpetobam -i ${outdir}/"merged_chr2_${base}.bedpe" -g $chrSizes > ${outdir}/"merged_${base}.bam"

