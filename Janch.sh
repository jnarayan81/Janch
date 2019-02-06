#!/bin/bash

echo -e "This is a bash scrip to check the clipping using BAM and REFERENCE file !"

#
# split the BAM fil
# /home/urbe/Tools/bamtools/bin/bamtools split -in improved3.fasta.sorted.bam -reference
# bedtools genomecov -ibam improved3.fasta.sorted.bam -d > cov.txt

#USAGE: ./Janch.sh -k yes -d /media/urbe/MyCDrive/JitDATA/FALCON -b /media/urbe/MyCDrive/JitDATA/FALCON/improved3.fasta.sorted.bam -g path/refGenome threshold

#Set the location
#Location of the tools
R=/home/urbe/anaconda3/bin/R
samtools=/home/urbe/anaconda3/bin/samtools
bamtools=/home/urbe/Tools/bamtools/bin/bamtools
bedtools=/home/urbe/anaconda3/bin/bedtools
#Location of scripts
perlScript=./scriptBase

Red=`tput setaf 1`
Green=`tput setaf 2`
Reset=`tput sgr0`

# location to getopts.sh file
source ./scriptBase/getopt.sh
#source scriptBase/*

USAGE="-k KEEP -d DIRECTORY -g GENOME -b BAM -t THRESHOLD -w WINDOW [-a START_DATE_TIME ]"
parse_options "${USAGE}" ${@}

echo "${Green}--:LOCATIONS:--${Reset}"
echo "${Green}Directory set:${Reset} ${DIRECTORY}"
echo "${Green}Genome name provided:${Reset} ${GENOME}"
echo "${Green}BAM file used:${Reset} ${BAM}"

#Parameters accepted -- write absolute path of the BAM file
bamLoc=${BAM}
ref=${GENOME}
filter=${THRESHOLD}
bin=${WINDOW}

bamFile=$(basename "$bamLoc" .bam)
dir=${bamLoc%/*}
dir=${DIRECTORY} #Stote the location

echo -e "\n${Green}Create a genome size file with${Reset} $ref"
$samtools faidx $dir/$ref
$samtools faidx $dir/$ref | cut -f 1,2 $dir/$ref.fai > $dir/$ref.sizes

echo -e "\n${Green}Make window size${Reset} $bin ${Green} for${Reset} $ref ${Green}genome${Reset}"
$bedtools makewindows -g $dir/$ref.sizes -w $bin > $dir/$ref.sizes.$bin.bp.windows.bed.win

#split the BAM file
echo -e "\n${Green}Splitting the BAM file${Reset} $bamLoc"
$bamtools split -in $dir/$bamLoc -reference

#Create an array with all the filer/dir inside ~/myDir
echo -e "\n${Green}Storing all the BAMs in array${Reset}"
arr=($dir/*.bam)

echo -e "\n${Green}Working on indivisual chromosome based BAM file${Reset}"
#Iterate through array using a counter
#Run till array-1 becuase one of the them will be ref BAM
for ((i=0; i<${#arr[@]}; i++)); do
    #do something to each element of array
    fname=$(basename "${arr[$i]}")
    if [ "$fname" != "$bamFile.bam" ] #to ignore the original BAM file
    then
	echo -e "FILE: $fname\nLOCATION: ${arr[$i]}"
	echo "Working on ${arr[$i]}"
	#Create depth file
	#$bedtools genomecov -ibam ${arr[$i]} -bga > $dir/coverage.$fname.tmp
	#We can grep on the fly, no need to create this huge file
	$samtools depth -a ${arr[$i]} > $dir/coverage.$fname.tmp

	echo "Converting BAM to SAM"
	$samtools view -h -o $dir/out_$fname.sam ${arr[$i]}

	#Read SAM and Extract soft and hard clipped reads in two files "softClipped.txt" and "hardClipped.txt"
	perl $perlScript/extractClipped2.pl $dir/out_$fname.sam $dir

	#Concat both clipped format
	cat $dir/hardClipped.txt $dir/softClipped.txt > $dir/bothClipped.$fname.txt

	#Create in Bed file
	perl $perlScript/clipped2bed.pl $dir/bothClipped.$fname.txt > $dir/bothClipped.$fname.bed

	#Check the overlapps with window : Check for 0 and 1 based system
	$bedtools coverage -a $dir/$ref.sizes.$bin.bp.windows.bed.win -b $dir/bothClipped.$fname.bed > $dir/clipCov.$fname.bed

	#Extract the name of the chromosomes
	echo "Extracting chr name"
	chrName=$(awk 'NR == 1 {print $1}' $dir/coverage.$fname.tmp);

	#Grep the chromosome - this is not needed as we are dealing chromosome wise
	more $dir/clipCov.$fname.bed | grep "$chrName" > $dir/clipCov2.$chrName.bed

	#Ignore with 0 / perfect region / no clipped region
	awk '$7 != 0 { print }' $dir/clipCov2.$chrName.bed  > $dir/clipCov.$chrName.clean.bed

	#Extract the chr of interest
	cat $dir/coverage.$fname.tmp | grep "$chrName" > $dir/cov.$chrName.bed

	#score the clipping
	perl $perlScript/covScore.pl $dir/cov.$chrName.bed $dir/clipCov.$chrName.clean.bed > $dir/final.cov.$chrName.clip

	#Plot with R
	echo "Plotting with R"
	awk {'print $1"\t"$8'} $dir/final.cov.$chrName.clip > $dir/chr.tmp
	Rscript $perlScript/plotDensity.R $dir/chr.tmp $dir/final.cov.$chrName.clip.png
	#$chrName -n $chrName

	#Remove sam
  #If file is empty
  find $dir -size  0 -print0 |xargs -0 rm --
	rm -rf $dir/coverage.$fname.tmp
	rm -rf $dir/out_$fname.sam
	rm -rf $dir/*.bed
	rm -rf $dir/*.txt
	rm -rf ${arr[$i]}
	    fi
    done

cat $dir/*.clip > final.all.clip
awk {'print $1"\t"$8'} final.all.clip > all.clip
Rscript $perlScript/plotDensity.R all.clip all.png

echo "Janch DONE ..."
