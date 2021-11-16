#!/bin/bash

MCPATH=/home/jcascitt/rnacache
RMPATH=/home/jcascitt/eco/bin/rapmap
SALMONPATH=/home/jcascitt/opt/bin/salmon
KALLPATH=/home/jcascitt/eco/bin/kallisto
BT2PATH=/home/jcascitt/eco/bin/bowtie2
SAMPATH=/home/jcascitt/eco/bin/samtools
VMTPATH=/home/jcascitt/eco/bin/vmtouch

OUTPATH=/tmp/jcascitt/si

##################################

RNADIR=$(dirname "$(readlink -f "$0")")
cd $RNADIR

DBDIR=$(realpath ${2})

$VMTPATH -e $RNADIR $DBPATH $OUTPATH

if [ -n "$1" ] && [ "$1" = "build" ]; then
	REFPATH=$(realpath ${3})

	mkdir ${DBDIR}
	ln -s $REFPATH ${DBDIR}/ref
	
	$VMTPATH -vt $REFPATH
	{ time $MCPATH/rnacache build ${DBDIR}/16.db $REFPATH -winlen 46 -winstride 31 -sketchlen 16 2>&1 ; } 2>> ${DBDIR}/mcr16_build_time
	$VMTPATH -e ${DBDIR}/16.db

	$VMTPATH -vt $REFPATH
	{ time $RMPATH quasiindex -t $REFPATH -i ${DBDIR}/rap_idx 2>&1 ; } 2>> ${DBDIR}/rap_build_time
	$VMTPATH -e ${DBDIR}/rap_idx
	
	$VMTPATH -vt $REFPATH
	{ time $SALMONPATH index -t $REFPATH -i ${DBDIR}/sal_idx -p 128 2>&1 ; } 2>> ${DBDIR}/sal_build_time
	$VMTPATH -e ${DBDIR}/sal_idx
	
	$VMTPATH -vt $REFPATH
	{ time $KALLPATH index -i ${DBDIR}/kal_idx $REFPATH 2>&1 ; } 2>> ${DBDIR}/kal_build_time
	$VMTPATH -e ${DBDIR}/kal_idx
	
	$VMTPATH -vt $REFPATH
	{ time ${BT2PATH}-build --threads 128 $REFPATH ${DBDIR}/bt2_idx 2>&1 ; } 2>> ${DBDIR}/bt2_build_time
	$VMTPATH -e ${DBDIR}/bt2_idx*
	
	echo "Indices built: ${DBDIR}" >> info
fi

if [ -n "$1" ] && [ "$1" = "benchmark" ]; then
	i_dir=${3}

	for t in 128 64 32 16 8 4 2 1; do

		cd ${RNADIR}/${i_dir}
		mkdir ${OUTPATH}/${i_dir}

		$VMTPATH -e $RNADIR $DBPATH $OUTPATH

		$VMTPATH -vt ${DBDIR}/16.db reads*.fq
		{ time $MCPATH/rnacache query ${DBDIR}/16.db reads1.fq reads2.fq -pairfiles -insertsize 250 -sketchlen 16 -threads $t -no-map -no-summary -no-query-params 2>&1 >/dev/null ;  } 2>> mcr16_${t}_time
		$VMTPATH -e ${DBDIR}/16.db

		$VMTPATH -vt ${DBDIR}/rap_idx reads*.fq
		{ time $RMPATH quasimap -i ${DBDIR}/rap_idx -1 reads1.fq -2 reads2.fq -t $t 2>&1 >/dev/null ; } 2>> rap_${t}_time

		$VMTPATH -vt ${DBDIR}/rap_idx reads*.fq
		{ time $RMPATH quasimap -i ${DBDIR}/rap_idx -1 reads1.fq -2 reads2.fq -t $t -s 2>&1 >/dev/null ; } 2>> ras_${t}_time
		$VMTPATH -e ${DBDIR}/rap_idx

		$VMTPATH -vt ${DBDIR}/sal_idx reads*.fq
		{ time $SALMONPATH quant -i ${DBDIR}/sal_idx -l A -1 reads1.fq -2 reads2.fq -p $t -o ${OUTPATH}/${i_dir}/sal_out --skipQuant -q 2>&1 ; } 2>> sal_${t}_time
		rm -rf ${OUTPATH}/${i_dir}/sal_out
		$VMTPATH -e ${DBDIR}/sal_idx

		$VMTPATH -vt ${DBDIR}/kal_idx reads*.fq
		{ time $KALLPATH quant -i ${DBDIR}/kal_idx -o ${OUTPATH}/${i_dir}/kal_out -t $t reads1.fq reads2.fq 2>&1 ; } 2>> kal_${t}_time
		rm -rf ${OUTPATH}/${i_dir}/kal_out
		$VMTPATH -e ${DBDIR}/kal_idx
		
		if (($t >= 16)); then
			$VMTPATH -vt ${DBDIR}/bt2_idx* reads*.fq
			{ time $BT2PATH -p $t -x ${DBDIR}/bt2_idx -1 reads1.fq -2 reads2.fq 2>&1 >/dev/null ; } 2>> bt2_def_${t}_time

			$VMTPATH -vt ${DBDIR}/bt2_idx* reads*.fq
			{ time $BT2PATH -p $t -k 20 -x ${DBDIR}/bt2_idx -1 reads1.fq -2 reads2.fq 2>&1 >/dev/null ; } 2>> bt2_20_${t}_time
		
			$VMTPATH -vt ${DBDIR}/bt2_idx* reads*.fq
			{ time $BT2PATH -k 200 -p $t -x ${DBDIR}/bt2_idx -1 reads1.fq -2 reads2.fq 2>&1 >/dev/null ; } 2>> bt2_200_${t}_time
			$VMTPATH -e ${DBDIR}/bt2_idx*
		fi

		$VMTPATH -e $OUTPATH
		rm -rf ${OUTPATH}/${i_dir}

		echo "Queries completed: ${i_dir} ${t}" >> ../info
	done
fi

if [ -n "$1" ] && [ "$1" = "run" ]; then
	for i_dir in "${@:3}"; do
	
		cd ${RNADIR}/${i_dir}
		mkdir ${OUTPATH}/${i_dir}

		$VMTPATH -e $RNADIR $DBPATH $OUTPATH

		$VMTPATH -vt ${DBDIR}/16.db reads*.fq
		{ time $MCPATH/rnacache-bam query ${DBDIR}/16.db reads1.fq reads2.fq -pairfiles -insertsize 250 -sketchlen 16 -remove-overpopulated-features -with-bam-out ${OUTPATH}/${i_dir}/mcr16_bam.bam -no-map -threads 112 -bam-threads 16 2>&1 ; } 2>> mcr16_bam_time
		$VMTPATH -e ${OUTPATH}/${i_dir}/mcr16_bam.bam

		$VMTPATH -vt ${DBDIR}/16.db reads*.fq
		{ time $MCPATH/rnacache query ${DBDIR}/16.db reads1.fq reads2.fq -pairfiles -insertsize 250 -sketchlen 16 -remove-overpopulated-features -accuracy -out ${OUTPATH}/${i_dir}/mcr16_out 2>&1 ; } 2>> mcr16_time
		tail -n 14 ${OUTPATH}/${i_dir}/mcr16_out >> mcr16_eval
		$VMTPATH -e ${OUTPATH}/${i_dir}/mcr16_out

		$VMTPATH -vt ${DBDIR}/16.db reads*.fq
		{ time $MCPATH/rnacache query ${DBDIR}/16.db reads1.fq reads2.fq -pairfiles -insertsize 250 -sketchlen 16 -remove-overpopulated-features -accuracy -cov-min 0 -hit-min 1 -hit-cutoff 0 2>&1 > ${OUTPATH}/${i_dir}/mcr16_nofilter_out ; } 2>> mcr16_nofilter_time
		tail -n 12 ${OUTPATH}/${i_dir}/mcr16_nofilter_out >> mcr16_nofilter_eval
		$VMTPATH -e ${OUTPATH}/${i_dir}/mcr16_nofilter_out

		$VMTPATH -vt ${DBDIR}/16.db reads*.fq
		{ time $MCPATH/rnacache query ${DBDIR}/16.db reads1.fq reads2.fq -pairfiles -insertsize 250 -sketchlen 16 -remove-overpopulated-features -accuracy -cov-min 0 2>&1 > ${OUTPATH}/${i_dir}/mcr16_onlineonly_out ; } 2>> mcr16_onlineonly_time
		tail -n 12 ${OUTPATH}/${i_dir}/mcr16_onlineonly_out >> mcr16_onlineonly_eval
		$VMTPATH -e ${OUTPATH}/${i_dir}/mcr16_onlineonly_out ${DBDIR}/16.db 

		$VMTPATH -vt ${DBDIR}/rap_idx reads*.fq
		{ time $RMPATH quasimap -i ${DBDIR}/rap_idx -1 reads1.fq -2 reads2.fq -t 112 2>&3 3>&- | $SAMPATH view -b -@ 16 - > ${OUTPATH}/${i_dir}/rap.bam 2>&1 3>&- ; } 3>&1 2>> rap_time
		$VMTPATH -e ${OUTPATH}/${i_dir}/rap.bam

		$VMTPATH -vt ${DBDIR}/rap_idx reads*.fq
		{ time $RMPATH quasimap -i ${DBDIR}/rap_idx -1 reads1.fq -2 reads2.fq -t 112 -s 2>&3 3>&- | $SAMPATH view -b -@ 16 - > ${OUTPATH}/${i_dir}/ras.bam 2>&1 3>&- ; } 3>&1 2>> ras_time
		$VMTPATH -e ${DBDIR}/rap_idx ${OUTPATH}/${i_dir}/ras.bam

		$VMTPATH -vt ${DBDIR}/sal_idx reads*.fq
		{ time $SALMONPATH quant -i ${DBDIR}/sal_idx -l A -1 reads1.fq -2 reads2.fq -p 112 -o ${OUTPATH}/${i_dir}/sal_out -z --skipQuant 2>&3 3>&- | $SAMPATH view -b -@ 16 - > ${OUTPATH}/${i_dir}/sal.bam 2>&1 3>&- ; } 3>&1 2>> sal_time
		$VMTPATH -e ${DBDIR}/sal_idx ${OUTPATH}/${i_dir}/sal.bam

		$VMTPATH -vt ${DBDIR}/kal_idx reads*.fq
		{ time $KALLPATH quant --pseudobam -i ${DBDIR}/kal_idx -o ${OUTPATH}/${i_dir}/kal_out -t 128 reads1.fq reads2.fq 2>&1 ; } 2>> kal_time
		mv ${OUTPATH}/${i_dir}/kal_out/pseudoalignments.bam ${OUTPATH}/${i_dir}/kal.bam
		$VMTPATH -e ${DBDIR}/kal_idx ${OUTPATH}/${i_dir}/kal_out ${OUTPATH}/${i_dir}/kal.bam

		$VMTPATH -vt ${DBDIR}/bt2_idx* reads*.fq
		{ time $BT2PATH -k 200 -p 120 --reorder -x ${DBDIR}/bt2_idx -1 reads1.fq -2 reads2.fq 2>&3 3>&- | $SAMPATH view -b -@ 8 - > ${OUTPATH}/${i_dir}/bt2_200.bam 2>&1 3>&- ; } 3>&1 2>> bt2_200_time
		$VMTPATH -e ${OUTPATH}/${i_dir}/bt2_200.bam

		$VMTPATH -vt ${DBDIR}/bt2_idx* reads*.fq
		{ time $BT2PATH -k 20 -p 120 --reorder -x ${DBDIR}/bt2_idx -1 reads1.fq -2 reads2.fq 2>&3 3>&- | $SAMPATH view -b -@ 8 - > ${OUTPATH}/${i_dir}/bt2_20.bam 2>&1 3>&- ; } 3>&1 2>> bt2_20_time
		$VMTPATH -e ${OUTPATH}/${i_dir}/bt2_20.bam

		$VMTPATH -vt ${DBDIR}/bt2_idx* reads*.fq
		{ time $BT2PATH -p 120 --reorder -x ${DBDIR}/bt2_idx -1 reads1.fq -2 reads2.fq 2>&3 3>&- | $SAMPATH view -b -@ 8 - > ${OUTPATH}/${i_dir}/bt2_def.bam 2>&1 3>&- ; } 3>&1 2>> bt2_def_time
		$VMTPATH -e ${OUTPATH}/${i_dir}/bt2_def.bam ${DBDIR}/bt2_idx*

		$VMTPATH -vt ${OUTPATH}/${i_dir}/*.bam
		num_reads=$(($(wc -l reads1.fq | awk '{print $1}')/4))
		${RNADIR}/evaluatemp.py $num_reads . ${OUTPATH}/${i_dir}/*.bam

		$VMTPATH -e $OUTPATH
		rm -rf ${OUTPATH}/${i_dir}
		
		echo "${i_dir} Queries completed." >> ../info
	done
fi

if [ -n "$1" ] && [ "$1" = "real" ]; then
	for i_dir in "${@:3}"; do
	
		cd ${RNADIR}/${i_dir}
		mkdir ${OUTPATH}/${i_dir}

		$VMTPATH -e $RNADIR $DBPATH $OUTPATH

		$VMTPATH -vt ${DBDIR}/16.db reads*.fq
		{ time $MCPATH/rnacache query ${DBDIR}/16.db reads1.fq reads2.fq -pairfiles -insertsize 250 -sketchlen 16 -remove-overpopulated-features -accuracy -out ${OUTPATH}/${i_dir}/mcr16_out 2>&1 ; } 2>> mcr16_time
		tail -n 14 ${OUTPATH}/${i_dir}/mcr16_out >> mcr16_eval
		$VMTPATH -e ${OUTPATH}/${i_dir}/mcr16_out ${DBDIR}/16.db 

		$VMTPATH -vt ${DBDIR}/rap_idx reads*.fq
		{ time $RMPATH quasimap -i ${DBDIR}/rap_idx -1 reads1.fq -2 reads2.fq -t 112 2>&3 3>&- | $SAMPATH view -b -@ 16 - > ${OUTPATH}/${i_dir}/rap.bam 2>&1 3>&- ; } 3>&1 2>> rap_time
		$VMTPATH -e ${OUTPATH}/${i_dir}/rap.bam

		$VMTPATH -vt ${DBDIR}/rap_idx reads*.fq
		{ time $RMPATH quasimap -i ${DBDIR}/rap_idx -1 reads1.fq -2 reads2.fq -t 112 -s 2>&3 3>&- | $SAMPATH view -b -@ 16 - > ${OUTPATH}/${i_dir}/ras.bam 2>&1 3>&- ; } 3>&1 2>> ras_time
		$VMTPATH -e ${DBDIR}/rap_idx ${OUTPATH}/${i_dir}/ras.bam

		$VMTPATH -vt ${DBDIR}/sal_idx reads*.fq
		{ time $SALMONPATH quant -i ${DBDIR}/sal_idx -l A -1 reads1.fq -2 reads2.fq -p 112 -o ${OUTPATH}/${i_dir}/sal_out -z --skipQuant 2>&3 3>&- | $SAMPATH view -b -@ 16 - > ${OUTPATH}/${i_dir}/sal.bam 2>&1 3>&- ; } 3>&1 2>> sal_time
		$VMTPATH -e ${DBDIR}/sal_idx ${OUTPATH}/${i_dir}/sal.bam

		$VMTPATH -vt ${DBDIR}/kal_idx reads*.fq
		{ time $KALLPATH quant --pseudobam -i ${DBDIR}/kal_idx -o ${OUTPATH}/${i_dir}/kal_out -t 128 reads1.fq reads2.fq 2>&1 ; } 2>> kal_time
		mv ${OUTPATH}/${i_dir}/kal_out/pseudoalignments.bam ${OUTPATH}/${i_dir}/kal.bam
		$VMTPATH -e ${DBDIR}/kal_idx ${OUTPATH}/${i_dir}/kal_out ${OUTPATH}/${i_dir}/kal.bam

		$VMTPATH -vt ${DBDIR}/bt2_idx* reads*.fq
		{ time $BT2PATH -k 200 -p 120 --reorder -x ${DBDIR}/bt2_idx -1 reads1.fq -2 reads2.fq 2>&3 3>&- | $SAMPATH view -b -@ 8 - > ${OUTPATH}/${i_dir}/bt2_200.bam 2>&1 3>&- ; } 3>&1 2>> bt2_200_time
		$VMTPATH -e ${OUTPATH}/${i_dir}/bt2_200.bam

		$VMTPATH -vt ${DBDIR}/bt2_idx* reads*.fq
		{ time $BT2PATH -k 20 -p 120 --reorder -x ${DBDIR}/bt2_idx -1 reads1.fq -2 reads2.fq 2>&3 3>&- | $SAMPATH view -b -@ 8 - > ${OUTPATH}/${i_dir}/bt2_20.bam 2>&1 3>&- ; } 3>&1 2>> bt2_20_time
		$VMTPATH -e ${OUTPATH}/${i_dir}/bt2_20.bam

		$VMTPATH -vt ${DBDIR}/bt2_idx* reads*.fq
		{ time $BT2PATH -p 120 --reorder -x ${DBDIR}/bt2_idx -1 reads1.fq -2 reads2.fq 2>&3 3>&- | $SAMPATH view -b -@ 8 - > ${OUTPATH}/${i_dir}/bt2_def.bam 2>&1 3>&- ; } 3>&1 2>> bt2_def_time
		$VMTPATH -e ${OUTPATH}/${i_dir}/bt2_200.bam ${DBDIR}/bt2_idx*

		$VMTPATH -vt ${OUTPATH}/${i_dir}/*.bam
		num_reads=$(($(wc -l reads1.fq | awk '{print $1}')/4))
		${RNADIR}/agreement2.py $num_reads ${OUTPATH}/${i_dir}/mcr16_out ${OUTPATH}/${i_dir}/*.bam > concordance.json

		$VMTPATH -e $OUTPATH
		rm -rf ${OUTPATH}/${i_dir}
		
		echo "${i_dir} Queries completed." >> ../info
	done
fi

$VMTPATH -e $RNADIR $DBPATH $OUTPATH $REFPATH

