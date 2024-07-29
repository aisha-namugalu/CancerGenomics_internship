#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=VQSRvietnam

# Set variables

REFERENCE_GENOME="/etc/ace-data/CancerGenomicsWG/Project3/AishaNamugalu/hg38.fa"
# Step 1: Index the reference genome
#        echo "Indexing the reference genome..."
#        wa index  $REFERENCE_GENOME


        samples=$(cat $1)
# Step 3: Align reads to the reference genome
#for i in $samples
#do
#        echo "Aligning reads to the reference genome for $i"
#        bwa mem -R "@RG\tID:$i\tSM:$i\tPL:ILLUMINA"\
#	$REFERENCE_GENOME Vietnam_Samples/${i}_1.fastq Vietnam_Samples/${i}_2.fastq | samtools view -Shb >${i}.bam
#done

# Step 4: Sort and indexing the aligned BAM file
#for i in $samples
#do
#      echo "Sorting the aligned BAM file for $i"
#      samtools sort ${i}.bam -o ${i}_sorted.bam
#      samtools index ${i}_sorted.bam

#done

# Step 5: Mark Duplicates
#for i in $samples
#do
 #   echo "Marking duplicates for $i"
  #  gatk MarkDuplicates \
   #     -I ${i}_sorted.bam \
    #    -O ${i}_marked.bam \
     #   -M ${i}_marked_dup_metrics.txt \
      #  --REMOVE_DUPLICATES false \
       # --CREATE_INDEX true
#done
#create sequence dictionary
#gatk CreateSequenceDictionary \
#--REFERENCE $REFERENCE_GENOME \
#--OUTPUT /etc/ace-data/CancerGenomicsWG/Project3/AishaNamugalu/hg38.fa.dict

#Regenerating index for the anotation vcf
#mkdir -p ~/gatk_temp && gatk IndexFeatureFile \
#         --input /etc/ace-data/CancerGenomicsWG/VariantCalling/GATK/reference/All_20180418.vcf.gz \
#         --tmp-dir  ~/gatk_temp
# Step 9: Base Quality Score Recalibration (BQSR)
k1="/etc/ace-data/CancerGenomicsWG/VariantCalling/GATK/reference/1000G_phase1.snps.high_confidence.hg38.vcf"
k2="/etc/ace-data/CancerGenomicsWG/VariantCalling/GATK/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf"
k3="/etc/ace-data/CancerGenomicsWG/VariantCalling/GATK/reference/Homo_sapiens_assembly38.known_indels.vcf"
#mkdir -p temp
# Perform BQSR for each sample
#for i in $samples
#do
#	echo "Performing Base Quality Score Recalibration (BQSR) for $i"
#    gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' BaseRecalibrator \
#        -I ${i}_sorted.bam \
#        -R $REFERENCE_GENOME \
#        --known-sites $k1 \
#        --known-sites $k2 \
#        --known-sites $k3 \
#        -O ${i}_bqsr.table
# done
# Apply BQSR for each sample
#for i in $samples
#do
#    echo "Applying BQSR for $i"
#    gatk ApplyBQSR \
#        -R $REFERENCE_GENOME \
 #       -I ${i}_sorted.bam \
 #       --bqsr-recal-file ${i}_bqsr.table \
#        -O ${i}_recal.bam
#done

# Build BAM index for each sample
#for i in $samples
#do
#    echo "Building BAM index for $i"
#    gatk BuildBamIndex \
#        -I ${i}_recal.bam
#done
	

# Step 11: gVCF Calling
#for i in $samples
#do
#        echo "Calling gVCF for $i"
#        gatk HaplotypeCaller \
#        -R $REFERENCE_GENOME \
#        -I ${i}_recal.bam \
#        -O ${i}_g.vcf \
#        -ERC GVCF --dbsnp $k2
#done

# Step 12: Genotype gVCF
#for i in $samples
#do
#        echo "Genotyping gVCF for $i"
#        gatk GenotypeGVCFs \
#        -R $REFERENCE_GENOME \
#        -V ${i}_g.vcf \
#        -O ${i}_g.vcf
#done
# Step 13: Variant Quality Score Recalibration (VQSR)
for i in $samples
do
	echo "Performing Variant Quality Score Recalibration (VQSR) for $i"
	gatk VariantRecalibrator \
	-R $REFERENCE_GENOME \
	-V ${i}_g.vcf \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $k1 \
        --resource:omni,known=false,training=true,truth=true,prior=12.0 $k2 \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 $k3 \	
        -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
        -mode SNP \
        -O ${i}_output.recal \
        --tranches-file ${i}_output.tranches \
        --rscript-file ${i}_output.plots.R
# Check if the .recal file was created
	if [ ! -f "${i}_output.recal" ]; then
        echo "Error: ${i}_output.recal was not created."
        exit 1
    fi
	gatk ApplyVQSR \
        -R $REFERENCE_GENOME \
        -V ${i}_g.vcf \
        --recal-file ${i}_output.recal \
        --tranches-file ${i}_output.tranches \
        -mode SNP \
        -O ${i}_RECALIBRATED_VCF
done

# Step 14: ANNOVAR Annotation
  #  echo "Annotating with ANNOVAR for $SAMPLE_NAME..."
 #   table_annovar.pl $RECALIBRATED_VCF /path/to/annovar/humandb/ \
#       # -buildver hg19 \
      #  -out $ANNOTATED_VCF \
     #   -remove \
    #    -protocol refGene,cytoBand,exac03,avsnp150,dbnsfp35a \
    #    -operation g,r,f,f,f \
   #     -nastring . \
  #      -vcfinput

    # Step 15: Filter Variants
 #   echo "Filtering variants for $SAMPLE_NAME..."
#    gatk VariantFiltration \
     #   -R $REFERENCE_GENOME \
    #    -V $ANNOTATED_VCF \
   #     --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0" \
  #      --filter-name "basic_snp_filter" \
 #       -O $FILTERED_VCF

    # Step 16: Haplotype Calling
#    echo "Calling haplotypes for $SAMPLE_NAME..."
#    gatk HaplotypeCaller \
   #     -R $REFERENCE_GENOME \
  #      -I $BQSR_BAM \
 #       -O $FILTERED_VCF
#
 #   echo "Processing for sample $SAMPLE_NAME completed successfully."
#done

#echo "Variant calling pipeline for all samples completed successfully."
#done
