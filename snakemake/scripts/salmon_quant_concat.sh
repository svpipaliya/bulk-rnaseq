# First, just run a simple ‘for loop’ that will print out the names of all the files we want to use:
for sample in `cat samples.txt`; do echo ${sample}; done

# For salmon we want to pull the TPM column (column 4) as the best estimate for transcript abundance
cd /project/6001383/pip17/sars_cov_2/data/PRJEB43380
mkdir 06-Counts
mkdir 06-Counts/tmp
for sample in `cat samples.txt`; do \
    echo ${sample}
    tail -n +2 04_salmon_alignment_quant/${sample}.out/${sample}.quant.sf | cut -f4 > 06-Counts/tmp/${sample}.count
done
 
# Next, we need to get the columns for the final table. Because all of these files are sorted in the exact same order (by gene ID), we can just use the columns from any of the files
tail -n +2 04_salmon_alignment_quant/${sample}.out/quant.genes.sf | cut -f1 > 06-Counts/tmp/geneids.txt
head 06-Counts/tmp/geneids.txt

# We want to combine all of these columns together using the ‘paste’ command, and put it in a temporary file:
paste 06-Counts/tmp/geneids.txt 06-Counts/tmp/*.count > 06-Counts/tmp/tmp.out
 
# Create a header of sample names and combine it with the temp file. Samples.txt files is taken and pipe that to the sort (to ensure they are in the same order) and then 'pase' command with the '-s' option, which takes a column of values and transposes them into a row, seperated by the tab character
cat <(cat samples.txt | sort | paste -s) 06-Counts/tmp/tmp.out > 06-Counts/rnaseq_salmon_counts.txt
rm -rf 06-Counts/tmp
head 06-Counts/salmon_combined_counts.txt

# End
 
