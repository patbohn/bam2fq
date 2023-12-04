# bam2fq

This is a minimal tool to convert bam output from dorado basecalling to fastq.gz files, separating simplex from duplex files. 

(note that in this version simplex files will contain those that generated a duplex offspring, i.e. -1)

Since the reads are loaded into memory, I also decided to calculate some basic statistics that are written to csv(.gz) files. 
Specifically, the tool generates a per-read statistics file that contains the read_id, read length, mean qscore (calculated properly by averaging mean error likelihood), 
their duplex state and the estimated poly A/T tag introduced in dorado 0.4.3. Another file contains the per base Qscore distribution. 
