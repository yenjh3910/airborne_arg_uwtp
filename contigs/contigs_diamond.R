# contigs_diamond

tmp <- read.table("../../airborne_arg_uwtp_result/contigs_diamond/SARG/ARP1_contigs.SARG.dmnd")
colnames(tmp) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                   "qstart","qend","sstart","send","evalue","bitscore") 
