# Transcriptional plasticity of mesophotic corals among natural populations and transplants of Montastraea cavernosa in the Gulf of Mexico and Belize

This is a collection of bioinformatics pipelines and scripts originally developed by Misha Matz (matz@utexas.edu) for tag-based RNA-seq, and modified by Michael Studivan (studivanms@gmail.com) for use with Eli Meyer's lab protocol.

#------------------------------
countreads_raw.pl
	Recognizes .fq files instead of .fastq
	
countreads_trim.pl
	Recognizes .trim files instead of .fq
	Note: will also produce "read count" values for the job submission scripts with 'trim' in the filename

#------------------------------
rnaseq_clipper_MS.pl
	Modified with the adaptor sequence used (NNNNGGG)
	Removes PCR duplicates sharing same degenerate header and transcript sequence
	Please note this is a rough draft and will be cleaned up in time
