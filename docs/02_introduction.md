In brief this workflow can be used to perform: 

+ Basecalling of a directory of pod5 or fast5 signal data
+ Basecalling in Duplex mode
+ Modified basecalling
+ Basecalling in real time
+ Output basecalled sequences in various formats: FASTQ, CRAM or Unaligned BAM
+ If a reference is provided a sorted and indexed BAM or CRAM will be output
for basecalling a directory of `pod5` or `fast5` signal data with `dorado`
and aligning it with `minimap2` to produce a sorted, indexed CRAM.
