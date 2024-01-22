Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | wf-basecalling-report.html | Report summarising the work done by the basecalling workflow | per-sample |
| Simplex alignment file of passed reads | {{ alias }}.pass.simplex.{{ format }} | BAM or CRAM file of simplex reads for the sample that pass QC filtering. | per-sample |
| Duplex alignment file of passed reads | {{ alias }}.pass.duplex.{{ format }} | BAM or CRAM file of duplex reads for the sample that pass QC filtering. Created if duplex basecalling is requested. | per-sample |
| Simplex alignment file index of passed reads | {{ alias }}.pass.simplex.{{ format }}.{{ index_format }} | The index of the resulting BAM or CRAM file with the simplex reads that pass QC filtering. | per-sample |
| Duplex alignment file index of passed reads | {{ alias }}.pass.duplex.{{ format }}.{{ index_format }} | The index of the resulting BAM or CRAM file with the duplex reads that pass QC filtering. Created if duplex basecalling is requested. | per-sample |
| Simplex alignment file of failed reads | {{ alias }}.fail.simplex.{{ format }} | BAM or CRAM file of simplex reads for the sample that fail QC filtering. | per-sample |
| Duplex alignment file of failed reads | {{ alias }}.fail.duplex.{{ format }} | BAM or CRAM file of duplex reads for the sample that fail QC filtering. Created if duplex basecalling is requested. | per-sample |
| Simplex alignment file index of failed reads | {{ alias }}.fail.simplex.{{ format }}.{{ index_format }} | The index of the resulting BAM or CRAM file with the simplex reads that fail QC filtering. | per-sample |
| Duplex alignment file index of failed reads | {{ alias }}.fail.duplex.{{ format }}.{{ index_format }} | The index of the resulting BAM or CRAM file with the duplex reads that fail QC filtering. Created if duplex basecalling is requested. | per-sample |
