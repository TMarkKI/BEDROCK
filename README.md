BEDROCK - BEDmethyl Recognition, Observation, and Characterisation Kit

Command to run the programme is as follows: python <path/to/BEDROCK/main.py> --spreadsheet <path/to/spreadsheet/tsv/format> --chr-list <path/to/chromosome/list> --ref <path/to/reference/file> --fai <path/to/reference/index> --annotation <path/to/gene/annotation> --outdir <output/path>

You can retrieve a chromosome list using: samtools idxstats <path/to/bam> | awk '$3 > 0 {print $1}'

Required sample spreadsheet with the following columns:
1. Sample Name
2. Sample Type (Can be the species of the queried sample)
3. Absolute path to Bedmethyl file
4. Absolute path to the depth file
