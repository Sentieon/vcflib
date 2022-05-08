# Sentieon's vcflib
An open-source Python library for parsing and manipulation of variant call
format (VCF) files.

## Highlights
- Transparently supports reading and writing uncompressed or bgzip-compressed
  and tabix-indexed VCF files.
- VCF and gVCF (with the `.g.vcf.gz` suffix) files are supported.
- Output files write tribble (`.vcf.idx`) and tabix (`vcf.gz.tbi`) indexes on
  the fly, without a separate pass through the data for indexing.
- Parallelization across genomic regions is supported with a `Sharder` class.
- Support for Python2.7 and Python3.

## Example Usuage - VCF filtering
A simple script that filters variants with a DP < 10 from an input VCF is
provided as an example at `example/filter_dp.py`.

```
PYTHONPATH=$(pwd) python example/filter_dp.py --input_vcf <VCF> --output_vcf <VCF>
```

## Contributing
Thank you for your interest in contributing to vcflib, we welcome community
contributions! The source files for vcflib are currently maintained in an
internal codebase, and pull requests cannot be merged directly, although they
will be reviewed by our team.
