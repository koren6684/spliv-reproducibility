# Koren (2018) data requirements

## Citation and acquisition

The replication targets the Koren (2018) American Journal of Agricultural
Economics application on crop yields and conflict (full bibliographic title,
DOI, and the author-confirmed public data URL are at https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/Q3UISS). The restricted
Stata panel file must be obtained from the original study's replication archive
or another author-authorized source; it is not redistributed here.

## Restrictions

No third-party data are included in this repository. Redistribution rights and
the correct DOI/acquisition URL must be confirmed by the author before release.

## Required file and layout

Set `SPLIV_KOREN_DATA` to the full path of:

`crop.dat.af.rnr3.dta`

The file may be anywhere on the local machine. The runner does not search old
absolute project paths or substitute another file.

## Required variables

`acled_inc_sum`, `maize_yield`, `wheat_yield`, `spi6`, `gid`, `year`, and
`sparsebare` are required by the staged Table 3 workflow.

## Expected sample sizes

The source reference table records 72,169 complete observations and 6,680
`gid` clusters for both maize and wheat Table 3 specifications. The check
script reports actual counts and should be treated as authoritative if the
authorized source revision differs.

## Checksums

The data are not staged. A locally available source file observed during audit
had SHA-256 `2c69c63c94f3a338db1f77d99155c313d79022f3319ee6385cf7739a17e7ade5`;
the author should confirm that this checksum corresponds to the redistributable
archive before relying on it.

## Validate local data

From this directory:

```bash
SPLIV_KOREN_DATA=/absolute/path/to/crop.dat.af.rnr3.dta Rscript check_data.R
```
