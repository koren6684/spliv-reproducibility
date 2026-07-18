# Lelkes, Sood, and Iyengar (2017) data requirements

## Citation and acquisition

The replication targets the Lelkes, Sood, and Iyengar (2017) broadband and
affective-polarization application (full bibliographic title, DOI, and the
author-confirmed public data URL are at https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/W1YJZQ). Obtain the original study's
authorized replication archive or a provider-approved archive; it is not
redistributed here.

## Restrictions

No third-party data are included. The source archive may contain individual,
county, provider, and demographic files whose redistribution terms must be
checked with the original authors/data provider before release.

## Required layout and filenames

Set `SPLIV_LELKES_DATA` to a directory containing both:

- `mergeddataset.RData` (must contain object `merged`)
- `mergedcounty.RData` (must contain object `countyproviders`)

The runner does not search legacy absolute project paths or substitute data.

## Required variables

The merged object must contain `infeels`, `outfeels`, `providers`, `Total`,
`year`, `region`, `percent_black`, `percent_white`, `percent_male`, `lowed`,
`unemploymentrate`, `density`, `HHINC`, and `state`.

## Expected sample size

The source reference table records 114,803 complete Table 1 observations and
48 state clusters. The check script reports actual counts for the authorized
archive and takes precedence if the source revision differs.

## Checksums

`mergeddataset.RData` was not available in the local audit tree. A local
`mergedcounty.RData` had SHA-256
`153560ef3a30236b69e2ed8728d99c38238e479ece06f2ac46ee597cf2b21301`; confirm
that checksum against the authorized archive before use.

## Validate local data

From this directory:

```bash
SPLIV_LELKES_DATA=/absolute/path/to/lelkes-data-directory Rscript check_data.R
```
