# Input data layout

The pipeline expects input data here by default:

```text
data/
  clusters.csv
  id_lists.csv
  proteins/
    proteins_visit_0.fst
  metabolites/
    nmr_threephases.csv
```

These files are intentionally ignored by git. They may contain restricted UK
Biobank or participant-level data and should not be pushed to GitHub.

You can also keep the files elsewhere and pass paths with command-line options
or environment variables; see the main `README.md`.

