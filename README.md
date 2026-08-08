# HIV Notification Temporal Trend Analysis

This project contains the code to analyse trends in HIV notifications in
Australia from 2001 to 2024, and to generate the tables and figures used in the
associated paper. The analysis identifies change points in notification trends 
using joinpoint (segmented) regression, characterises CD4 count at diagnosis over time, 
and maps Primary Health Network (PHN)–level notification rates and their 
geographic inequality. Throughout, two groups are analysed separately:
notifications with a **first-ever diagnosis in Australia (FED)** and those
**previously diagnosed overseas (PDO)**.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21847208.svg)](https://doi.org/10.5281/zenodo.21847208)

## Aims

The aim of this analysis is to quantify how HIV notifications and the geographic distribution of 
notifications across Australian PHNs have changed over 2001–2024. For the total population and for key
subgroups, the code produces trend estimates, change points, PHN-level rates,
measures of geographic inequality, and a set of publication-ready tables and
figures.

## Outputs

Running the analysis reproduces the following paper outputs:

- **Figure** – "rainbow" plots of CD4 count at diagnosis over time (proportion
  in each CD4 category and monthly-average number of notifications).
- **Figure** – median CD4 count at diagnosis by year.
- **Figure** – joinpoints and segmented trends of monthly notifications.
- **Figure** – PHN-level notification-rate heatmaps (per 100,000), FED vs PDO.
- **Figure** – geographic inequality across PHNs: Gini coefficient over time and
  Lorenz curves (2001 vs 2024).
- **Tables** – trend estimates, change points, national and PHN-level rate trends.

## Maintainers and developers

1. [Rongxing Weng](https://github.com/RongxingW); ORCiD ID: [0000-0003-1792-2186](https://orcid.org/0000-0003-1792-2186)

Affiliation: The Kirby Institute, UNSW Sydney, NSW, Australia

For any inquiries, please contact Rweng@kirby.unsw.edu.au or flag an issue.
The code will be updated as required to correct issues or to improve or add
features. Please check for updated versions periodically.

## Project structure

The analysis is contained in two RMarkdown scripts. `0_Setupmodel.Rmd` creates
the project folder skeleton, and `1_timeseriesAnalysis.Rmd` runs the full
analysis top to bottom. Reusable functions live in the `code/` folder and are
sourced by the main script. Project-specific inputs and outputs live under
`projects/`.

```bash
├── projects/
│   └── 2000_2024_changing_point/   # specific project name
│       ├── data                    # model input data (NOT tracked; see Data availability)
│       │   ├── AU_pop_size_2001_2024.csv       # national population by year (ABS)
│       │   ├── SA2_pop_size_2001_2024.csv      # SA2-level population by year (ABS)
│       │   ├── phn_sa2_20250901.csv            # SA2-to-PHN concordance with allocation ratios
│       │   └── ...
│       ├── output                  # outputs: tables (.csv)
│       └── figures                 # outputs: figures (.png)
├── code/
│   ├── Joinpoint_analysis.R        # segmented/joinpoint trend regression
│   ├── Rainbow_plot_and_stats.R    # CD4-at-diagnosis "rainbow" plots (proportion + number)
│   └── PHN_timeseries.R            # PHN rates, NB rate trends, heatmaps, Gini/Lorenz
├── 0_Setupmodel.Rmd                # creates the project folder skeleton
├── 1_timeseriesAnalysis.Rmd        # main analysis (run this)
├── LICENSE
└── README.md
```

The `code/` folder contains the functions used by the main RMarkdown script;
each `.R` file begins with a header describing its purpose, inputs, and outputs.
Inputs and outputs for the analysis are stored within the
`projects/2000_2024_changing_point/` folder.

## Using the code

To use the code, clone or download this repository to a convenient location on
your computer. You will need the following software and associated packages:

1. R, a free statistical program. (Optional) RStudio, a user interface for R.
2. R packages used by the analysis:

   `tidyverse`, `segmented`, `MASS`, `patchwork`, `ggtext`, `broom`, `grid`, `gtable`


### Steps to run

1. Run `0_Setupmodel.Rmd` to create
   `projects/2000_2024_changing_point/{data, output, figures}`.
2. Place the cleaned HIV notification CSV and the population / PHN-concordance
   CSVs in the project's `data/` folder.
3. Open `1_timeseriesAnalysis.Rmd`, check the file paths in the `Initialise` and
   `Input` chunks (in particular `HIVDataFolder`, which points to a local copy
   of the restricted data), and run the chunks top to bottom.

## Data availability

The HIV notification data used by this analysis are **not** publicly available
and are **not included** in this repository, due to privacy and ethics
restrictions. Only aggregate outputs (tables and figures) are
shareable. [Population denominators](https://www.abs.gov.au/statistics/people/population/regional-population/latest-release) and the [SA2-to-PHN concordance](https://www.health.gov.au/resources/publications/primary-health-networks-phn-2023-statistical-area-level-2-2021?language=en) are derived 
from publicly available data from the Australian Bureau of Statistics (ABS) and 
the Australian Government Department of Health, Disability and Ageing.

## Publication

The following publication is associated with this project and used the code in
this repository to generate the results and figures.

(...)

## License

Code in this repository is released under the MIT License (see `LICENSE`).

## Disclaimer

The code has been made publicly available for transparency and replication
purposes and in the hope that it will be useful. We take no responsibility for
results generated with the code or their interpretation, but are happy to assist
with its use and application.
