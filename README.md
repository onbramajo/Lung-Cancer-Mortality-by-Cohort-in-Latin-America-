Birth-cohort patterns of lung cancer mortality in seven Latin American countries

Code for the paper "Birth-cohort patterns of lung cancer mortality in seven Latin American countries: a cross-national analysis" (Octavio Bramajo).

The script lung_cancer_cohorts_latam.R reproduces all results, figures and tables in the paper and supplement.

Data

The data are not included here. Both sources are public and can be downloaded for free.

Deaths. WHO Mortality Database, raw ICD-10 data files. https://www.who.int/data/data-collection-tools/who-mortality-database Download the files Morticd10_part1 to Morticd10_part6, unzip them and put them in data/.

Population. UN World Population Prospects 2024, population by single age, males and females (estimates). https://population.un.org/wpp/downloads Put these two files in data/:

WPP2024_POP_F01_2_POPULATION_SINGLE_AGE_MALE.xlsx
WPP2024_POP_F01_3_POPULATION_SINGLE_AGE_FEMALE.xlsx
Folder structure
project/
├── lung_cancer_cohorts_latam.R
├── README.md
├── data/
│   ├── Morticd10_part1 ... Morticd10_part6
│   ├── WPP2024_POP_F01_2_POPULATION_SINGLE_AGE_MALE.xlsx
│   └── WPP2024_POP_F01_3_POPULATION_SINGLE_AGE_FEMALE.xlsx
└── output/          (created by the script)
    ├── figures/
    └── tables/
How to run

Set the working directory to the project folder and run the whole script:

r
setwd("path/to/project")
source("lung_cancer_cohorts_latam.R")

The smoothing step (section 5) is the slowest part and can take several minutes.

R packages: tidyverse, readxl, ungroup, metR, splines (part of base R), broom.

r
install.packages(c("tidyverse", "readxl", "ungroup", "metR", "broom"))
What the script does
Section	Content
1	Reads WHO deaths for C33-C34 and C76-C80, by country, year, sex and 5-year age group
2	Reads WPP 2024 population, computes death rates
3	Data quality: % of deaths from ill-defined causes (R00-R99) and % of cancer deaths coded C76-C80
4	Age-standardised death rates, ages 40+ (WHO world standard; pooled Latin American 2010 population as sensitivity)
5	Smoothing to single ages with ungroup::pclm2D and rates of mortality improvement (ROMI)
6	Tests A and B for non-linear cohort effects (data to 2019)
7	Two-factor age-period and age-cohort models: QAIC, cohort curvature, drift and period effects (all years and without 2020-2023)
8	Figures
9	Tables

Countries and first year used: Uruguay, Argentina, Chile, Colombia (1997), Mexico (1998), Cuba (2001), Brazil (2006). Last year is 2023 for all. Uruguay has no data for 2011.

Outputs

Figures (output/figures/, PNG and PDF):

File	Paper
Figure_1	Age-standardised lung cancer death rates, ages 40+
Figure_2	ROMI Lexis surfaces, lung cancer
Figure_3	Cohort curvature, lung cancer
Figure_S1	Cohort effects before removing the linear trend
Figure_S2	Cohort curvature, lung cancer + C76-C80
Figure_S3	ROMI Lexis surfaces, lung cancer + C76-C80
Figure_S4	Cohort curvature, lung cancer, without 2020-2023
Figure_S5	Period effects, lung cancer

Tables (output/tables/):

File	Paper
Table_S1.csv	Data quality
Table_S2.csv	Lung cancer deaths and age-standardised rates
Table_S3.csv	Age-period vs age-cohort models (QAIC) and drift
Table_S4.csv	Tests for non-linear cohort effects

Intermediate results (output/): age-standardised rates, ROMI values, test results, model comparison, cohort curvature, drift and period effects, as CSV files.

Contact

Octavio Bramajo
