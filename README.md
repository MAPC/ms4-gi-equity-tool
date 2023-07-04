ms4-gi-equity-tool
==============================



Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │   └── 1-run-ms4-parcel-model.ipynb    <- runs src > features > parcel_model.py
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data
    │   │   └── make_dataset.py <- Sets data paths, variables, and fields
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   └── suitability_criteria.py <- basic functions for "suitability" modeling
    │   │   └── ms4_funcs.py            <- specific functions for this project
    │   │   └── ms4-comm-vis.R          <- R script for generating  "community visibility" layer
    │   │   └── parcel_model.py         <- final python script for creating data for each muni
    │   │
    │   ├── models         <- Scripts generated via modelbuilder from ArcGIS Pro
    │   │   ├── pler_calc.py            <- calculates phosphorus load estimate based on land cover
    │   │   └── row_segmentation.py     <- generates right-of-way segments
    │   │
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io


--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
