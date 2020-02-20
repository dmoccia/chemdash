# chemdash
Chemical Dataset Exploration Tool Built on Plotly Dash

**This README is under construction so use at your own risk**

The chemdash app is a simple analysis application built upon the plotly dash framework.  As input the application takes a _csv_ file requiring both **smiles** and **Compound_id** fields, other fields will be retained and displayed in the data table portion of the application.  A UMAP (similar to t-SNE/PCA) plot will be generated as well as a Molecular Grid Map (MGM).   Depending on the column selection in the data table, a distribution plot for the given column will also be displayed.  Selection of data points in either the UMAP or MGM plot will result in filtering of the data table and highlighting within the distribution plot.  For more information on this demo app please check out [my blog post](https://cognitivedataworks.com/blog).

## Installation

The included yml file should allow you to get started quite easily if you have anaconda installed.  

```
conda env create -f environment.yml
```

Alternatively you should be able to create an environment off the rdkit channel and then install the remaining dependencies.  


## Execution

Download or clone the repository locally, make sure your new environment is active and navigate to the directory containing _chemdash.py_ and run:

```
python chemdash.py
```

### 1000 Compounds
Takes about 30 seconds to load the application

### 10000 Compounds 
Takes ~20 minutes to load the application.  If you refer to the blog post above you will see that this is due to the Molecular Grid Map, if you alter the code slightly you should be able to turn that feature off and 10k compounds will likely load in a few minutes.  


