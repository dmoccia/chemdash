# chemdash
Chemical Dataset Exploration Tool Built on Plotly Dash


The chemdash app is a simple analysis application built upon the plotly dash framework.  As input the application takes a _csv_ file requiring both **smiles** and **Compound_id** fields, other fields will be retained and displayed in the data table portion of the application.  A UMAP (similar to t-SNE/PCA) plot will be generated as well as a Molecular Grid Map (MGM).   Depending on the column selection in the data table, a distribution plot for the given column will also be displayed.  Selection of data points in either the UMAP or MGM plot will result in filtering of the data table and highlighting within the distribution plot.  For more information on this demo app please check out [my blog post](https://cognitivedataworks.com/2020/02/building-a-gan-compound-explorer-in-plotly-dash).

## Installation & Execution

Download the repository locally, then install the dependencies from the top of the git repo:

### Docker

This works on platforms that support docker:
_Note: this works on Macs with Apple processors, but be sure to keep platform set to linux/amd64 so Intel-based packages work._
```
docker buildx build --platform linux/amd64 -t chemdash:latest .
```

To start chemdash from the CLI:
```
docker run --rm --name chemdash -p 8000:8000 -v /tmp:/tmp -v /var/tmp:/var/tmp -d chemdash:latest dataset_800.csv  --port 8000 --debug 
```

Then open a browser while watching the logs in the terminal:
```
open http://0.0.0.0:8000
docker logs -f chemdash
```

To stop it from the CLI:
```
docker kill chemdash
```

### Non-docker install and run:

A non-docker install can be useful for debugging.

Install non-python dependencies on ubuntu (note 22.04 is required, since some modules in 24.04, with python 3.12, have problems.)
```
./install-deps-ubuntu-22.04.sh
```

Install non-python dependencies on Mac:
```
./install-deps-osx.sh
```

Install the python dependencies in a virtualenv:
``
python3.10 -m venv venv
. ./venv/bin/activate
pip install -r requirements.txt
```

To run, make sure your new environment is active and navigate to the directory containing _chemdash.py_ and run it.
```
. venv/bin/activate
python3 ./chemdash.py dataset_800.csv --port 8000 --debug
```

### PyCharm

The repo contains two PyCharm run configurations, both of which use the venv, each with differen test datasets. 
You make either of those as above, and then specify either `chemdash venv` or `chemdash docker` to run the respective way.

NOTE: PyCharm debugging works best with the `chemdash venv` option.


## Data Loading

### 1000 Compounds
Takes about 30 seconds to load the application

### 10000 Compounds 
Takes ~20 minutes to load the application.  If you refer to the blog post above you will see that this is due to the Molecular Grid Map, if you alter the code slightly you should be able to turn that feature off and 10k compounds will likely load in a few minutes.  


## Usage

When the app is running, open your browser to http://0.0.0.0:8000


