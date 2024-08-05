FROM ubuntu:22.04
# NOTE: ubuntu 24.04 uses python 3.12 which is not working w/ some packages, so we stick w. 22.04 and python 3.10.

RUN mkdir -p /opt/chemdash
WORKDIR /opt/chemdash

# Do the slow parts first, so rebuilds are fast when we make changes to the python tree...

# Non-python deps.
COPY install-deps-ubuntu-22.04.sh /opt/chemdash/install-deps-ubuntu-22.04.sh
RUN ./install-deps-ubuntu-22.04.sh

# Python deps.
COPY ./requirements.txt /opt/chemdash/requirements.txt
ENV PIP_ROOT_USER_ACTION=ignore
ENV PIP_NO_CACHE_DIR=off

# These two have a dependency resolution order problem.  Handle manually though they are also in requirements.txt
RUN python3 -m pip install -r requirements.txt
RUN python3 -m pip install -e git+https://github.com/gatagat/lap.git@v0.4.0#egg=lap

# Application code.  This updates most frequently but is a cheap rebuild. 
COPY . /opt/chemdash

# In production we override with gunicorn, etc. 
ENTRYPOINT ["/opt/chemdash/chemdash.py"]
CMD ["--help"]

