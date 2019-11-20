FROM continuumio/miniconda3
RUN conda update conda
RUN conda install -c conda-forge jupyterlab
RUN conda install -c conda-forge scikit-learn
RUN conda install -c conda-forge statsmodels
RUN conda install -c conda-forge matplotlib
RUN conda install -c conda-forge scipy
RUN conda install -c conda-forge cmocean
RUN conda install -c conda-forge altair vega_datasets
RUN conda install -c conda-forge geopandas
RUN conda install -c conda-forge xarray
RUN conda install -c conda-forge cartopy
RUN conda install -c conda-forge nodejs
RUN conda install -c conda-forge descartes
RUN conda install -c conda-forge dask
RUN conda install -c conda-forge gsw
RUN pip install regionmask
RUN jupyter labextension install jupyterlab_vim 
#RUN curl -sL https://deb.nodesource.com/setup_12.x | bash - \
#	&& apt-get install -y nodejs
#RUN conda install -c conda-forge netcdf4
CMD jupyter lab --notebook-dir=/opt/skipjack-bigeye-separation --ip='0.0.0.0' --port=8877 --no-browser --allow-root
