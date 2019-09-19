FROM continuumio/miniconda3
RUN conda install jupyterlab
RUN conda update jupyterlab
RUN conda install scikit-learn
RUN conda install statsmodels
RUN conda install netCDF4
RUN conda install matplotlib
RUN conda install scipy
RUN conda install -c conda-forge cmocean
RUN conda install -c conda-forge basemap
CMD jupyter lab --notebook-dir=/opt/skipjack-bigeye-separation --ip='0.0.0.0' --port=8877 --no-browser --allow-root
