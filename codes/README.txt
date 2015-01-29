-- INTRODUCTION --

The Python scripts in this folder processes raw data into intermediate products that are used by a separate script to produce the figures for the paper. This separate script can be found in ../images/plots.py. The raw data can be downloaded as described below. You may also find the data flow chart (DataFlow.pdf) helpful.

-- PYTHON SET UP --

These scripts have been tested with Python 3.3.2 on CentOS 6.6. The modules (and their versions) are listed in requirements.txt. You can install the exact versions into your Python environment by running:

pip install -r requirements.txt

It is highly recommended that a virtual environment is set up for this. For more information, see http://docs.python-guide.org/en/latest/dev/virtualenvs/.

-- USING THE SCRIPTS --

Some paramters in the scripts, e.g. the directories of the data, may need to be modified. The entire family of scripts provided here, barring changes to external factors such as changes to the URLs from which the raw data are downloaded, should be able to, in principle, fully reproduce the results.

As many of these scripts are written for broader and more varied purposes, they may produce data that are redundant to this paper. For example, in the regime assignment scripts, you will find references to "sub-CD regimes"; these were directions in the research that have been discarded and hence they can be ignored. For the same reason, the scripts may execute their tasks in a roundabout way, but they should work just fine.

See DataFlow.pdf on how the scripts relate to different raw and intermediate datasets.

-- DATA FLOW --

The data flow chart (DataFlow.pdf) plots the relationship between the raw data, the intermediate products and the final diagrams; and how the scripts connect these different components. Note that this chart indicates the flow of data. It allows you, for example, to find out which raw data sources does each diagram use.

However, the scripts that produces them and the intermediate data may need more than just those indicated in the flow chart. For example, the reconstruct_longer.py script that produces the second intermediate data for ED Fig. 1 needs GPCP 2.2 to execute, even though it is not stated in the flow chart for ED Fig. 1.

-- RAW DATA SOURCES --

The cloud regimes can be obtained in two ways. One, the daily cloud regimes for the extended tropics can be downloaded directly at http://mwac.its.monash.edu/mwac/pub/listPubCollections.jspx. Two, the scripts regime_assign.py, regime_netcdf.py, and regime_assign_3hr.py allow you to produce the regimes from the raw ISCCP D1 data (you will also need the centroids files in the data folder). The first method is easy; however, you will not be able to get the three-hourly cloud regimes that is used to produce daytime-averaged rainfall values from TRMM 3B42. (Note that this three-hourly cloud regimes is NOT the IR-only cloud regimes available for download.)

TRMM 3B42 (TMPA): use download_trmm3b42.py (takes days!)
TRMM 3A25: use download_trmm3a25.py, then convert to HDF5 using h4toh5*
GPCP 2.2: download directly from NOAA PSD at http://www.esrl.noaa.gov/psd/data/gridded/data.gpcp.html
GPCP 1DD: use download_gpcp1dd.py
ERA-Interim ω500: download directly at http://apps.ecmwf.int/datasets/ (select all four times of the day)
ISCCP D1: use download_isccp.py

*: 3A25 is recorded in HDF4. In principle, the module netCDF4 should be able to open such files, provided the software netcdf-4 is compiled with the appropriate flags. However, as my installation is NOT configured this manner, I cannot read HDF4. I resolved this by converting to HDF5 first. This tool can be found here: http://www.hdfgroup.org/h4toh5/.

-- CONTACT --

If you need help on using the scripts, please feel free to contact me.

Jackson Tan (jackson.tan@nasa.gov)
27 January 2015
