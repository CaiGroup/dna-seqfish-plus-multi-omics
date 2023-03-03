## 3D visualization 
### 3D visualization is performed using mayavi package under python
To install mayavi 4.8.1 package, please follow the document:
https://docs.enthought.com/mayavi/mayavi/installation.html#installing-with-conda-forge
Conda-forge installation is recommended, python 3.8 was used.

To run H4K20me3_MajSat_volume_3D_plot.py file:
Organize all IF intensity file, segmentation file and detected DNA bin position information in example_data folder
change config.ini file to corresponding cell need to be visualized 

then run code under mayavi environment:
mayavi2 -x H4K20me3_MajSat_volume_3D_plot.py


