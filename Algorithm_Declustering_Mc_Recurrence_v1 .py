#!/usr/bin/env python
# coding: utf-8

# In[1]:


### Shivaprakash, April 1, 2022

###--------------------------------------------------------------###
#  OQ engine HMTK script for catalogue declustering,...
#  Magnitude of completeness, and G-R parameter calculation
###--------------------------------------------------------------###       

#### Import the required libraries

get_ipython().run_line_magic('matplotlib', 'inline')

# Python dependences

import os
import numpy as np   # Numpy - Python's numerical library
import matplotlib.pyplot as plt  # Matplotlib - Python's plotting library
from copy import deepcopy   # Python module for copying objects


# In[2]:


#Declustering algorithm: using the Gardner & Knopoff (1974) Distance & Time Windows

###   Set-up the tools   ###

# Import HMTK I/O Tools
from openquake.hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter


# HMTK Declustering Tools
from openquake.hmtk.seismicity.declusterer.dec_afteran import Afteran
from openquake.hmtk.seismicity.declusterer.dec_gardner_knopoff import GardnerKnopoffType1
from openquake.hmtk.seismicity.declusterer.distance_time_windows import GardnerKnopoffWindow, GruenthalWindow, UhrhammerWindow


# HMTK Completeness Tools
from openquake.hmtk.seismicity.completeness.comp_stepp_1971 import Stepp1971

print('Import Okay')


# In[3]:



###---> Import the catalogue

catalog_sk = 'eq_cat_isc_usgs_skorea_1900_2022_input.csv'
parser = CsvCatalogueParser(catalog_sk)
catalog  = parser.read_file()

print('Catalogue contains %s events' % catalog.get_number_events())


# In[4]:


###---> Set up the declustering algorithm
# Step 1 - set-up the tool
gardner_knopoff = GardnerKnopoffType1()


# Create the configuration file
declust_config = {'time_distance_window': GardnerKnopoffWindow(),
                  'fs_time_prop': 1.0}
print(declust_config)


# In[5]:


# Run the declustering
print('Running declustering ...')
cluster_index, cluster_flag = gardner_knopoff.decluster(catalog , declust_config)


# In[6]:


print('Done!')
print('%s Clusters found' % np.max(cluster_index))
print('%s Non-poissionian events identified' % np.sum(cluster_flag != 0))


# In[7]:


data = np.column_stack([catalog.get_decimal_time(),
                        catalog.data['magnitude'],
                        catalog.data['longitude'],
                        catalog.data['latitude'], cluster_index, cluster_flag])
print('      Time    Magnitude    Long.    Lat.   Cluster No. Index (-1 = foreshock, 0 = mainshock, 1 = afterschock)')
for row in data:
    print('%14.8f  %6.2f  %8.3f  %8.3f  %6.0f  %6.0f' %(row[0], row[1], row[2], row[3], row[4], row[5]))    


# In[8]:


catalog.get_number_events()    


# In[9]:


import copy
from copy import deepcopy

# Copying the catalogue and saving it under a new name "catalogue_dec"(declustered catalogue) 
catalog_dec = deepcopy(catalog )


# In[10]:


# Logical indexing: Chossing the outputs for the main events: Cluster_flag = 0 
mainshock_flag = cluster_flag == 0 


# In[11]:


# Filtering the foreshocks and aftershocks in the copy of the catalogue 
catalog_dec.purge_catalogue(mainshock_flag)


# In[13]:


# Printing the number of events considered main shocks
print('Declustering: Oaky')
print("Number of events in original catalogue: %g" % catalog.get_number_events())
print('Number of mainshocks: %g' % catalog_dec.get_number_events())


# In[14]:


## Saving the catalogue

# Selecting path and name for the output file 

import os

output_catalog_dec = 'output/South_Korea_catalogue_declustered_v1.csv'

if os.path.exists(output_catalog_dec):
    os.remove(output_catalog_dec)

# Call the method and save the output file under the name "cat_csv"
cat_csv = CsvCatalogueWriter(output_catalog_dec) 

# Write the purged catalogue
cat_csv.write_file(catalog_dec)
print("Catalogue successfully written to %s" % output_catalog_dec)


# In[18]:


catalog_dec.get_number_events()


# In[20]:


## Modified by Prakash 

## Calculating the Mc

comp_config = {'magnitude_bin': 0.1, 'time_bin': 5., 'increment_lock': True }
from openquake.hmtk.seismicity.completeness.comp_stepp_1971 import Stepp1971
completeness_algorithm = Stepp1971 ()
completeness_table = completeness_algorithm.completeness (catalog_dec , comp_config )
print(comp_config)
print(completeness_table)
#plot_stepp_1972.create_stepp_plot(comp_config, completeness_table)


# In[22]:


## calculating the recurrence 


mle_config = {'magnitude_interval': 0.1 , 'Average Type': 'Weighted', 
              'reference_magnitude': None }

from openquake.hmtk.seismicity.occurrence.b_maximum_likelihood import BMaxLikelihood
recurrence = BMaxLikelihood ()

bval , sigmab , aval , sigmaa = recurrence.calculate ( catalog_dec , mle_config , 
                                                      completeness = completeness_table )


# In[23]:


bval


# In[24]:


aval


# In[ ]:





# In[ ]:




