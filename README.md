# ParaView Plugins
This repository contains some sample plugins for ParaView, written in Python. Plugins can be readers, writers or filters. The following plugins are currently available:
1. [Incompact3d](https://www.incompact3d.com/) Reader
2. Incompact3d Statistics Data Reader
3. Analyse Incompact3d filter
4. Compressible Flow filter


To add a python plugin to ParaView, follow the simple steps below after downloading the python file (e.g. CompressibleFlow.py):
1. Open Tools > Manage Plugins...
2. Select Load New and open the downloaded python file 
3. The plugin is loaded and you can enable Auto Load under the plugin name if you want it available for all your future sessions.
4. Once the plugin is loaded:
   1. If the plugin is a reader, it will show up in file open options (along with other file formats in the list). 
   2. In case of a filter, you will see the filer in the list of Filters (or Ctrl+Space to search for it).
   3. Writers will appear as options in Save Data 

Detailed information about how to use the Incompact3d reader is documented in the repo's [Wiki](https://github.com/rajesh-ae/ParaView_plugins/wiki) . 
