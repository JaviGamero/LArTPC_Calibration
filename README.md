# LAr-TPC calibration via light

## Configuring a python virtual environment
With both ROOT and python (a working version with root) correctly installed, you need to install the python package `venv` to create the virtual environment: 
```bash
python3.11 -m pip install venv <env_name>
```

Afterwards, you have to activate it and add root settings for python: 
```bash
source <env_name>/bin/activate
source <path_to_installed_root_main_folder>/bin/thisroot.sh
```

With this done, you will be able to use root inside your environment.  

In this repository you may find a file called ***requirements.txt*** with the python necessary packages. Inside your virtual environment, use the next command to install it: 
```bash
pip install -r requirements.txt
``` 

Also, if you want to use jupyter notebook's technology, with everything installed use: 
```bash
root --notebook
```