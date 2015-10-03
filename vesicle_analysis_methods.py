import MDAnalysis as mda
import time
import os
import numpy
from multiprocessing import Pool
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob

# Set global variables

systems_dir = "/nfs/homes/ikenney/Projects/vesicles/scripts/analysis/systems/"
fig_dir = "/nfs/homes/ikenney/Projects/vesicles/scripts/analysis/figs/"
hostname = "spuddafett"
data_dir = "/nfs/homes/ikenney/Projects/vesicles/scripts/analysis/data/"
critical_components = ["load_tpr_xtc","load_pdb","load_gro","load_pdb_xtc","load_gro_xtc"]

function_dict = {}

# Initial run will identify missing systems and missing data

def integrity(systems,N=40):
    tobuild = [] # Elements will have _build() run on them
    tofill = [] # Elements will have _fill() run on them
    
    # Detect new systems and missing data
    for x in [y.split("/")[-1] for y in glob.glob(systems_dir+"*")]:
        if x not in systems.index.tolist(): # Missing systems 
            print "\'{}\' not found in data set. Build? (y/n)".format(x)
            userinput = raw_input("Build:\t")
            if userinput.lower() == 'y':
                tobuild.append(x)
                print("Added to build list.\n")
        else: # Missing data
            for y in systems.columns.tolist():
                if systems[y].isnull().loc[x] and y in critical_components:
                    print("{} missing from {}. Fill in data? (y/n)".format(y,x))
                    userinput = raw_input("Fill:\t")
                    if userinput.lower() == 'y':
                        tofill.append([y,x])
                        print("Added to fill list.\n")
    newsys = systems
    
    for s in tobuild:
        newsys = newsys.append(_build(s,systems,N=N))
    for s in tofill:
        newsys = newsys.append(_fill(s,systems,N=N))
    pd.to_pickle("backup.df")
    return newsys.sort(columns=['sizes'])

### POPULATION FUNCTIONS ###

def populate_rad_gyr(sysname,table): # Handle radius of gyration calculations
    table['rgyr'].loc[sysname] = _radius_gyration_performance([table['tops'],table['traj']])
    table['rgyr_mean'].loc[sysname] = np.mean(table['rgyr'].loc[sysname])
    table['rgyr_std'].loc[sysname] = np.std(table['rgyr'].loc[sysname])
    return True

def populate_tpr_xtc(sysname,table,N=40): # Handle loading tpr and xtc
    column = "load_tpr_xtc"
    table[column].loc[sysname] = _tpr_xtc_load(sysname,table,N)
    table[column+'_mean'].loc[sysname] = np.mean(table[column].loc[sysname])
    table[column+'_std'].loc[sysname] = np.std(table[column].loc[sysname])
    return True

def populate_gro(sysname,table,N=40): # Handle loading gro
    column = "load_gro"
    table[column].loc[sysname] = _gro_load(sysname,table,N)
    table[column+'_mean'].loc[sysname] = np.mean(table[column].loc[sysname])
    table[column+'_std'].loc[sysname] = np.std(table[column].loc[sysname])
    return True

def populate_pdb(sysname,table,N=40): # Handle loading pdb 
    column = "load_pdb"
    table[column].loc[sysname] = _pdb_load(sysname,table,N)
    table[column+'_mean'].loc[sysname] = np.mean(table[column].loc[sysname])
    table[column+'_std'].loc[sysname] = np.std(table[column].loc[sysname])
    return True

def populate_pdb_xtc(sysname,table,N=40): # Handle loading pdb xtc
    column = "load_pdb_xtc"
    table[column].loc[sysname] = _pdb_xtc_load(sysname,table,N)
    table[column+'_mean'].loc[sysname] = np.mean(table[column].loc[sysname])
    table[column+'_std'].loc[sysname] = np.std(table[column].loc[sysname])
    return True

def populate_gro_xtc(sysname,table,N=40): # Handle loading gro and xtc
    column = "load_gro_xtc"
    table[column].loc[sysname] = _gro_xtc_load(sysname,table,N)
    table[column+'_mean'].loc[sysname] = np.mean(table[column].loc[sysname])
    table[column+'_std'].loc[sysname] = np.std(table[column].loc[sysname])
    return True

### Computation Functions ###

function_dict['rgyr'] = _radius_gyration_performance
def _radius_gyration_performance(a):
    u = mda.Universe(a[0][-1],a[1][-1])
    vals = []
    for frame in u.trajectory:
        start = time.time()
        u.atoms.radius_of_gyration()
        vals.append(time.time()-start)
    return vals

function_dict['load_tpr_xtc'] = _tpr_xtc_load
def _tpr_xtc_load(sysname,table,N=40):
    def load():
        start = time.time()
        mda.Universe(table['tops'].loc[sysname],table['traj'].loc[sysname])
        return time.time() - start
    return [load() for _ in range(N)]

function_dict['load_gro'] = _gro_load
def _gro_load(sysname,table,N=40):
    def load():
        start = time.time()
        mda.Universe(table['gros'].loc[sysname])
        return time.time() - start
    return [load() for _ in range(N)]

function_dict['load_pdb'] = _pdb_load
def _pdb_load(sysname,table,N=40):
    def load():
        start = time.time()
        mda.Universe(table['pdbs'].loc[sysname])
        return time.time() - start
    return [load() for _ in range(N)]

function_dict['load_gro_xtc'] = _gro_xtc_load
def _gro_xtc_load(sysname,table,N=40):
    def load():
        start = time.time()
        mda.Universe(table['gros'].loc[sysname],table['traj'].loc[sysname])
        return time.time() - start
    return [load() for _ in range(N)]

function_dict['load_pdb_xtc'] = _pdb_xtc_load
def _pdb_xtc_load(sysname,table,N=40):
    def load():
        start = time.time()
        mda.Universe(table['pdbs'].loc[sysname],table['traj'].loc[sysname])
        return time.time() - start
    return [load() for _ in range(N)]

def _build(system_name,systems,N=40):
    base_dir = systems_dir + system_name
    print("Building system in: " + base_dir)
    buildrow = pd.DataFrame(index=[system_name],columns=systems.columns.tolist())
    if os.path.exists(base_dir+"/nvt/nvt.tpr"):
        buildrow['tops'] = base_dir+"/nvt/nvt.tpr"
    else:
        print("Missing path:\t" + base_dir+"/nvt/nvt.tpr")
        
    if os.path.exists(base_dir+"/nvt/analysis.xtc"):
        buildrow['traj'] = base_dir+"/nvt/analysis.xtc"
    else:
        print("Missing path:\t" + base_dir+"/nvt/analysis.xtc")   
        
    if os.path.exists(base_dir+"/emin/emin.gro"):
        buildrow['gros'] = base_dir+"/emin/emin.gro"
    else:
        print("Missing path:\t" + base_dir+"/emin/emin.gro") 
        
    if os.path.exists(base_dir+"/emin/emin.pdb"):
        buildrow['pdbs'] = base_dir+"/emin/emin.pdb"
    else:
        print("Missing path:\t" + base_dir+"/emin/emin.pdb") 
        
    if buildrow['gros'].notnull().loc[system_name]:
        populate_gro(system_name,buildrow,N=N)
        
    if buildrow['pdbs'].notnull().loc[system_name]:
        populate_pdb(system_name,buildrow,N=N)
            
    if buildrow['tops'].notnull().loc[system_name] and buildrow['traj'].notnull().loc[system_name]:
        populate_tpr_xtc(system_name,buildrow,N=N)
            
    if buildrow['gros'].notnull().loc[system_name] and buildrow['traj'].notnull().loc[system_name]:
        populate_gro_xtc(system_name,buildrow,N=N)
            
    if buildrow['pdbs'].notnull().loc[system_name] and buildrow['traj'].notnull().loc[system_name]:
        populate_pdb_xtc(system_name,buildrow,N=N)
            
    buildrow['sizes'] = mda.Universe(buildrow['tops'].ix[-1]).atoms.n_atoms
    return buildrow

def _fill(system_info,system,N=40):
    system_name = system_info[1]
    column = system_info[0]
    
    base_dir = systems_dir + system_name

    print("Filling in data from: " + base_dir)

    comps = system.split("_")[1:-1]
    continuation = True
    
    for c in comps:
        if system[column].isnull().loc[system_name]:
            print("Missing components to fill data table: {}, {}".format(system_name,column))
            continuation = False

    if not continuation:
        exit(1)

    function_dict[column](system_name,system,N=N)
        
    return True
