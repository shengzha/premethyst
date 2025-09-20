import h5py
import os
import numpy as np
import argparse
import sys


if len(sys.argv) < 3:
    print("\npremethyst cgmaps2h5 [CGmap.gz folder] [output h5 file prefix]\n")
    sys.exit(1)

# to run, type in the command line:
# python /home/groups/ravnica/src/sciMET/sciMET_cellCalls2h5.py [input folder path] [output prefix]

# Define and parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_folder")
parser.add_argument("output_prefix")
args = parser.parse_args()

# Input folder and output file paths
folder_path = args.input_folder
hdf5_path = args.output_prefix + ".h5"

# defaults from facet or premethyst
dataset = 1
compression = "gzip"
compression_opts = 9 # 6 from facet
amethyst_version = "amethyst2.0.0"
overwrite = True


def read_cgmap_gz(file_path):
    try:
        cov = np.genfromtxt(file_path, delimiter='\t', dtype=[('chr', 'S10'), ('nuc', 'S1'), ('pos', int), ('context', 'S3'), ('context2', 'S3'), ('pct', float),('c',int),('t',int)])
        if cov.size == 0:
            print(f"Failed to read data from {file_path}.")
            return None

        cg_idx = np.char.strip(cov['context2']) == b'CG'
        columns_to_select = ['chr','pos', 't', 'c']
        cg = np.sort(cov[cg_idx][columns_to_select], order=['chr', 'pos'])
        
        if all(cg_idx):
            return zip(['CG'], [cg])
        else:
            ch = np.sort(cov[~cg_idx][columns_to_select], order=['chr', 'pos'])
            return zip(['CG', 'CH'], [cg, ch])
    except:
        print(f"Failed to read data from {file_path}.")
        return None


# Ingest source files to HDF5 file
mode = "w" if overwrite else "a"
with h5py.File(hdf5_path, mode) as f:
    # Specify file format version
    metadata_group = f.create_group("metadata")
    metadata_group.create_dataset("version", data=amethyst_version)


    for file_name in os.listdir(folder_path):
        if not file_name.endswith(".CGmap.gz"):
            continue

        file_path = os.path.join(folder_path, file_name)
        cell_id = file_name.replace('.CGmap.gz', '')
        print(f"Reading {cell_id}...")

        context_data = read_cgmap_gz(file_path)
        if context_data is None: # empty
            continue

        for (context, data) in context_data:
            target = f"/{context}/{cell_id}/{dataset}"
            f.create_dataset(target, data=data, compression=compression, compression_opts=compression_opts)

print("Files combined and written to HDF5.")