#!/usr/bin/env python
# coding: utf-8

# Imports
import os
import synapseclient
import synapseutils

# Use "os" library to read your api keys, or you can hardcode it if that's easier
AUTHID = os.environ.get('SYNAPSE_USERNAME')
AUTHPW = os.environ.get('SYNAPSE_PASSWORD')

# Log-in
syn = synapseclient.Synapse() 
syn.login(AUTHID, AUTHPW)

# Create some constants to store the paths to the data
# This is the path to the directory that contains the data in the structure you want to upload
DIRECTORY_FOR_MY_PROJECT = "/path/to/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq"

# You will create a manifest file that will will contain the above structure
PATH_TO_MANIFEST_FILE = "manifest-for-upload.tsv"

# Step 1: Find the synapse ID of the project
my_project_id = "syn63660123"

# Step 2: Create a manifest TSV file to upload data in bulk
# Note: When this command is run it will re-create your directory structure within
# Synapse. Be aware of this before running this command.
# If folders with the exact names already exists in Synapse, those folders will be used.
synapseutils.generate_sync_manifest(
    syn=syn,
    directory_path=DIRECTORY_FOR_MY_PROJECT,
    parent_id=my_project_id,
    manifest_path=PATH_TO_MANIFEST_FILE,
)

# Step 3: After generating the manifest file, we can upload the data in bulk
synapseutils.syncToSynapse(
    syn=syn, manifestFile=PATH_TO_MANIFEST_FILE, sendMessages=False
)

# DONE!
# -----
