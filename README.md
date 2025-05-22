# PISA (Protein Interfaces, Surfaces and Assemblies)

## Description
Note: This version of PISA is identical to that used in the interactive tool [PISA](https://www.ccp4.ac.uk/html/pisa.html) released in the latest version of the CCP4 Software Suite, version 9.0. 

PISA was developed for the investigation of macromolecular interfaces (e.g. between proteins and nucleic acids) and to predict the most probable quaternary structures or assemblies.

In the last version of PISA (CCP4 v 9.0) there is now an option to not calculate and predict assemblies, but only perform assembly interface analysis. [PDBe](https://pdbe.org) uses this option in PISA to calculate macromolecular interaction data which is integrated with functional annotations from [PDBe-KB](https://pdbe-kb.org) and displayed on PDBe and PDBe-KB pages.

## Compiling and running PISA-lite

1.) Download the source code
```
git clone https://github.com/PDBe-KB/pisa-lite
cd pisa-lite
``` 

2.) **Either** export the environment variable `$SRCDIR` to point the compiler to the correct directory. E.g.:

```
export SRCDIR=/path/to/your/pisa-lite
```

**Or** edit the `compile.sh` file in the root directory of the repository (`.../pisa-lite/compile.sh`) directly

Add the complete path to the pisa-lite directory:
```
srcdir=/path/to/your/pisa-lite
```

4.) Compile PISA-lite
```
chmod +x compile.sh
./compile.sh
```

## Configuration

Running PISA-lite requires a configuration file. An example file is [available here](https://github.com/PDBe-KB/pisa-lite/blob/main/setup/pisa_cfg_tmp).

The paths for the following items have to be configured:

* DATA_ROOT
* SRS_DIR
* MOLREF_DIR
* PISTORE_DIR

## Usage 

After compiling, the executable for PISA will be `.../pisa-lite/build/pisa`. 

To run PISA analysis for an input file (PDB or mmCIF file), the following command: 

```
pisa [name] -analyse [coorfile] [options] [cfg]
```

[name]: (required) Session name. The results will be stored in a directory named after the session name. 

[coorfile] : (required) Path to coordinates file (PDB or mmCIF)

[options] : Optional processing keys (see below)

[cfg]: (required) Path to configuration file (it must be provided unless it is already set as environmental variable in PISA_CONF_FILE). Note that this always has to be the last parameter. A default file (`pisa_cfg`) is provided in the repository root. 

###  Options

Type `./build/pisa` to view all the options and their explanations.

## XML file

PISA-lite generates the results in binary format. The data has to be converted into human-readable XML format using the following command:

```
pisa name -xml {interfaces|assemblies} [cfg] > outputfile.xml
```

Note that the `name` parameter has to match an existing session name. For example, if PISA-lite was used to analyse a PDB structure and the name of that session was `analyise_3bow`, then to convert the binary to XML, the name parameter has to be `analyise_3bow` as well.

Also note that `interfaces` and `assemblies` are the exact terms expected by the command.

### Interfaces:

The xml output file for interfaces has the following information:

#### For each interface:

* id           :         Interface ID
* int_area     :         Area of interface  (A^2)
* Int_solv_en  :         Solvation energy (kcal/mol)
* pvalue       :         Probability that solvation energy gain for interface atom may be greater than binding energy  
* stab_en      :         Stabilisation energy (Kcal/mol)

#### For each bond type:

* h-bonds        :     Hydrogen bonds
* salt-bridges     :     Salt bridges 
* ss-bonda         :     Disulfide bonds
* cov-bonds        :     Covalent Bonds
* other-bonds      :     Other interface contacts within distance cutoff 4.0 A 

##### Details (for each bond of each type):

* n_bonds         :     number of bonds 
* res-1            :     name of residue 1
* chain-1          :     name chain 1 (auth) 
* label_asym_id-1  :     name chain (asym) 1
* orig_label_asym_id-1     : Asym chain identifier in the original model file
* pdbx_sifts_xref_db_acc-1 : Uniprot accession number
* pdbx_sifts_xref_db_num-1 : Sequence position of the UniProt entry that corresponds to the residue mapping
* seqnum-1                 : Sequence number atom 1 (auth)
* label_seqnum-1           : Sequence number atom 1 (label)
* atname-1                : atom 1 name 
* Inscode-1                : insertion code of 1st linked atom
* res-2                    : name of Residue 2
* chain-2                  : name chain 2 (auth)
* label_asym_id-2          : name chain (asym)
* orig_label_asym_id-2     : Asym chain identifier in the original model file
* pdbx_sifts_xref_db_acc-2 : Uniprot accession number
* pdbx_sifts_xref_db_num-2 : Sequence position of the UniProt entry that corresponds to the residue mapping
* seqnum-2                 : Sequence number atom 1
* label_seqnum-2           : Sequence number atom 2 (label)
* atname-2                 : atom 1 name 
* Inscode-2                : Insertion code of 2nd linked atom
* Dist                     : Bond distance (A)

#### Residue list (per interface ID)

##### For each Residue :

* ser_no          :   residue numbering 
* name            :   residue name 
* seq_num         :   sequence number 
* label_seq_num   :   sequence number (label)
* ins_code        :   insertion code residue 
* bonds           :     
* solv_en         :   Solvation energy effect (kcal/mol)
* asa             :   accessible surface area (A^2)
* bsa             :   Buried surface area  (A^2)

### Assemblies:

The xml output file for assemblies has the following information:

* diss_energy     :  maximal free energy of dissociation  (kcal/mol)
* diss_energy0    :  ground-level free energy of dissociation (kcal/mol)
* asa             :  accessible surface area  (A^2)
* bsa             :  buried surface area (A^2)
* entropy         :  entropy change at dissociation 
* Entropy_0       :  ground-level entropy change at dissociation 
* diss_area       :  dissociation Interface Area (A^2)
* int_energy      :  solvation Energy Gain (kcal/mol)
* n_diss          :  number of dissociating parts
* Formula         :  assembly formula
* composition     :  assembly composition

## Run in Docker

A Docker image is available for PISA-lite. To run PISA-lite in Docker, use the following command:

``` bash

# Run the analysis
docker run -v $PWD:/data pdbegroup/pisa-lite pisa test-session -analyse /data/example_data/6gve.cif /data/pisa_cfg

# Get the interfaces
docker run -v $PWD:/data pdbegroup/pisa-lite pisa test-session -xml interface /data/pisa_cfg

```

The above command will create a session named `test-session`, run PISA on the provided example data (6gve.cif) and store the results in the session directory. We mount the current directory to the `/data` directory in the Docker container. So the sessions would be stored in the `sessions` directory in the current directory.

## Acknowledgements

This code is based on CCP4 PISA and has been developed through a collaboration between EMBL-EBI and CCP4. 
