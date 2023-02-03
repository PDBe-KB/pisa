# PISA-lite

## Description

PISA-lite (Protein Interfaces, Surfaces and Assemblies) is a lightweight version of the interactive tool [CCP4 PISA](https://www.ccp4.ac.uk/html/pisa.html). CCP4 PISA was developed to investigate macromolecular interfaces (e.g. between proteins and nucleic acids) and to identify the most probable quaternary structures or assemblies. 

PISA-lite differs from CCP4 PISA in that it does not calculate assemblies, but performs interface analysis significantly faster. [PDBe](https://pdbe.org) uses PISA-lite to calculate macromolecular interaction data which is integrated with functional annotations from [PDBe-KB](https://pdbe-kb.org) and displayed on PDBe and PDBe-KB pages.

## Highlights

1. Exploration of macromolecular interfaces 
2. Identification of probable assemblies
3. Database searches of structurally similar interfaces and assemblies

## Compiling and running PISA-lite

1.) Download the source code
```
git clone https://github.com/PDBe-KB/pisa-lite

cd pisa-lite
``` 

2.) Edit the `compile.sh` file in the root directory of the repository (`.../pisa-lite/compile.sh`)

Add the complete path to the pisa-lite directory:
```
srcdir=/complete/path/to/pisa-lite
```

4.) Compile PISA-lite
```
chmod +x compile.sh
./compile.sh
```

## Configuration

Running PISA-lite

## Usage 

After compiling, the executable for PISA will be `.../pisa-lite/build/pisa`. 

To run PISA analysis for an input file (PDB or mmCIF file), the following command: 

```
pisa name -analyse coorfile [options] [cfg]
```
```
name: Session name (required). The results will be stored in a directory named after the session name. 

coorfile : coordinates file (PDB or mmCIF)

[cfg]: Configuration file (it must be provided unless it is already set as environmental variable in PISA_CONF_FILE ). See description bellow. 

[options] : Optional processing keys (see below)

```
######  Options:

```
--lig-exclude = 'list'
```

Entries:

'list' (default). Exclude ligands in a list when calculating the oligomeric state. The list must be coma separated, no spaces, case-sensitive list of ligand names.
'all'       Exclude all ligands
'agents'    Exclude ligands listed in agents.dat ( path pointing to agents.dat specified in configuration file, see below )

```
--lig=auto = fixed
```
Auto processing for ligands (Fixed default). Ligands may be forcefully fixed on protein surface subject to task complexity and data available.

Entries:

fixed:  fix ligands on protein surface
free:   ligands are not fixed on the protein surface
```
pisa name -list {interfaces|monomers|assemblies|stock} [cfg]  
```
Outputs lists of interfaces, monomers, assemblies and assembly stock
```
pisa name -view spec serial_no [cfg]   
```
View an interface, monomer, assembly or a dissociate in Rasmol 

spec={interface|monomer|assembly|dissociate}
serial_no   Serial number shown in the corresponding list.
```
pisa name -mg spec serial_no [cfg]
 ```          
View an interface, monomer, assembly or a dissociate in CCP4-MG
```
pisa name -download spec serial_no [cfg] > output_file

pisa name -pdb spec serial_no [cfg] > output_file

pisa name -cif spec serial_no [cfg] > output_file
```
Download an interface, monomer, assembly or a dissociate. Flags -download and -pdb are equivalent and output PDB format. Flag -cif will output in mmCIF format.

```
pisa name -detail spec serial_no [cfg]
```
Get details of an interface, monomer or assembly
```
pisa name -350n assembly_serial_no [cfg]
```
Generate Remark 350
Special use:
assembly_serial_no=0: generate remark for complexe made by the content of asu (i.e. the bare content of input file), without symmetry analysis
assembly_serial_no<0: generate remarks 300 and 350 for the whole crystal split number -assembly_serial_no. E.g., use assembly_serial_no=-1.

```
pisa name -350n -1 [cfg]
```
Remark 350 for automatic annotation

```
pisa name -xml {interfaces|assemblies} [cfg] > outputfile.xml
```
Generates XML output

```
pisa name -erase [cfg]
```
Erases session data

```
pisa -cfg-template [cfg]
```
A configuration file template is generated using this command.  


## Dependencies 

For running PISA with command line (offline/locally), the following keys must be specified in the configuration file:

```
DATA_ROOT    Path pointing to the main directory where the session is 
SRS_DIR           Path pointing to directory that contains SRS files 
MOLREF_DIR   Path pointing to directory that contains Molref files 
PISTORE_DIR  Path pointing to store directory containing the files: 
agents.dat
asm_params.dat
intfstats.dat
syminfo_pisa.lib
rcsb_symops.dat
HELP_DIR         Path pointing to directory containing Html-formatted help files 
```

## XML file

```
pisa name -xml {interfaces|assemblies} [cfg] > outputfile.xml

```
###### Interfaces:

The xml output file for interfaces has the following information:

For each interface:

```
id           :         Interface ID

int_area     :         Area of interface  (A^2)

Int_solv_en  :         Solvation energy (kcal/mol)

pvalue       :         Probability that solvation energy gain for interface atom may be greater than binding energy  

stab_en      :         Stabilisation energy (Kcal/mol)
```

For each bond type:

###### Type:

```
h-bonds        :     Hydrogen bonds

salt-bridges     :     Salt bridges 

ss-bonda         :     Disulfide bonds

cov-bonds        :     Covalent Bonds

other-bonds      :     Other interface contacts within distance cutoff 4.0 A 

```

Details (for each bond of each type):

```
n_bonds         :     number of bonds 

res-1            :     name of residue 1

chain-1          :     name chain 1 (auth) 

label_asym_id-1  :     name chain (asym) 1

orig_label_asym_id-1     : Asym chain identifier in the original model file

pdbx_sifts_xref_db_acc-1 : Uniprot accession number

pdbx_sifts_xref_db_num-1 : Sequence position of the UniProt entry that corresponds to the residue mapping

seqnum-1                 : Sequence number atom 1 (auth)

label_seqnum-1           : Sequence number atom 1 (label)

atname-1                : atom 1 name 

Inscode-1                : insertion code of 1st linked atom

res-2                    : name of Residue 2

chain-2                  : name chain 2 (auth)

label_asym_id-2          : name chain (asym)

orig_label_asym_id-2     : Asym chain identifier in the original model file

pdbx_sifts_xref_db_acc-2 : Uniprot accession number

pdbx_sifts_xref_db_num-2 : Sequence position of the UniProt entry that corresponds to the residue mapping

seqnum-2                 : Sequence number atom 1

label_seqnum-2           : Sequence number atom 2 (label)

atname-2                 : atom 1 name 

Inscode-2                : Insertion code of 2nd linked atom

Dist                     : Bond distance (A)
```
###### Residue list (per interface ID)

For each Residue :

```
_ser_no_          :   residue numbering 

_name_            :   residue name 

_seq_num_         :   sequence number 

_label_seq_num_   :   sequence number (label)

_ins_code_        :   insertion code residue 

_bonds_           :     

_solv_en_         :   Solvation energy effect (kcal/mol)

asa             :   accessible surface area (A^2)

bsa             :   Buried surface area  (A^2)
```
###### Assemblies:

The xml output file for assemblies has the following information:

```
diss_energy     :  maximal free energy of dissociation  (kcal/mol)

diss_energy0    :  ground-level free energy of dissociation (kcal/mol)

asa             :  accessible surface area  (A^2)

bsa             :  buried surface area (A^2)

entropy         :  entropy change at dissociation 

Entropy_0       :  ground-level entropy change at dissociation 

diss_area       :  dissociation Interface Area (A^2)

int_energy      :  solvation Energy Gain (kcal/mol)

n_diss          :  number of dissociating parts

Formula         :  assembly formula

composition     :  assembly composition



