# PISA-lite


## Description

PISA (Protein Interfaces, surfaces and assemblies) v2 is the latest version of the interactive tool PISA used in CCP4. The program was developed for the investigation of macromolecular interfaces (such as proteins, DNA/RNA and ligands) and for the identification of probable quaternary structures or assemblies. 

## Highlights

1. Exploration of macromolecular interfaces 
2. Identification of probable assemblies
3. Database searches of structurally similar interfaces and assemblies

## Compiling and running Pisa v2

```
git clone https://github.com/PDBe-KB/pisa-lite

cd pisa-lite
``` 

1. Modify source directory in the file compile.ssh  
2. Create 'build' directory in the source directory 
3. For compilation run ./compile.ssh 
4. To run the code use executable ./build/pisa . 

##Usage 

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
Options:

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
```
```
pisa name -pdb spec serial_no [cfg] > output_file
```
```
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

```

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

_id_           :         Interface ID

_int_area_     :         Area of interface  (A^2)

_Int_solv_en_  :         Solvation energy (kcal/mol)

_pvalue_       :         Probability that solvation energy gain for interface atom may be greater than binding energy  

_stab_en_      :         Stabilisation energy (Kcal/mol)


For each bond type:

###### Type:

_h-bonds_        :     Hydrogen bonds

_salt-bridges_     :     Salt bridges 

_ss-bonda_         :     Disulfide bonds

_cov-bonds_        :     Covalent Bonds

_other-bonds_      :     Other interface contacts within distance cutoff 4.0 A 

Details (for each bond of each type):

_n_bonds_         :     number of bonds 

_res-1_            :     name of residue 1

_chain-1_          :     name chain 1 (auth) 

_label_asym_id-1_  :     name chain (asym) 1

_orig_label_asym_id-1_     : Asym chain identifier in the original model file

_pdbx_sifts_xref_db_acc-1_ : Uniprot accession number

_pdbx_sifts_xref_db_num-1_ : Sequence position of the UniProt entry that corresponds to the residue mapping

_seqnum-1_                 : Sequence number atom 1 (auth)

_label_seqnum-1_           : Sequence number atom 1 (label)

_atname-1_                : atom 1 name 

_Inscode-1_                : insertion code of 1st linked atom

_res-2_                    : name of Residue 2

_chain-2_                  : name chain 2 (auth)

_label_asym_id-2_          : name chain (asym)

_orig_label_asym_id-2_     : Asym chain identifier in the original model file

_pdbx_sifts_xref_db_acc-2_ : Uniprot accession number

_pdbx_sifts_xref_db_num-2_ : Sequence position of the UniProt entry that corresponds to the residue mapping

_seqnum-2_                 : Sequence number atom 1

_label_seqnum-2_           : Sequence number atom 2 (label)

_atname-2_                 : atom 1 name 

_Inscode-2_                : Insertion code of 2nd linked atom

_Dist_                     : Bond distance (A)

###### Residue list (per interface ID)

For each Residue :

_ser_no_          :   residue numbering 

_name_            :   residue name 

_seq_num_         :   sequence number 

_label_seq_num_   :   sequence number (label)

_ins_code_        :   insertion code residue 

_bonds_           :     

_solv_en_         :   Solvation energy effect (kcal/mol)

_asa_             :   accessible surface area (A^2)

_bsa_             :   Buried surface area  (A^2)

###### Assemblies:

The xml output file for assemblies has the following information:

_diss_energy_     :  maximal free energy of dissociation  (kcal/mol)

_diss_energy0_    :  ground-level free energy of dissociation (kcal/mol)

_asa_             :  accessible surface area  (A^2)

_bsa_             :  buried surface area (A^2)

_entropy_         :  entropy change at dissociation 

_Entropy_0_       :  ground-level entropy change at dissociation 

_diss_area_       :  dissociation Interface Area (A^2)

_int_energy_      :  solvation Energy Gain (kcal/mol)

_n_diss_          :  number of dissociating parts

_Formula_         :  assembly formula

_composition_     :  assembly composition



