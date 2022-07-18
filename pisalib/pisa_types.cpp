// $Id: pisa_types.cpp $
// =================================================================
//
//    15.02.19   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -------------------------------------------------------------------
//
//  **** Module  :  PISA_Types <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2019
//
// =================================================================
//

#include "pisa_types.h"

namespace pisa  {

  mmdb::cpstr  BuildDate  = "15-03-2019 ";
  mmdb::cpstr  auth_email = "eugene.krissinel@stfc.ac.uk";

  ///   PISA configuration file must be either specified in the command
  /// or provided by environment variable pisa_conf_file. All further
  /// configuration is read from this configuration file. It must be
  /// an absolute path, i.e. starting from '/'.
  mmdb::cpstr  conf_file = "PISA_CONF_FILE";

  mmdb::cpstr  ccp4_env       = "CCP4"                ;
  mmdb::cpstr  conf_path_ccp4 = "/share/pisa/pisa.cfg";

  ///  prefix for session directories
  mmdb::cpstr  default_session_prefix  = "pisa_";

  ///  configuration file tags
  mmdb::cpstr  data_root_tag        = "DATA_ROOT"          ;
  mmdb::cpstr  srs_dir_tag          = "SRS_DIR"            ;
  mmdb::cpstr  molref_dir_tag       = "MOLREF_DIR"         ;
  mmdb::cpstr  pistore_dir_tag      = "PISTORE_DIR"        ;
  mmdb::cpstr  help_dir_tag         = "HELP_DIR"           ;
  mmdb::cpstr  rasmol_com_tag       = "RASMOL_COM"         ;
  mmdb::cpstr  jmol_com_tag         = "JMOL_COM"           ;
  mmdb::cpstr  ccp4mg_com_tag       = "CCP4MG_COM"         ;
  mmdb::cpstr  session_prefix_tag   = "SESSION_PREFIX"     ;
  mmdb::cpstr  php_uri_tag          = "PHP_URI"            ;
  mmdb::cpstr  help_uri_tag         = "HELP_URI"           ;
  mmdb::cpstr  dnl_url_tag          = "DOWNLOAD_URL"       ;
  mmdb::cpstr  jsrview_uri_tag      = "JSRVIEW_URI"        ;
  mmdb::cpstr  pdb_dir_tag          = "PDB_DIR"            ;
  mmdb::cpstr  pdb_dir_format_tag   = "PDB_DIR_FORMAT"     ;
  mmdb::cpstr  erase_time_tag       = "ERASE_TIME"         ;
  mmdb::cpstr  expiry_time_tag      = "EXPIRY_TIME"        ;

  mmdb::cpstr  pdb_plain_pdb_key      = "PDB_PLAIN_PDB"      ;
  mmdb::cpstr  pdb_plain_pdb_gz_key   = "PDB_PLAIN_PDB_GZ"   ;
  mmdb::cpstr  pdb_plain_mmcif_key    = "PDB_PLAIN_mmCIF"    ;
  mmdb::cpstr  pdb_plain_mmcif_gz_key = "PDB_PLAIN_mmCIF_GZ" ;
  mmdb::cpstr  pdb_split_pdb_key      = "PDB_SPLIT_PDB"      ;
  mmdb::cpstr  pdb_split_pdb_gz_key   = "PDB_SPLIT_PDB_GZ"   ;
  mmdb::cpstr  pdb_split_mmcif_key    = "PDB_SPLIT_mmCIF"    ;
  mmdb::cpstr  pdb_split_mmcif_gz_key = "PDB_SPLIT_mmCIF_GZ" ;


  ///  Store files
  mmdb::cpstr  syminfo_file         = "syminfo_pisa.lib"   ;
  mmdb::cpstr  asm_params_file      = "asm_params.dat"     ;
  mmdb::cpstr  rcsb_symops_file     = "rcsb_symops.dat"    ;
  mmdb::cpstr  agents_file          = "agents.dat"         ;
  mmdb::cpstr  intfstats_file       = "intfstats.dat"      ;

  ///  session files
  mmdb::cpstr  params_file          = "_pisa_params"       ;
  mmdb::cpstr  interface_file       = "_pisa_interfaces"   ;
  mmdb::cpstr  structure_file       = "_pisa_structure"    ;
  mmdb::cpstr  assembly_file        = "_pisa_assemblies"   ;
  mmdb::cpstr  rasmol_file          = "_pisa_rasmol.pdb"   ;
  mmdb::cpstr  ccp4mg_file          = "_pisa_ccp4mg.py"    ;

  ///  User-defined data names
  mmdb::cpstr  hydrogen_udd         = "hydrogen_udd"       ;

  ///  "Exclude all" key for ligand specification
  mmdb::cpstr lig_excl_all_key      = "(all)"              ;
  mmdb::cpstr lig_excl_agents_key   = "(agents)"           ;


/*
  ///  XML markup for assembly lists

  mmdb::cpstr xml_pdb_entry         = "pdb_entry"          ;
  mmdb::cpstr xml_pdb_code          = "pdb_code"           ;


  ///  XML markup for interfaces

  mmdb::cpstr xml_ints_open          = "pisa_interfaces"   ;
  mmdb::cpstr xml_ints_nints         = "n_interfaces"      ;
  mmdb::cpstr xml_interface          = "interface"         ;
  mmdb::cpstr xml_interface_id       = "id"                ;
  mmdb::cpstr xml_interface_type     = "type"              ;
  mmdb::cpstr xml_interface_nocc     = "n_occ"             ;
  mmdb::cpstr xml_interface_pval     = "pvalue"            ;
  mmdb::cpstr xml_interface_sten     = "stab_en"           ;
  mmdb::cpstr xml_interface_css      = "css"               ;
  mmdb::cpstr xml_interface_overlap  = "overlap"           ;
  mmdb::cpstr xml_interface_xrel     = "x-rel"             ;
  mmdb::cpstr xml_interface_fixed    = "fixed"             ;
  mmdb::cpstr xml_molecule_id        = "id"                ;
  mmdb::cpstr xml_molecule_visual_id = "visual_id"         ;
  mmdb::cpstr xml_molecule_class     = "class"             ;
  mmdb::cpstr xml_symop_no           = "symop_no"          ;
  mmdb::cpstr xml_symop              = "symop"             ;
  mmdb::cpstr xml_cell_i             = "cell_i"            ;
  mmdb::cpstr xml_cell_j             = "cell_j"            ;
  mmdb::cpstr xml_cell_k             = "cell_k"            ;
  mmdb::cpstr xml_interface_natoms   = "int_natoms"        ;
  mmdb::cpstr xml_interface_nres     = "int_nres"          ;
  mmdb::cpstr xml_interface_area     = "int_area"          ;
  mmdb::cpstr xml_interface_se       = "int_solv_en"       ;
  mmdb::cpstr xml_interface_res      = "residues"          ;
  mmdb::cpstr xml_residue            = "residue"           ;
  mmdb::cpstr xml_residue_serno      = "ser_no"            ;
  mmdb::cpstr xml_residue_name       = "name"              ;
  mmdb::cpstr xml_residue_seqnum     = "seq_num"           ;
  mmdb::cpstr xml_residue_inscode    = "ins_code"          ;
  mmdb::cpstr xml_residue_bonds      = "bonds"             ;
  mmdb::cpstr xml_residue_asa        = "asa"               ;
  mmdb::cpstr xml_residue_bsa        = "bsa"               ;
  mmdb::cpstr xml_residue_se         = "solv_en"           ;
  mmdb::cpstr xml_mol_weight         = "weight"            ;
  mmdb::cpstr xml_bond_nbonds        = "n_bonds"           ;
  mmdb::cpstr xml_bond_bond          = "bond"              ;
  mmdb::cpstr xml_bond_chain1        = "chain-1"           ;
  mmdb::cpstr xml_bond_label_asym_id1 = "label_asym_id-1"           ; 
  mmdb::cpstr xml_bond_orig_label_asym_id1 = "orig_label_asym_id-1"           ;
  mmdb::cpstr xml_bond_pdbx_sifts_xref_db_num1 = "pdbx_sifts_xref_db_num-1"           ; //GDL:uniprot
  mmdb::cpstr xml_bond_label_asym_id2 = "label_asym_id-2"           ;                                                               mmdb::cpstr mmdb::cpstr  xml_bond_orig_label_asym_id2 = "orig_label_asym_id-2"           ;
 mmdb::cpstr xml_bond_pdbx_sifts_xref_db_num2 = "pdbx_sifts_xref_db_num-2"           ; 
  mmdb::cpstr xml_bond_res1          = "res-1"             ;
  mmdb::cpstr xml_bond_seqnum1       = "seqnum-1"          ;
  mmdb::cpstr xml_bond_inscode1      = "inscode-1"         ;
  mmdb::cpstr xml_bond_atname1       = "atname-1"          ;
  mmdb::cpstr xml_bond_chain2        = "chain-2"           ;
  mmdb::cpstr xml_bond_res2          = "res-2"             ;
  mmdb::cpstr xml_bond_seqnum2       = "seqnum-2"          ;
  mmdb::cpstr xml_bond_inscode2      = "inscode-2"         ;
  mmdb::cpstr xml_bond_atname2       = "atname-2"          ;
  mmdb::cpstr xml_bond_dist          = "dist"              ;
  mmdb::cpstr xml_bond_hbonds        = "h-bonds"           ;
  mmdb::cpstr xml_bond_sbridges      = "salt-bridges"      ;
  mmdb::cpstr xml_bond_dsbonds       = "ss-bonds"          ;
  mmdb::cpstr xml_bond_covbonds      = "cov-bonds"         ;


  ///  XML markup for assembly lists

  mmdb::cpstr xml_asm_open           = "pisa_multimers"    ;
  mmdb::cpstr xml_pisa_results       = "pisa_results"      ;
  mmdb::cpstr xml_name               = "name"              ;
  mmdb::cpstr xml_status             = "status"            ;
  mmdb::cpstr xml_status_desc        = "status_description";
  mmdb::cpstr xml_status_note        = "status_note"       ;
  mmdb::cpstr xml_total_asm          = "total_asm"         ;
  mmdb::cpstr xml_assessment         = "assessment"        ;
  mmdb::cpstr xml_mult_state         = "multimeric_state"  ;
  mmdb::cpstr xml_asm_set            = "asm_set"           ;
  mmdb::cpstr xml_asu_complex        = "asu_complex"       ;
  mmdb::cpstr xml_asm_set_serno      = "ser_no"            ;
  mmdb::cpstr xml_assembly           = "assembly"          ;
  mmdb::cpstr xml_asm_ser_no         = "serial_no"         ;
  mmdb::cpstr xml_asm_id             = "id"                ;
  mmdb::cpstr xml_asm_size           = "size"              ;
  mmdb::cpstr xml_asm_mmsize         = "mmsize"            ;
  mmdb::cpstr xml_asm_freesize       = "freesize"          ;
  mmdb::cpstr xml_asm_diss_energy    = "diss_energy"       ;
  mmdb::cpstr xml_asm_diss_energy_0  = "diss_energy_0"     ;
  mmdb::cpstr xml_asm_asa            = "asa"               ;
  mmdb::cpstr xml_asm_bsa            = "bsa"               ;
  mmdb::cpstr xml_asm_entropy        = "entropy"           ;
  mmdb::cpstr xml_asm_entropy_0      = "entropy_0"         ;
  mmdb::cpstr xml_asm_diss_area      = "diss_area"         ;
  mmdb::cpstr xml_asm_int_dg         = "int_energy"        ;
  mmdb::cpstr xml_asm_nuc            = "n_uc"              ;
  mmdb::cpstr xml_asm_ndiss          = "n_diss"            ;
  mmdb::cpstr xml_asm_sym_num        = "symNumber"         ;
  mmdb::cpstr xml_asm_score          = "score"             ;
  mmdb::cpstr xml_asm_formula        = "formula"           ;
  mmdb::cpstr xml_asm_composition    = "composition"       ;
  mmdb::cpstr xml_orig_chains        = "all_chains_at_identity";

  mmdb::cpstr xml_molecule           = "molecule"          ;
  mmdb::cpstr xml_molecule_name      = "chain_id"          ;
  mmdb::cpstr xml_molecule_rxx       = "rxx"               ;
  mmdb::cpstr xml_molecule_rxy       = "rxy"               ;
  mmdb::cpstr xml_molecule_rxz       = "rxz"               ;
  mmdb::cpstr xml_molecule_tx        = "tx"                ;
  mmdb::cpstr xml_molecule_ryx       = "ryx"               ;
  mmdb::cpstr xml_molecule_ryy       = "ryy"               ;
  mmdb::cpstr xml_molecule_ryz       = "ryz"               ;
  mmdb::cpstr xml_molecule_ty        = "ty"                ;
  mmdb::cpstr xml_molecule_rzx       = "rzx"               ;
  mmdb::cpstr xml_molecule_rzy       = "rzy"               ;
  mmdb::cpstr xml_molecule_rzz       = "rzz"               ;
  mmdb::cpstr xml_molecule_tz        = "tz"                ;

  mmdb::cpstr xml_molecule_frxx      = "rxx-f"             ;
  mmdb::cpstr xml_molecule_frxy      = "rxy-f"             ;
  mmdb::cpstr xml_molecule_frxz      = "rxz-f"             ;
  mmdb::cpstr xml_molecule_ftx       = "tx-f"              ;
  mmdb::cpstr xml_molecule_fryx      = "ryx-f"             ;
  mmdb::cpstr xml_molecule_fryy      = "ryy-f"             ;
  mmdb::cpstr xml_molecule_fryz      = "ryz-f"             ;
  mmdb::cpstr xml_molecule_fty       = "ty-f"              ;
  mmdb::cpstr xml_molecule_frzx      = "rzx-f"             ;
  mmdb::cpstr xml_molecule_frzy      = "rzy-f"             ;
  mmdb::cpstr xml_molecule_frzz      = "rzz-f"             ;
  mmdb::cpstr xml_molecule_ftz       = "tz-f"              ;

  mmdb::cpstr xml_molecule_sym_id    = "symId"             ;

*/

  ///  TRconst is gas constant multiplied by temperature (300K)
  /// in kcal/mol: 300.0*1.9872/1000.0 = 0.59616
  const mmdb::realtype TRconst = 0.59616;

}  // namespace pisa
