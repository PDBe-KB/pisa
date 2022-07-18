// $Id: pisa_defs.h $
// =================================================================
//
//    29.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -------------------------------------------------------------------
//
//  **** Module  :  PISA_Defs <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2016
//
// =================================================================
//

#ifndef  __PISA_Defs__
#define  __PISA_Defs__

// THIS FILE SHOULD BE NEVER USED IN HEADERS, USE PISA_TYPES.H INSTEAD

///  Platform. This should be __DEC also for SUN and LINUX

#define __DEC

///  Base download URL (used only in QtPISA)

#define base_pdb_url "http://files.rcsb.org/download/"

///  Rasmol colours

#define  COLOUR_domain1     "[150,150,150]"
#define  COLOUR_surface1    "[0,0,230]"
#define  COLOUR_interface1  "[230,0,0]"
#define  COLOUR_domain2     "[230,230,230]"
#define  COLOUR_surface2    "[0,230,230]"
#define  COLOUR_interface2  "[0,230,0]"


///  XML markup for assembly lists

#define xml_pdb_entry         "pdb_entry"
#define xml_pdb_code          "pdb_code"


///  XML markup for interfaces

#define xml_ints_open          "pisa_interfaces"
#define xml_ints_nints         "n_interfaces"
#define xml_interface          "interface"
#define xml_interface_id       "id"
#define xml_interface_type     "type"
#define xml_interface_nocc     "n_occ"
#define xml_interface_pval     "pvalue"
#define xml_interface_sten     "stab_en"
#define xml_interface_css      "css"
#define xml_interface_overlap  "overlap"
#define xml_interface_xrel     "x-rel"
#define xml_interface_fixed    "fixed"
#define xml_asmunit_id         "id"
#define xml_asmunit_visual_id  "visual_id"
#define xml_asmunit_class      "class"
#define xml_symop_no           "symop_no"
#define xml_symop              "symop"
#define xml_cell_i             "cell_i"
#define xml_cell_j             "cell_j"
#define xml_cell_k             "cell_k"
#define xml_interface_natoms   "int_natoms"
#define xml_interface_nres     "int_nres"
#define xml_interface_area     "int_area"
#define xml_interface_se       "int_solv_en"
#define xml_interface_res      "residues"
#define xml_residue            "residue"
#define xml_residue_serno      "ser_no"
#define xml_residue_name       "name"
#define xml_residue_seqnum     "seq_num"
#define xml_residue_label_seqnum   "label_seq_num"
#define xml_residue_inscode    "ins_code"
#define xml_residue_bonds      "bonds"
#define xml_residue_asa        "asa"
#define xml_residue_bsa        "bsa"
#define xml_residue_se         "solv_en"
#define xml_mol_weight         "weight"
#define xml_bond_nbonds        "n_bonds"
#define xml_bond_bond          "bond"
#define xml_bond_chain1        "chain-1"
#define xml_bond_label_asym_id1 "label_asym_id-1" //GDL:added label
#define xml_bond_atom_site_id1 "atom_site_id-1" //GDL:added atom site id
#define xml_bond_orig_label_asym_id1 "orig_label_asym_id-1" //GDL:added orig_label_asym
#define xml_bond_pdbx_sifts_xref_db_acc_1 "pdbx_sifts_xref_db_acc-1" //GDL:added uniprot acc
#define xml_bond_pdbx_sifts_xref_db_num_1 "pdbx_sifts_xref_db_num-1" //GDL:added uniprot num flag
#define xml_bond_pdbx_sifts_xref_db_name_1 "pdbx_sifts_xref_db_name-1" //GDL:added orig_label_asym
#define xml_bond_res1          "res-1"
#define xml_bond_seqnum1       "seqnum-1"
#define xml_bond_label_seqnum1 "label_seqnum-1"
#define xml_bond_inscode1      "inscode-1"
#define xml_bond_atname1       "atname-1"
#define xml_bond_chain2        "chain-2"
#define xml_bond_atom_site_id2 "atom_site_id-2" //GDL:added atom site id
#define xml_bond_label_asym_id2 "label_asym_id-2" //GDL:added label
#define xml_bond_orig_label_asym_id2 "orig_label_asym_id-2" //GDL:added remapping_orig_label_asym label
#define xml_bond_pdbx_sifts_xref_db_acc_2 "pdbx_sifts_xref_db_acc-2" //GDL:added uniprot acc 
#define xml_bond_pdbx_sifts_xref_db_num_2 "pdbx_sifts_xref_db_num-2" //GDL:added orig_label_asym
#define xml_bond_pdbx_sifts_xref_db_name_2 "pdbx_sifts_xref_db_name-2" //GDL:added orig_label_asym
#define xml_bond_label_seqnum2 "label_seqnum-2"
#define xml_bond_res2          "res-2"
#define xml_bond_seqnum2       "seqnum-2"
#define xml_bond_inscode2      "inscode-2"
#define xml_bond_atname2       "atname-2"
#define xml_bond_dist          "dist"
#define xml_bond_hbonds        "h-bonds"
#define xml_bond_sbridges      "salt-bridges"
#define xml_bond_dsbonds       "ss-bonds"
#define xml_bond_covbonds      "cov-bonds"
#define xml_bond_otherbonds    "other-bonds"


///  XML markup for assembly lists

#define xml_asm_open           "pisa_multimers"
#define xml_pisa_results       "pisa_results"
#define xml_name               "name"
#define xml_status             "status"
#define xml_status_desc        "status_description"
#define xml_status_note        "status_note"
#define xml_total_asm          "total_asm"
#define xml_assessment         "assessment"
#define xml_mult_state         "multimeric_state"
#define xml_cryst_split        "asm_set"
#define xml_asu_complex        "asu_complex"
#define xml_cryst_split_serno  "ser_no"
#define xml_assembly           "assembly"
#define xml_asm_ser_no         "serial_no"
#define xml_asm_id             "id"
#define xml_asm_size           "size"
#define xml_asm_mmsize         "mmsize"
#define xml_asm_freesize       "freesize"
#define xml_asm_diss_energy    "diss_energy"
#define xml_asm_diss_energy_0  "diss_energy_0"
#define xml_asm_asa            "asa"
#define xml_asm_bsa            "bsa"
#define xml_asm_entropy        "entropy"
#define xml_asm_entropy_0      "entropy_0"
#define xml_asm_diss_area      "diss_area"
#define xml_asm_int_dg         "int_energy"
#define xml_asm_nuc            "n_uc"
#define xml_asm_ndiss          "n_diss"
#define xml_asm_sym_num        "symNumber"
#define xml_asm_score          "score"
#define xml_asm_formula        "formula"
#define xml_asm_composition    "composition"
#define xml_asm_interfaces     "interfaces"
#define xml_asm_nints          "n_interfaces"
#define xml_asm_interface      "interface"
#define xml_asm_intf_id        "id"
#define xml_asm_intf_diss      "dissociates"
#define xml_orig_chains        "all_chains_at_identity"

#define xml_asmunit            "molecule"
#define xml_asmunit_name       "chain_id"

//GDL:commenting out info we don't need in the xml
//#define xml_asmunit_rxx        "rxx"
//#define xml_asmunit_rxy        "rxy"
//#define xml_asmunit_rxz        "rxz"
//#define xml_asmunit_tx         "tx"
//#define xml_asmunit_ryx        "ryx"
//#define xml_asmunit_ryy        "ryy"
//#define xml_asmunit_ryz        "ryz"
//#define xml_asmunit_ty         "ty"
//#define xml_asmunit_rzx        "rzx"
//#define xml_asmunit_rzy        "rzy"
//#define xml_asmunit_rzz        "rzz"
//#define xml_asmunit_tz         "tz"

#define xml_asmunit_frxx       "rxx-f"
#define xml_asmunit_frxy       "rxy-f"
#define xml_asmunit_frxz       "rxz-f"
#define xml_asmunit_ftx        "tx-f"
#define xml_asmunit_fryx       "ryx-f"
#define xml_asmunit_fryy       "ryy-f"
#define xml_asmunit_fryz       "ryz-f"
#define xml_asmunit_fty        "ty-f"
#define xml_asmunit_frzx       "rzx-f"
#define xml_asmunit_frzy       "rzy-f"
#define xml_asmunit_frzz       "rzz-f"
#define xml_asmunit_ftz        "tz-f"

#define xml_asmunit_sym_id     "symId"


#endif

