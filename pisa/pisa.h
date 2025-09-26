// $Id: pisa.h $
// =================================================================
//
//    04.03.16   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -------------------------------------------------------------------
//
//  **** Module  :  PISA <main program definitions>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2016
//
// =================================================================
//

#ifndef  __PISA__
#define  __PISA__


///  Job control tags
#define  job_help_tag                  "-help"
#define  job_logo_tag                  "-version"
#define  job_analyse_tag               "-analyse"
#define  job_lig_auto_tag              "--lig=auto"
#define  job_lig_fixed_tag             "--lig=fixed"
#define  job_lig_free_tag              "--lig=free"
#define  job_lig_exclude_tag           "--lig-exclude="
//GDL: add flag 'as is'
#define  job_as_is_on_tag              "--as-is"
//#define  job_lig_excl_all_tag          "(all)"
//#define  job_lig_excl_agents_tag       "(agents)"
#define  job_erase_tag                 "-erase"
#define  job_list_tag                  "-list"
#define  job_interfaces_tag            "interfaces"
#define  job_monomers_tag              "monomers"
#define  job_assemblies_tag            "assemblies"
#define  job_view_tag                  "-view"
#define  job_stock_tag                 "stock"
#define  job_mg_tag                    "-mg"
#define  job_download_tag              "-download"
#define  job_pdb_tag                   "-pdb"
#define  job_cif_tag                   "-cif"
#define  job_detail_tag                "-detail"
#define  job_dimer_special_tag         "-dimer_special"
#define  job_dimer_loop_tag            "-dimer_loop"
#define  job_interface_tag             "interface"
#define  job_monomer_tag               "monomer"
#define  job_assembly_tag              "assembly"
#define  job_dissociate_tag            "dissociate"
#define  job_350_tag                   "-350"
#define  job_350n_tag                  "-350n"
#define  job_xml_tag                   "-xml"
#define  job_cfg_template_tag          "-cfg-template"

///  Job control ids
enum JOB_CONTROL_ID  {
  JOB_None,
  JOB_Help,
  JOB_Analyse,
  JOB_Erase,
  JOB_Logo,
  JOB_List_Interfaces,
  JOB_List_Monomers,
  JOB_List_Assemblies,
  JOB_List_Stock,
  JOB_View_Interface,
  JOB_View_Monomer,
  JOB_View_Assembly,
  JOB_View_Dissociate,
  JOB_MG_Interface,
  JOB_MG_Monomer,
  JOB_MG_Assembly,
  JOB_MG_Dissociate,
  JOB_Download_Interface,
  JOB_Download_Monomer,
  JOB_Download_Assembly,
  JOB_Download_Dissociate,
  JOB_CIF_Interface,
  JOB_CIF_Monomer,
  JOB_CIF_Assembly,
  JOB_CIF_Dissociate,
  JOB_Detail_Interface,
  JOB_Detail_Monomer,
  JOB_Detail_Assembly,
  JOB_Remark_350n,
  JOB_Remark_350,
  JOB_Make_XML_Interfaces,
  JOB_Make_XML_Assemblies,
  JOB_Cfg_Template,
  JOB_DimerSpecial,
  JOB_DimerLoop
};


#endif

