// $Id: pisa_types.h $
// =================================================================
//
//    15.03.19   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -------------------------------------------------------------------
//
//  **** Module  :  PISA_Types <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2019
//
// =================================================================
//

#ifndef  __PISA_Types__
#define  __PISA_Types__

#include "mmdb2/mmdb_mattype.h"

///  Platform. This should be __DEC also for SUN and LINUX
#define __DEC

namespace pisa  {

  enum PISALIB_VERSION  {
    MAJOR_VERSION = 2,  //!< major version
    MINOR_VERSION = 1,  //!< minor version
    MICRO_VERSION = 2   //!< micro version
  };

  extern mmdb::cpstr  BuildDate ;
  extern mmdb::cpstr  auth_email;

  ///   PISA configuration file must be either specified in the command
  /// or provided by environment variable pisa_conf_file. All further
  /// configuration is read from this configuration file. It must be
  /// an absolute path, i.e. starting from '/'.
  extern mmdb::cpstr  conf_file         ;

  extern mmdb::cpstr  ccp4_env          ;
  extern mmdb::cpstr  conf_path_ccp4    ;

  ///  Result codes
  enum RESULT_CODE  {
    RESULT_Stopped               = -2,
    RESULT_noResults             = -1,
    RESULT_Ok                    =  0,
    RESULT_Instructions          =  1,
    RESULT_Logo                  =  2,
    RESULT_ConfigurationError    =  3,
    RESULT_SessionDirFailure     =  4,
    RESULT_noSRS                 =  5,
    RESULT_noMolRefIndex         =  6,
    RESULT_noMolRefData          =  7,
    RESULT_noAgentsData          =  8,
    RESULT_notAnalysed           =  9,
    RESULT_noAtoms               = 10,
    RESULT_noCoordinateData      = 11,
    RESULT_noMacromolecules      = 12,
    RESULT_notReadyForContacts   = 13,
    RESULT_noSpaceGroup          = 14,
    RESULT_noRCSBSpaceGroup      = 15,
    RESULT_noSymOps              = 16,
    RESULT_noSymOp               = 17,
    RESULT_SurfaceCalcsFailed    = 18,
    RESULT_badCoordinateFile     = 19,
    RESULT_cantWriteResults      = 20,
    RESULT_cantReadResults       = 21,
    RESULT_SRSFailed             = 22,
    RESULT_InterfaceCalcFailure  = 23,
    RESULT_noRCSBData            = 24,
    RESULT_noRCSBassignment      = 25,
    RESULT_noAssemblyParameters  = 26,
    RESULT_noInterfaceStats      = 27,
    RESULT_unknownInterfaceType  = 28,
    RESULT_wrongRefParameter     = 29,
    RESULT_incompleteStatSet     = 30,
    RESULT_unknownIntfStatError  = 31,
    RESULT_SessionDoesntExist    = 32,
    RESULT_noSessionResults      = 33,
    RESULT_cantWriteFile         = 34,
    RESULT_SessionErased         = 35,
    RESULT_monomerNoOutOfRange   = 36,
    RESULT_interfaceNoOutOfRange = 37,
    RESULT_assemblyNoOutOfRange  = 38,
    RESULT_cantWriteRasmolData   = 39,
    RESULT_cantWriteMGData       = 40,
    RESULT_noDimerFound          = 41
  };


  ///  Assembler return codes
  enum ASSMB_RC  {
    ASSMB_Ok                    =  0,
    ASSMB_Void                  =  1,
    ASSMB_incompleteData        =  2,
    ASSMB_noSymOps              =  3,
    ASSMB_noDomains             =  4,
    ASSMB_noInterfaces          =  5,
    ASSMB_Overlap               =  6,
    ASSMB_noSymOp               =  7,
    ASSMB_improperSymOp         =  8,
    ASSMB_noNCSOp               =  9,
    ASSMB_brokenComposition     = 10,
    ASSMB_brokenComplementarity = 11,
    ASSMB_repeatedAssignment    = 12,
    ASSMB_assemblyOverflow      = 13,
    ASSMB_tooBigSystem          = 14,
    ASSMB_timeLimit             = 15
  };

  ///  PISA configuration status constants
  enum CONFIGURATION_STATUS  {
    CFG_Configured       =  0,
    CFG_Unconfigured     =  1,
    CFG_noConfFile       =  2,
    CFG_cantOpenFile     =  3,
    CFG_noDataRoot       =  4,
    CFG_noSessionPrefix  =  5,
    CFG_noSRSDir         =  6,
    CFG_noMolRefDir      =  7,
    CFG_noPIStoreDir     =  8,
    CFG_noRasmolCommand  =  9,
    CFG_noJmolCommand    = 10,
    CFG_noCCP4MGCommand  = 11
  };

  ///  prefix for session directories
  extern mmdb::cpstr  default_session_prefix;

  /// session directory return codes
  enum SDIR_RC  {
    SDIR_Created      = 0,
    SDIR_Exists       = 1,
    SDIR_Fail         = 2,
    SDIR_doesntExist  = 3,
    SDIR_noResults    = 4
  };

  ///  configuration file tags
  extern mmdb::cpstr  data_root_tag      ;
  extern mmdb::cpstr  srs_dir_tag        ;
  extern mmdb::cpstr  molref_dir_tag     ;
  extern mmdb::cpstr  pistore_dir_tag    ;
  extern mmdb::cpstr  help_dir_tag       ;
  extern mmdb::cpstr  dnl_url_tag        ;
  extern mmdb::cpstr  rasmol_com_tag     ;
  extern mmdb::cpstr  jmol_com_tag       ;
  extern mmdb::cpstr  ccp4mg_com_tag     ;
  extern mmdb::cpstr  session_prefix_tag ;
  extern mmdb::cpstr  php_uri_tag        ;
  extern mmdb::cpstr  help_uri_tag       ;
  extern mmdb::cpstr  jsrview_uri_tag    ;
  extern mmdb::cpstr  pdb_dir_tag        ;
  extern mmdb::cpstr  pdb_dir_format_tag ;
  extern mmdb::cpstr  erase_time_tag     ;
  extern mmdb::cpstr  expiry_time_tag    ;

  extern mmdb::cpstr  pdb_plain_pdb_key      ;
  extern mmdb::cpstr  pdb_plain_pdb_gz_key   ;
  extern mmdb::cpstr  pdb_plain_mmcif_key    ;
  extern mmdb::cpstr  pdb_plain_mmcif_gz_key ;
  extern mmdb::cpstr  pdb_split_pdb_key      ;
  extern mmdb::cpstr  pdb_split_pdb_gz_key   ;
  extern mmdb::cpstr  pdb_split_mmcif_key    ;
  extern mmdb::cpstr  pdb_split_mmcif_gz_key ;

  ///  Store files
  extern mmdb::cpstr  syminfo_file       ;
  extern mmdb::cpstr  asm_params_file    ;
  extern mmdb::cpstr  rcsb_symops_file   ;
  extern mmdb::cpstr  agents_file        ;
  extern mmdb::cpstr  intfstats_file     ;

  /// dictionary codes
  enum DICT_RC  {
    DICT_Ok            = 0,
    DICT_noSRS         = 1,
    DICT_noMolRefIndex = 2,
    DICT_noMolRefData  = 3,
    DICT_noAgentsData  = 4,
    DICT_noIntfStats   = 5
  };

  ///  Data analysis status
  typedef mmdb::word DATASTATUS;
  enum DATASTATUS_FLAG  {
    DS_Ok               = 0x00000000,
    DS_notAnalysed      = 0x00000001,
    DS_noAtoms          = 0x00000002,
    DS_noChains         = 0x00000004,
    DS_noMacromolecules = 0x00000008,
    DS_noCrystData      = 0x00000010,
    DS_noCellData       = 0x00000020,
    DS_noSpaceGroup     = 0x00000040,
    DS_noSymOp          = 0x00000080,
    DS_badSymOp         = 0x00000100,
    DS_badFileName      = 0x00000200,
    DS_noMolRef         = 0x00000400,
    DS_badFile          = 0x00010000
  };

  ///  session files
  extern mmdb::cpstr  params_file    ;
  extern mmdb::cpstr  interface_file ;
  extern mmdb::cpstr  structure_file ;
  extern mmdb::cpstr  assembly_file  ;
  extern mmdb::cpstr  rasmol_file    ;
  extern mmdb::cpstr  ccp4mg_file    ;

  ///  Domain classes
  enum DOMAIN_CLASS  {
    DCLASS_None    = 0,
    DCLASS_Protein = 1,
    DCLASS_DNA     = 2,
    DCLASS_RNA     = 3,
    DCLASS_Ligand  = 4
  };

  ///  User-defined data names
  extern mmdb::cpstr  hydrogen_udd  ;

  ///  Viewer target keys
  enum VIEWER_TARGET  {
    TARGET_NONE     = -1,
    TARGET_Download =  0,
    TARGET_CIF      =  1,
    TARGET_Rasmol   =  2,
    TARGET_CCP4MG   =  3
  };

  ///  Ligand treatment keys
  enum LIGAND_KEY  {
    LIGANDS_Auto    = 0,
    LIGANDS_FreeAll = 1,
    LIGANDS_FixAll  = 2
  };

  enum AS_IS_KEY  {
    AS_IS_off    = 0,
    AS_IS_on     = 1
  };



  enum ISORT_KEY  {
    ISORT_Off    = 0,
    ISORT_Area   = 1,
    ISORT_DeltaG = 2,
    ISORT_IType  = 3
  };

  /// Interface status
  typedef mmdb::word INTERFACE_STATUS;
  enum INTERFACE_STATUS_FLAG  {
    INTS_Ok            = 0x00,
    INTS_notCalculated = 0x01,
    INTS_noInterfaces  = 0x02,
    INTS_mmOverlap     = 0x04,
    INTS_ligandOverlap = 0x08
  };

  ///  "Exclude all" key for ligand specification
  extern mmdb::cpstr lig_excl_all_key                ;
  extern mmdb::cpstr lig_excl_agents_key             ;

  ///  TRconst is gas constant multiplied by temperature (300K)
  /// in kcal/mol: 300.0*1.9872/1000.0 = 0.59616
  extern const mmdb::realtype TRconst;

  enum PDBDIR_FORMAT  {
    PDB_Unspecified,
    PDB_PLAIN_PDB,
    PDB_PLAIN_PDB_GZ,
    PDB_PLAIN_mmCIF,
    PDB_PLAIN_mmCIF_GZ,
    PDB_SPLIT_PDB,
    PDB_SPLIT_PDB_GZ,
    PDB_SPLIT_mmCIF,
    PDB_SPLIT_mmCIF_GZ
  };

}  // namespace pisa

#endif
