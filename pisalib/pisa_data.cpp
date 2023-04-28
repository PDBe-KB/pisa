// $Id: pisa_data.cpp $
// =================================================================
//
//    14.03.19   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_data <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Data
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2019
//
// =================================================================
//

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#include<iostream>
using namespace std;

#include "pisa_data.h"
#include "pisa_types.h"
#include "pisa_defs.h"

#if defined(_MSC_VER) || defined (WIN32)
 #include <io.h>
 #include <direct.h>
#else
 #include <unistd.h>
 #include <dirent.h>
#endif

#include "mmdb2/mmdb_tables.h"
#include "ccp4srs/ccp4srs_defs.h"

namespace pisa  {

  // ========================  Data  ============================

  Data::Data ( mmdb::cpstr confPath )  {

    InitData();

    if (!confPath)                 readConfiguration ( confPath );
    else if (strcmp(confPath,"!")) readConfiguration ( confPath );

  }

  Data::~Data()  {

    FreeMemory();

    //  Release server configuration parameters
    if (syminfo_lib)    delete[] syminfo_lib;
    if (cfgPath)        delete[] cfgPath;
    if (dataRoot)       delete[] dataRoot;
    if (srsDir)         delete[] srsDir;
    if (MolRefDir)      delete[] MolRefDir;
    if (PIStoreDir)     delete[] PIStoreDir;
    if (helpDir)        delete[] helpDir;
    if (dnlUrl)         delete[] dnlUrl;
    if (rasmol_com)     delete[] rasmol_com;
    if (jmol_com)       delete[] jmol_com;
    if (ccp4mg_com)     delete[] ccp4mg_com;
    if (session_prefix) delete[] session_prefix;
    if (phpUri)         delete[] phpUri;
    if (helpUri)        delete[] helpUri;
    if (jsrviewUri)     delete[] jsrviewUri;
    if (pdbDir)         delete[] pdbDir;

    if (SRS)       delete   SRS;
    if (molRef)    delete   molRef;
    if (agentsRef) delete[] agentsRef;
    if (proSurf)   delete   proSurf;

    FreeIntfStats();

  }

  void Data::FreeMemory()  {

    if (sessionID)  delete[] sessionID;
    if (crDir)      delete[] crDir;
    if (query)      delete   query;

    query       = NULL;
    crDir       = NULL;
    sessionID   = NULL;

  }

  void Data::FreeIntfStats()  {
    if (intfStats)  {
      for (int i=0;i<nIntfStats;i++)
        if (intfStats[i])  delete intfStats[i];
      delete[] intfStats;
      intfStats = NULL;
    }
    nIntfStats = 0;
  }


  void Data::InitData()  {

    //  data references
    SRS              = NULL;   // common SRS resource
    molRef           = NULL;   // common molecular reference
    agentsRef        = NULL;   // crystallization agents reference
    proSurf          = NULL;   // surface calculation module
    intfStats        = NULL;   // interface statistics for radar plots
    nIntfStats       = 0;      // number of interface parameters

    //  Program configuration
    confStatus       = CFG_Unconfigured;  // configuration status
    syminfo_lib      = NULL;   // path to CCP4 symmetry information file
    cfgPath          = NULL;   // path to Configuration File
    dataRoot         = NULL;   // general data root
    srsDir           = NULL;   // path to directory with SRS files
    MolRefDir        = NULL;   // path to directory with MolRef files
    PIStoreDir       = NULL;   // path to directory with pistore files
    helpDir          = NULL;   // path to directory with help files
    dnlUrl           = NULL;   // Url to download files
    rasmol_com       = NULL;   // rasmol-launching command
    jmol_com         = NULL;   // jmol-launching command
    ccp4mg_com       = NULL;   // CCP4-MG launch command

    //  Session configuration
    session_prefix = NULL;      // prefix for session directories
    mmdb::CreateCopy ( session_prefix,default_session_prefix );
    sessionID      = NULL;      // unique session ID
    crDir          = NULL;      // session directory

    //  Server configuration
    phpUri         = NULL;       // URI to directory with php scripts
    helpUri        = NULL;       // URI to directory with help files
    mmdb::CreateCopy ( phpUri,"http://www.pisa.com/php/" );
    jsrviewUri     = NULL;       // URI to jsrview directory
    mmdb::CreateCopy ( jsrviewUri,"http://www.pisa.com/jsrview/" );
    pdbDir         = NULL;       // path to PDB directory
    mmdb::CreateCopy ( pdbDir,"" );
    pdbFormat      = PDB_Unspecified; // file format in PDB directory
    expiry_time    = 3*3600;     // time for cleaning the storage, secs
    erase_time     = 48*3600;    // time for unconditional cleaning, secs

    //  Structure data
    query       = new QueryData();
    nCellLayers = 3;           // number of cell layers to check
    nCellOut    = 5;           // num. of cell layers to output in symops

    MMDBErrLine[0] = char(0);
    MMDBErrLineNo  = 0;

//    resultCode = RESULT_noResults;

    //  Interface calculation parameters
    sphCodeNo    = 36;     // spherical code number
    probe_radius = 1.4;    // radius of probe sphere (1.4 A for water)

  }

  bool Data::getConfPath ( mmdb::io::RFile f,
                           mmdb::pstr      S,
                           int             lenS,
                           mmdb::cpstr     tag,
                           mmdb::pstr    & path,
                           bool            dir )  {
  mmdb::pstr p;

    if (!strcmp(S,tag))  {

      f.ReadLine ( S,lenS );
      mmdb::CutSpaces  ( S,mmdb::SCUTKEY_BEGEND );

  #if defined(_MSC_VER) || defined (WIN32)
      p = S;
      while (*p)  {
        if (*p=='/')  *p = '\\';
        p++;
      }
  #endif

      if (dir && (S[strlen(S)-1]!=mmdb::io::_dir_sep_c))
        strcat ( S,mmdb::io::_dir_sep );

      // check for environmental variables
      if (S[0]!='$')
         mmdb::CreateCopy ( path,S );
      else  {
        p = mmdb::FirstOccurence ( S,mmdb::io::_dir_sep_c );
        if (p)  {
          *p = char(0);
          mmdb::CreateCopCat ( path,getenv(&(S[1])),
                               mmdb::io::_dir_sep,&(p[1]) );
        } else
          mmdb::CreateCopy ( path,getenv(&(S[1])) );
      }

      return true;

    }

    return false;

  }

  bool Data::getConfLine ( mmdb::io::RFile f, mmdb::pstr S, int lenS,
                           mmdb::cpstr   tag, mmdb::pstr & line )  {
    if (!strcmp(S,tag))  {
      f.ReadLine ( S,lenS );
      mmdb::CutSpaces  ( S,mmdb::SCUTKEY_BEGEND );
      mmdb::CreateCopy ( line,S );
      return true;
    }
    return false;
  }

  bool Data::getConfInt ( mmdb::io::RFile f,
                          mmdb::pstr      S,
                          int             lenS,
                          mmdb::cpstr     tag,
                          int           & iv )  {
  mmdb::pstr     endptr;
  mmdb::realtype v;
    if (!strcmp(S,tag))  {
      f.ReadLine ( S,lenS );
      mmdb::CutSpaces  ( S,mmdb::SCUTKEY_BEGEND );
      v = strtod ( S,&endptr );
      if (endptr!=S)  {
        iv = mmdb::mround ( v );
        return true;
      }
    }
    return false;
  }

  bool Data::getConfBool ( mmdb::io::RFile f,
                           mmdb::pstr      S,
                           int             lenS,
                           mmdb::cpstr     tag,
                           bool          & B )  {
    if (!strcmp(S,tag))  {
      f.ReadLine ( S,lenS );
      mmdb::CutSpaces  ( S,mmdb::SCUTKEY_BEGEND );
      if (!strcasecmp(S,"YES"))      B = true;
      else if (!strcasecmp(S,"NO"))  B = false;
      return true;
    }
    return false;
  }

  bool Data::getConfOnOff  ( mmdb::io::RFile  f, mmdb::pstr S, int lenS,
                                     mmdb::cpstr tag, bool & B )  {
    if (!strcmp(S,tag))  {
      f.ReadLine ( S,lenS );
      mmdb::CutSpaces  ( S,mmdb::SCUTKEY_BEGEND );
      if (!strcasecmp(S,"ON"))        B = true;
      else if (!strcasecmp(S,"OFF"))  B = false;
      return true;
    }
    return false;
  }

  void Data::set_ccp4mg_path ( mmdb::cpstr mg_path )  {
    mmdb::CreateCopy ( ccp4mg_com,mg_path );
  }


  void Data::readConfiguration ( mmdb::cpstr confPath )  {
  //  Getting parameters from the configuration file
  mmdb::io::File f;
  mmdb::pstr     L;
  char           S[1000];
  int            m;

    L = NULL;
    mmdb::CreateCopy ( L,"" );

    mmdb::CreateCopy ( cfgPath,confPath );
    confStatus = CFG_Unconfigured;  // configuration status

    if (!cfgPath)
      mmdb::CreateCopy ( cfgPath,getenv(conf_file) );
                                    // may be provided by environment

    if (!cfgPath)  {
      // may be provided by environment
      mmdb::CreateCopy ( cfgPath,getenv(ccp4_env) );
      if (cfgPath)
        mmdb::CreateConcat ( cfgPath,conf_path_ccp4 );
    }

    mmdb::CreateCopy ( dataRoot,getenv(data_root_tag) );
    if (dataRoot)  {
      if (dataRoot[strlen(dataRoot)-1]!=mmdb::io::_dir_sep_c)
        mmdb::CreateConcat ( dataRoot,mmdb::io::_dir_sep );
    }

    if (cfgPath)  {

      f.assign ( cfgPath,true );
      if (f.reset(true))  {

        m = sizeof(S);

        while (!f.FileEnd())  {

          f.ReadLine ( S,m );
          mmdb::CutSpaces ( S,mmdb::SCUTKEY_BEGEND );
          if (!dataRoot)
            getConfPath ( f,S,m,data_root_tag   ,dataRoot  ,true  );
          getConfPath ( f,S,m,srs_dir_tag       ,srsDir    ,true  );
          getConfPath ( f,S,m,molref_dir_tag    ,MolRefDir ,true  );
          getConfPath ( f,S,m,pistore_dir_tag   ,PIStoreDir,true  );
          getConfPath ( f,S,m,help_dir_tag      ,helpDir   ,true  );
          getConfLine ( f,S,m,dnl_url_tag       ,dnlUrl           );
          getConfPath ( f,S,m,rasmol_com_tag    ,rasmol_com,false );
          getConfPath ( f,S,m,jmol_com_tag      ,jmol_com  ,false );
          getConfPath ( f,S,m,ccp4mg_com_tag    ,ccp4mg_com,false );
          getConfLine ( f,S,m,session_prefix_tag,session_prefix   );
          getConfPath ( f,S,m,php_uri_tag       ,phpUri     ,true );
          getConfPath ( f,S,m,help_uri_tag      ,helpUri   ,true  );
          getConfPath ( f,S,m,jsrview_uri_tag   ,jsrviewUri ,true );
          getConfPath ( f,S,m,pdb_dir_tag       ,pdbDir     ,true );

          if (getConfInt(f,S,m,expiry_time_tag,expiry_time))
            expiry_time *= 3600;  // now in secs
          if (getConfInt(f,S,m,erase_time_tag,erase_time))
            erase_time  *= 3600;  // now in secs

          if (getConfLine(f,S,m,pdb_dir_format_tag,L))  {
            if (!strcasecmp(L,pdb_plain_pdb_key))
              pdbFormat = PDB_PLAIN_PDB;
            else if (!strcasecmp(L,pdb_plain_pdb_gz_key))
              pdbFormat = PDB_PLAIN_PDB_GZ;
            else if (!strcasecmp(L,pdb_plain_mmcif_key))
              pdbFormat = PDB_PLAIN_mmCIF;
            else if (!strcasecmp(L,pdb_plain_mmcif_gz_key))
              pdbFormat = PDB_PLAIN_mmCIF_GZ;
            else if (!strcasecmp(L,pdb_split_pdb_key))
              pdbFormat = PDB_SPLIT_PDB;
            else if (!strcasecmp(L,pdb_split_pdb_gz_key))
              pdbFormat = PDB_SPLIT_PDB_GZ;
            else if (!strcasecmp(L,pdb_split_mmcif_key))
              pdbFormat = PDB_SPLIT_mmCIF;
            else if (!strcasecmp(L,pdb_split_mmcif_gz_key))
              pdbFormat = PDB_SPLIT_mmCIF_GZ;
          }

        }

        f.shut();

        mmdb::CreateCopCat  ( syminfo_lib,PIStoreDir,syminfo_file );
        SetSyminfoLib ( syminfo_lib );

        confStatus = CFG_Configured;

      } else
        confStatus = CFG_cantOpenFile;

    } else
      confStatus = CFG_noConfFile;

    if (L)  delete[] L;

  /*
    printf ( "\n CONFIGURATION:\n"
             " dataRoot       = '%s'\n"
             " SRSDir         = '%s'\n"
             " MolRefDir      = '%s'\n"
             " PIStoreDir     = '%s'\n"
             " helpDir        = '%s'\n"
             " rasmol_com     = '%s'\n"
             " jmol_com       = '%s'\n"
             " ccp4mg_com     = '%s'\n"
             " session_prefix = '%s'\n\n",
             dataRoot,srsDir,MolRefDir,PIStoreDir,helpDir,
             rasmol_com,jmol_com,ccp4mg_com,session_prefix
           );
  */

  }


  mmdb::pstr Data::getConfigurationTemplate ( mmdb::pstr & S )  {
  char N[100];

    mmdb::CreateCopCat ( S,
"# ------------------------------------------------------------------- #\n"
"#                                                                     #\n"
"#          This is configuratrion file for PISA software.             #\n"
"#                                                                     #\n"
"#   When used in command-prompt mode, this file must be specified     #\n"
"#       as last argument in thecommand line, or pointed out by        #\n"
"#              PISA_CONF_FILE environmental variable.                 #\n"
"#                                                                     #\n"
"#    This file may be also used to configure the QtPISA graphical     #\n"
"#   application, by either reading it using the \"Load CFG\" button   #\n"
"#  in PISA Settings Dialog, or by running QtPISA from command prompt  #\n"
"#    with this file as the only command-line argument. QtPISA needs   #\n"
"#           to be configure only once after installation.             #\n"
"#                                                                     #\n"
"# ------------------------------------------------------------------- #\n"
"\n"
"\n"
"#  ",data_root_tag," must point on the directory to contain session\n"
"#  directories.\n",
data_root_tag,"\n" );
    if (dataRoot && dataRoot[0])
         mmdb::CreateConcat ( S,dataRoot,"\n" );
    else mmdb::CreateConcat ( S,"/path/to/dataroot\n" );

    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  ",srs_dir_tag," must point on the directory containing SRS files.\n",
srs_dir_tag,"\n" );
    if (srsDir && srsDir[0])
          mmdb::CreateConcat ( S,srsDir,"\n" );
    else  mmdb::CreateConcat ( S,"/path/to/srsdir\n" );

    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  ",molref_dir_tag," must point on the directory containing MolRef files.\n",
molref_dir_tag,"\n" );
    if (MolRefDir && MolRefDir[0])
          mmdb::CreateConcat ( S,MolRefDir,"\n" );
    else  mmdb::CreateConcat ( S,"/path/to/molref\n" );

    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  ",pistore_dir_tag," must point on the directory containing files:\n" );
    mmdb::CreateConcat ( S,
"#          ",agents_file,"\n" );
    mmdb::CreateConcat ( S,
"#          ",asm_params_file,"\n" );
    mmdb::CreateConcat ( S,
"#          ",intfstats_file,"\n" );
    mmdb::CreateConcat ( S,
"#          ",syminfo_file," and\n" );
    mmdb::CreateConcat ( S,
"#          ",rcsb_symops_file,"\n",
pistore_dir_tag,"\n" );
    if (PIStoreDir && PIStoreDir[0])
          mmdb::CreateConcat ( S,PIStoreDir,"\n" );
    else  mmdb::CreateConcat ( S,"/path/to/pisastore\n" );

    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  ",help_dir_tag," must point on the directory containing HTML-formatted\n"
"#  help files.\n",
help_dir_tag,"\n" );
    if (helpDir && helpDir[0])
          mmdb::CreateConcat ( S,helpDir,"\n" );
    else  mmdb::CreateConcat ( S,"/path/to/pisahelp\n" );

    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  ",dnl_url_tag," is used only in QtPISA. It must give a valid base URL\n"
"#  for downloading PDB's cif-formatted files. The full downlad URL is\n"
"#  formed by appending lower-case PDB code with extension '.cif' to the\n"
"#  base download URL.\n",
dnl_url_tag,"\n" );
    if (dnlUrl && dnlUrl[0])
          mmdb::CreateConcat ( S,dnlUrl,"\n" );
    else  mmdb::CreateConcat ( S,base_pdb_url,"\n" );

    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  ",rasmol_com_tag," must give the rasmol launch command line\n",
rasmol_com_tag,"\n" );
    if (rasmol_com && rasmol_com[0])
          mmdb::CreateConcat ( S,rasmol_com,"\n" );
    else  mmdb::CreateConcat ( S,"/path/to/rasmol\n" );

    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  ",jmol_com_tag," must give path to jmol.jar\n",
jmol_com_tag,"\n" );
    if (jmol_com && jmol_com[0])
          mmdb::CreateConcat ( S,jmol_com,"\n" );
    else  mmdb::CreateConcat ( S,"/path/to/jmol.jar\n" );

    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  ",ccp4mg_com_tag," must give the ccp4mg launch command line\n",
ccp4mg_com_tag,"\n" );
    if (ccp4mg_com && ccp4mg_com[0])
         mmdb::CreateConcat ( S,ccp4mg_com,"\n" );
#ifdef Q_OS_WIN
    else mmdb::CreateConcat ( S,"/path/to/ccp4mg.exe\n" );
#endif
#ifdef Q_OS_LINUX
    else mmdb::CreateConcat ( S,"/path/to/ccp4mg\n" );
#endif
#ifdef Q_OS_MAC
    else mmdb::CreateConcat ( S,"/path/to/QtMG.app\n" );
#endif
    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  ",session_prefix_tag," is prefix for the names of session directories\n"
"#  created in DATA_PATH (\"pisrv_\" is used by default). Be sure to\n"
"#  have unique prefixes for each configuration file that is invoked\n"
"#  from a different user login or apache service. Session directories\n"
"#  are regularly removed from DATA_PATH, and SESSION_PREFIX allows\n"
"#  one to avoid permission conflicts between different services.\n" );
    mmdb::CreateConcat ( S,
session_prefix_tag,"\n",
session_prefix,"\n" );


    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  =========================================================================\n"
"#     The following configuration parameters are needed only for jspisa,\n"
"#   which is a web-server application. They may be ignored or even removed\n"
"#        from configuration files for command-prompt pisa and qtpisa.\n"
"#  =========================================================================\n"
"\n"
"\n"
"#  ",pdb_dir_tag," is used only by jspisa (web-server). It must\n"
"#  give absolute path to PDB directory.\n" );
    mmdb::CreateConcat ( S,
pdb_dir_tag,"\n",
pdbDir,"\n" );

    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  ",help_uri_tag," should point on the directory containing HTML-formatted\n"
"#  help files.\n",
help_uri_tag,"\n" );
    if (helpUri && helpUri[0])
          mmdb::CreateConcat ( S,helpUri,"\n" );
    else  mmdb::CreateConcat ( S,"http://www.domain.com/path/to/pisahelp\n" );

    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  ",pdb_dir_format_tag," is used only by jspisa (web-server). It must\n"
"#  specify whether the PDB dircetory has plain structures (all files\n"
"#  at root), or split (files found in sub-directories named by letters\n"
"#  2 and 3 of the PDB code, small case), and whether in PDB or mmCIF\n"
"#  format, gzipped or not. Permissible values include:\n"
"#      ",pdb_plain_pdb_key,"\n" );
    mmdb::CreateConcat ( S,
"#      ",pdb_plain_mmcif_key,"\n"
"#      ",pdb_plain_pdb_gz_key,"\n" );
    mmdb::CreateConcat ( S,
"#      ",pdb_plain_mmcif_gz_key,"\n"
"#      ",pdb_split_pdb_key,"\n" );
    mmdb::CreateConcat ( S,
"#      ",pdb_split_mmcif_key,"\n"
"#      ",pdb_split_pdb_gz_key,"\n" );
    mmdb::CreateConcat ( S,
"#      ",pdb_split_mmcif_gz_key,"\n",
pdb_dir_format_tag,"\n" );
    mmdb::CreateConcat ( S,pdb_plain_pdb_gz_key,"\n" );

    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  ",php_uri_tag," is used only by jspisa (web-server). It must\n"
"#  start with \"http://\" and point to directory with php scripts.\n" );
    mmdb::CreateConcat ( S,
php_uri_tag,"\n",
phpUri,"\n" );

    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  ",jsrview_uri_tag," is used only by jspisa (web-server). It must\n"
"#  start with \"http://\" and point to directory with javascript\n"
"#  support for jsrview api.\n" );
    mmdb::CreateConcat ( S,
jsrview_uri_tag,"\n",
jsrviewUri,"\n" );

    sprintf ( N,"%i",expiry_time );
    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  ",expiry_time_tag," is used only by jspisa (web-server). It gives\n"
"#  time (integer hours) for session to expire. Expired sessions are\n"
"#  erased when necessary unless they are locked.\n" );
    mmdb::CreateConcat ( S,
expiry_time_tag,"\n",
N,"\n" );

    sprintf ( N,"%i",erase_time );
    mmdb::CreateConcat ( S,
"\n"
"\n"
"#  ",erase_time_tag," is used only by jspisa (web-server). It gives\n"
"#  time (integer hours) for session to be erased unconditionally,\n"
"#  whether locked or not.\n" );
    mmdb::CreateConcat ( S,
erase_time_tag,"\n",
N,"\n" );

    return S;

  }


  CONFIGURATION_STATUS Data::ConfStatus()  {
    if (confStatus!=CFG_Configured)
                          return confStatus;
    if (!dataRoot)        return CFG_noDataRoot;
    if (!session_prefix)  return CFG_noSessionPrefix;
    if (!srsDir)          return CFG_noSRSDir;
    if (!MolRefDir)       return CFG_noMolRefDir;
    if (!PIStoreDir)      return CFG_noPIStoreDir;
//    if (!helpDir)     return CFG_noPIStoreDir;
    if (!rasmol_com)      return CFG_noRasmolCommand;
    if (!jmol_com)        return CFG_noJmolCommand;
    if (!ccp4mg_com)      return CFG_noCCP4MGCommand;
    return CFG_Configured;
  }


  DICT_RC Data::openMolRef ( bool indexOnly )  {
    if (!molRef)  {
      molRef = new MolRefIndex();
      molRef->assign ( MolRefDir );
      if (molRef->loadIndex())  {
        delete molRef;
        molRef = NULL;
        return DICT_noMolRefIndex;
      }
    }
    if ((!indexOnly) && (!molRef->isDataOpen()))  {
      if (molRef->openData())  {
        delete molRef;
        molRef = NULL;
        return DICT_noMolRefData;
      }
    }
    return DICT_Ok;
  }

  DICT_RC Data::openAgentsRef()  {
  mmdb::io::File f;
  mmdb::pstr     S,p;
  int            maxlen;

    if (!agentsRef)  {
      S = NULL;
      mmdb::CreateConcat  ( S,PIStoreDir,agents_file );
      f.assign ( S,true );
      if (f.reset(true))  {
        delete[] S;
        maxlen = 2000;
        S = new char[maxlen+2];
        while (!f.FileEnd())  {
          f.ReadNonBlankLine ( S,maxlen );
          mmdb::CutSpaces ( S,mmdb::SCUTKEY_BEGIN );
          if (S[0]!='#')  {  // filter out comments
            p = mmdb::FirstOccurence ( S,',' );
            if (p)  {
              *p = char(0);
              mmdb::DelSpaces ( S );
              mmdb::CreateConcat ( agentsRef,",",S );
            }
          }
        }
        mmdb::CreateConcat ( agentsRef,"," );
  //      printf ( " AgentsRef='%s'\n",AgentsRef );
        f.shut();
      }
      delete[] S;
    }

    if (agentsRef)  return DICT_Ok;
              else  return DICT_noAgentsData;

  }

  DICT_RC Data::readIntfStats()  {
  mmdb::io::File f;
  mmdb::pstr     S;
  bool           ok;

    if (!intfStats)  {
      S = NULL;
      mmdb::CreateConcat  ( S,PIStoreDir,intfstats_file );
      f.assign ( S,false,true );
      if (f.reset(true))  {
        URIStats.read ( f );
        f.ReadInt ( &nIntfStats );
        intfStats = new PDataStats[nIntfStats];
        ok = true;
        for (int i=0;i<nIntfStats;i++)  {
          intfStats[i] = new DataStats();
          ok = ok && intfStats[i]->read ( f );
        }
        f.shut();
        if (!ok)
          FreeIntfStats();
      }
      delete[] S;
    }

    if (intfStats)  return DICT_Ok;
              else  return DICT_noIntfStats;

  }


  DICT_RC Data::openDictionaries ( bool read_intf_stats )  {
  int     i;
  DICT_RC rc;

    // 1. If SBase is not open, open and initialize it

    if (!SRS)  {

      SRS = new ccp4srs::Manager();
      if (SRS->loadIndex(srsDir)!=ccp4srs::CCP4SRS_Ok)  {
        delete SRS;
        SRS = NULL;
        return DICT_noSRS;
      }

      // 1.2 For extra efficiency, pre-load all aminoacids and
      //     nucleotides. The names come from mmdb_tables.h

      for (i=0;i<mmdb::nAminoacidNames;i++)
        SRS->loadStructure ( mmdb::pstr(mmdb::AAProperties[i].name) );

      for (i=0;i<mmdb::nNucleotideNames;i++)
        SRS->loadStructure ( mmdb::NucleotideName[i] );

    }

    // 2. If MolRef is not open, open and initialize it

    rc = openMolRef ( false );
    if (rc!=DICT_Ok)
      return rc;

    // 3. If AgentsRef is not open, open it and read agents in

    rc = openAgentsRef();
    if (rc!=DICT_Ok)
      return rc;

    // 4. Read interface statistics if needed

    if (read_intf_stats)  rc = readIntfStats();
    return rc;

  }

  ccp4srs::CCP4SRS_RC Data::getEnergyTypes ( mmdb::PManager MMDB )  {
    if (openDictionaries(false)==DICT_Ok)
      return SRS->getEnergyTypes ( MMDB,NULL );
    return ccp4srs::CCP4SRS_Incomplete;
  }

  void Data::makeProSurf()  {
    if (!proSurf)  {
      proSurf = new ProSurf();
      proSurf->setCodeNo         ( sphCodeNo    );
      proSurf->setProbeRadius    ( probe_radius );
      proSurf->setExcludeSolvent ( true         );
    }
  }


  RESULT_CODE Data::readAssemblyParameters()  {
  mmdb::io::File f;
  mmdb::pstr     S;
  RESULT_CODE   rc;

    S = NULL;
    mmdb::CreateCopCat ( S,PIStoreDir,asm_params_file );
    f.assign ( S,false,true );
    if (f.reset(true))  {
      readEntropyFactors ( f );
      readBondFactors    ( f );
      readASPs           ( f );
      f.shut();
      rc = RESULT_Ok;
    } else
      rc = RESULT_noAssemblyParameters;

    /*
    mmdb::CreateCopCat ( S,PIStoreDir,asm_params_file,".new" );
    printf ( " S='%s'\n",S );
    f.assign ( S,false,true );
    if (f.rewrite())  {
      writeEntropyFactors ( f );
      writeBondFactors    ( f );
      writeASPs           ( f );
      f.shut();
    }
  */

    if (S)  delete[] S;

    return rc;

  }

  void Data::setSessionPrefix ( mmdb::cpstr prefix )  {
    mmdb::CreateCopy ( session_prefix,prefix );
  }

  void Data::makeCrDirPath ( mmdb::cpstr sessionName )  {
    mmdb::CreateCopy ( sessionID,sessionName );
    if (sessionID[0])
          mmdb::CreateCopCat ( crDir,dataRoot,session_prefix,sessionID,
                               mmdb::io::_dir_sep );
    else  mmdb::CreateCopy   ( crDir,dataRoot );
  }


  #if defined(_MSC_VER) || defined (WIN32)

  int Data::makeCrDir ( mmdb::cpstr sessionName )  {
  struct _finddata_t dir;
  long               hfile;
  int                rc;

    makeCrDirPath ( sessionName );
    /*
    mmdb::CreateCopy   ( sessionID,sessionName );
    mmdb::CreateCopCat ( crDir,dataRoot,session_prefix,sessionID );
    */

    hfile = _findfirst(crDir, &dir);
  //  printf("creating dir %s\n",crDir);
    if ( hfile != -1L )  {
      rc = SDIR_Exists;
      _findclose ( hfile );
    } else if (_mkdir(crDir))
      rc = SDIR_Fail;
    else
      rc = SDIR_Created;

    mmdb::CreateConcat ( crDir,mmdb::io::_dir_sep );

    return rc;

  }

  int Data::checkCrDir ( mmdb::cpstr sessionName )  {
  struct _finddata_t dir;
  int                rc;
  long               hfile;

    makeCrDirPath ( sessionName );
    /*
    mmdb::CreateCopy   ( sessionID,sessionName );
    mmdb::CreateCopCat ( crDir,dataRoot,session_prefix,sessionID );
    */

    crDir[strlen(crDir)-1] = char(0);
    hfile = _findfirst(crDir, &dir);
    
    mmdb::CreateConcat ( crDir,mmdb::io::_dir_sep );
    
    if (hfile!=-1L)  {
      rc = SDIR_Exists;
      _findclose ( hfile );
      //mmdb::CreateCopCat ( crDir,dataRoot,session_prefix,sessionID,mmdb::io::_dir_sep );
      readParams();
      if (!crDir)
        rc = SDIR_noResults;
    } else
      rc = SDIR_doesntExist;

    return rc;

  }

  #else

  int Data::makeCrDir ( mmdb::cpstr sessionName )  {
  DIR * dir;
  int   rc;

    makeCrDirPath ( sessionName );
    /*
    mmdb::CreateCopy   ( sessionID,sessionName );
    mmdb::CreateCopCat ( crDir,dataRoot,session_prefix,sessionID,
                         mmdb::io::_dir_sep );
    */

    dir = opendir(crDir);
    if (dir)  {
      rc = SDIR_Exists;
      closedir ( dir );
    } else if (mkdir(crDir,06777))
      rc = SDIR_Fail;
    else
      rc = SDIR_Created;

    return rc;

  }

  int Data::checkCrDir ( mmdb::cpstr sessionName )  {
  DIR  * dir;
  int    rc;

    makeCrDirPath ( sessionName );

   
    
    /*
    mmdb::CreateCopy   ( sessionID,sessionName );
    mmdb::CreateCopCat ( crDir,dataRoot,session_prefix,sessionID,
                         mmdb::io::_dir_sep );
    */

    dir = opendir(crDir);
    //cout<<"GDL:dir"<<"\t"<<crDir<<"\n";
    
    if (dir)  {
      rc = SDIR_Exists;
      closedir ( dir );
      readParams();
      if (!crDir)
        rc = SDIR_noResults;
    } else
      rc = SDIR_doesntExist;

    return rc;

  }

  #endif

  void Data::rmCrDir()  {
    if (crDir)  {
      removeDir ( crDir );
      delete[] crDir;
    }
    crDir = NULL;
  }

  RESULT_CODE Data::eraseSession ( mmdb::cpstr sessionName )  {

    if (confStatus!=CFG_Configured)
      return RESULT_ConfigurationError;

    makeCrDirPath ( sessionName );
    /*
    mmdb::CreateCopy   ( sessionID,sessionName );
    mmdb::CreateCopCat ( crDir,dataRoot,session_prefix,sessionID,
                         mmdb::io::_dir_sep );
    */

    if (crDir)  {
      removeDir ( crDir );
      delete[] crDir;
    }
    crDir = NULL;

    return RESULT_SessionErased;

  }


  //  ----------------------------------------------------------------

  mmdb::cpstr Data::makeFName ( mmdb::pstr & S, mmdb::cpstr FName,
                                int mod )  {
  // makes path to file FName in the session directory
  mmdb::pstr p;
    if (mod>=0)  {
      p = new char[strlen(FName)+100];
      sprintf ( p,"%s%i",FName,mod );
      if (crDir)  mmdb::CreateCopCat ( S,crDir,p );
            else  mmdb::CreateCopy   ( S,p );
      delete[] p;
    } else
      mmdb::CreateCopCat ( S,crDir,FName );
    return S;
  }

  mmdb::cpstr Data::makePDBFName ( mmdb::pstr & S,
                                   mmdb::cpstr pdbCode )  {
  mmdb::IDCode pdbid;

    strcpy ( pdbid,pdbCode );
    mmdb::LowerCase ( pdbid );

    switch (pdbFormat)  {
      default:
      case PDB_PLAIN_PDB      :
      case PDB_SPLIT_PDB      :
              mmdb::CreateCopCat ( S,"pdb",pdbid,".ent");
            break;
      case PDB_PLAIN_PDB_GZ   :
      case PDB_SPLIT_PDB_GZ   :
              mmdb::CreateCopCat ( S,"pdb",pdbid,".ent.gz");
            break;
      case PDB_PLAIN_mmCIF    :
      case PDB_SPLIT_mmCIF    :
              mmdb::CreateCopCat ( S,pdbid,".cif");
            break;
      case PDB_PLAIN_mmCIF_GZ :
      case PDB_SPLIT_mmCIF_GZ :
              mmdb::CreateCopCat ( S,pdbid,".cif.gz");
    }

    return S;

  }

  mmdb::cpstr Data::makePDBFPath ( mmdb::pstr & S,
                                   mmdb::cpstr pdbCode )  {
  mmdb::pstr   fname = NULL;
  mmdb::IDCode subdir;

    makePDBFName ( fname,pdbCode );

    switch (pdbFormat)  {
      default:
      case PDB_PLAIN_PDB      :
      case PDB_PLAIN_PDB_GZ   :
      case PDB_PLAIN_mmCIF    :
      case PDB_PLAIN_mmCIF_GZ :
              mmdb::CreateCopCat ( S,pdbDir,fname );
            break;
      case PDB_SPLIT_PDB      :
      case PDB_SPLIT_PDB_GZ   :
      case PDB_SPLIT_mmCIF    :
      case PDB_SPLIT_mmCIF_GZ :
              subdir[0] = pdbCode[1];
              subdir[1] = pdbCode[2];
              subdir[2] = char(0);
              mmdb::LowerCase ( subdir );
              mmdb::CreateCopCat ( S,pdbDir,"/",subdir,"/",fname );
    }

    return S;

  }


  int type_matrix[5][5] = {
    { -1,-1,-1,-1,-1 },
    { -1, 0, 2, 1, 3 },
    { -1, 2, 7, 5, 8 },
    { -1, 1, 5, 4, 6 },
    { -1, 3, 8, 6, 9 }
  };

  RESULT_CODE Data::queryIntfStats ( PInterface interface,
                                     int       refParamNo,
                                     mmdb::rvector  stats )  {
  int intfType;

    if (nIntfStats<=0)
      return RESULT_noInterfaceStats;

    intfType = type_matrix[interface->dclass1][interface->dclass2];
    if (intfType<0)
      return RESULT_unknownInterfaceType;

    stats[0] = interface->intArea;
    stats[1] = interface->intDeltaG;
    stats[2] = interface->stabEn;
    stats[3] = interface->PValue;
    stats[4] = interface->nHBonds;
    stats[5] = interface->nSBridges;
    stats[6] = interface->nDSBonds;

    switch (intfStats[intfType]->getStats(refParamNo,stats,&URIStats)) {
      case  0 :  return RESULT_Ok;
      case -1 :  return RESULT_wrongRefParameter;
      case -2 :  return RESULT_incompleteStatSet;
      default :  return RESULT_unknownIntfStatError;
    }

  }


  //  ----------------------------------------------------------------

  void Data::writeParams()  {
  mmdb::io::File f;
  mmdb::pstr     S;

    if (crDir)  {

      S = NULL;
      mmdb::CreateConcat ( S,crDir,params_file );
      f.assign ( S,false,true );
      if (S)  delete[] S;

      if (f.rewrite())  {
        writeParams ( f );
        f.shut();
      }

    }

  }

  void Data::writeParams ( mmdb::io::RFile f )  {
  int pdbf = pdbFormat;

    StreamWrite   ( f,query        );

    f.CreateWrite ( dataRoot       );
    f.CreateWrite ( session_prefix );
    f.CreateWrite ( phpUri         );
    f.CreateWrite ( jsrviewUri     );
    f.CreateWrite ( pdbDir         );
    f.WriteInt    ( &pdbf          );
    f.CreateWrite ( sessionID      );
    f.WriteInt    ( &sphCodeNo     );
    f.WriteReal   ( &probe_radius  );

  }


  void Data::readParams()  {
  mmdb::io::File f;
  mmdb::pstr     S;

    if (crDir)  {

      S = NULL;
      mmdb::CreateConcat ( S,crDir,params_file );
      f.assign ( S,false,true );
      if (S) delete[] S;

      if (f.reset(true))  {

        readParams ( f );
        f.shut();

      } else  {

        delete[] crDir;
        crDir = NULL;

      }

    }

  }

  void Data::readParams ( mmdb::io::RFile f )  {
  int  pdbf;

    StreamRead   ( f,query        );

    f.CreateRead ( dataRoot       );
    f.CreateRead ( session_prefix );
    f.CreateRead ( phpUri         );
    f.CreateRead ( jsrviewUri     );
    f.CreateRead ( pdbDir         );
    f.ReadInt    ( &pdbf          );
    f.CreateRead ( sessionID      );
    f.ReadInt    ( &sphCodeNo     );
    f.ReadReal   ( &probe_radius  );

    pdbFormat = (PDBDIR_FORMAT)pdbf;

  }


  //  ----------------------------------------------------------------

  int Data::writePIData()  {
  mmdb::io::File f;
  mmdb::pstr  FN;
  int   rc;

    rc = RESULT_cantWriteResults;

    if (query)  {
      FN = NULL;
      f.assign ( makeFName(FN,interface_file),false,true );
      if (f.rewrite())  {
        query->writePIData ( f );
        f.shut();
        rc = RESULT_Ok;
      }
      if (FN)  delete[] FN;
    }

    return rc;

  }

  void Data::writePIData ( mmdb::io::RFile f )  {
    query->writePIData ( f );
  }

  RESULT_CODE Data::readPIData()  {
  mmdb::io::File f;
  mmdb::pstr    FN;
  RESULT_CODE   rc;

    FN = NULL;
    f.assign ( makeFName(FN,interface_file),false,true );
    if (f.reset(true))  {
      if (!query)  query = new QueryData();
      query->readPIData ( f );
      f.shut();
      rc = RESULT_Ok;
    } else
      rc = RESULT_cantReadResults;

    if (FN)  delete[] FN;

    return rc;

  }

  void Data::readPIData ( mmdb::io::RFile f )  {
    if (!query)  query = new QueryData();
    query->readPIData ( f );
  }


  int Data::writeStructure()  {
  mmdb::io::File f;
  mmdb::pstr     FN;
  int            rc;

    rc = RESULT_cantWriteResults;

    if (query)  {
      FN = NULL;
      f.assign ( makeFName(FN,structure_file),false,true );
      if (f.rewrite())  {
        query->writeStructure ( f );
        f.shut();
        rc = RESULT_Ok;
      }
      if (FN)  delete[] FN;
    }

    return rc;

  }

  void Data::writeStructure ( mmdb::io::RFile f )  {
    query->writeStructure ( f );
  }

  RESULT_CODE Data::readStructure()  {
  mmdb::io::File f;
  mmdb::pstr    FN;
  RESULT_CODE   rc;

    FN = NULL;
    f.assign ( makeFName(FN,structure_file),false,true );
    if (f.reset(true))  {
      if (!query)  query = new QueryData();
      query->readStructure ( f );
      f.shut();
      rc = RESULT_Ok;
    } else
      rc = RESULT_cantReadResults;
    if (FN)  delete[] FN;

    return rc;

  }

  void  Data::readStructure ( mmdb::io::RFile f )  {
    if (!query)  query = new QueryData();
    query->readStructure ( f );
  }

  int Data::writeAssemblies()  {
  mmdb::io::File f;
  mmdb::pstr  FN;
  int   rc;

    rc = RESULT_cantWriteResults;

    if (query)  {
      FN = NULL;
      f.assign ( makeFName(FN,assembly_file),false,true );
      if (f.rewrite())  {
        query->writeAssemblies ( f );
        f.shut();
        rc = RESULT_Ok;
      }
      if (FN)  delete[] FN;
    }

    return rc;

  }

  void Data::writeAssemblies ( mmdb::io::RFile f )  {
    query->writeAssemblies ( f );
  }

  RESULT_CODE Data::readAssemblies()  {
  mmdb::io::File f;
  mmdb::pstr    FN;
  RESULT_CODE   rc;

    FN = NULL;
    f.assign ( makeFName(FN,assembly_file),false,true );
    if (f.reset(true))  {
      if (!query)  query = new QueryData();
      query->readAssemblies ( f );
      f.shut();
      rc = RESULT_Ok;
    } else
      rc = RESULT_cantReadResults;
    if (FN)  delete[] FN;

    return rc;

  }

  void Data::readAssemblies ( mmdb::io::RFile f )  {
    if (!query)  query = new QueryData();
    query->readAssemblies ( f );
  }

  RESULT_CODE Data::StartTextOutput ( mmdb::io::RFile f,
                                      mmdb::cpstr fileName )  {
  char S[100];
    f.assign ( fileName,true );
    if (f.rewrite())  {
      sprintf ( S,"\n PISA v.%i.%i.%i built %s "
                  " with SSM v.%i.%i.%i, SRS v.%i.%i.%i, "
                  "MMDB v.%i,%i.%i\n"
                  " Session %s",
                  MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION,
                  BuildDate,
                  ssm::MAJOR_VERSION,ssm::MINOR_VERSION,
                                     ssm::MICRO_VERSION,
                  ccp4srs::MAJOR_VERSION,ccp4srs::MINOR_VERSION,
                                         ccp4srs::MICRO_VERSION,
                  mmdb::MAJOR_VERSION,mmdb::MINOR_VERSION,
                                      mmdb::MICRO_VERSION,
                  sessionID
                );
      f.WriteLine ( S );
      f.LF();
      return RESULT_Ok;
    } else
      return RESULT_cantWriteFile;
  }


  //  ----------------------------------------------------------------

  void Data::printResultMessage ( int resCode )  {
  char S[100];

    switch (resCode)  {

      case RESULT_noResults          :
          printf ( "\n ... no results calculated.\n" );
        break;

      case RESULT_Ok                 :
  //        printf ( "\n ... normal processing.\n" );
        break;

      case RESULT_Instructions       :
          printf ( "\n ... instructions printed.\n" );
        break;

      case RESULT_Logo               :
          printf ( "\n ... version/logo printed.\n" );
        break;

      case RESULT_ConfigurationError :
          printf ( "\n"
                   " **** configuration error:\n"
                   "\n"
                   " <<< %s >>>\n\n",
                   getConfigurationMessage(ConfStatus(),S) );
        break;

      case RESULT_SessionDirFailure  :
          printf ( "\n"
                   " **** cannot create session directory:\n"
                   "\n"
                   "   %s\n\n",crDir );
        break;

      case RESULT_noSRS              :
          printf ( "\n *** no SRS library, quit\n" );
        break;

      case RESULT_noMolRefIndex      :
          printf ( "\n *** no MolRef index file, quit\n");
        break;

      case RESULT_noMolRefData       :
          printf ( "\n *** no MolRef data file, quit\n");
        break;

      case RESULT_noAgentsData       :
          printf ( "\n *** no agent.dat file in pistore, quit\n");
        break;

      case RESULT_notAnalysed        :
          printf ( "\n *** structure analysis failed, quit\n");
        break;

      case RESULT_noAtoms            :
          printf ( "\n *** no atom data found in input file, quit\n");
        break;

      case RESULT_noCoordinateData   :
          printf ( "\n *** no chains found in input file, "
                   "quit\n");
        break;

      case RESULT_notReadyForContacts :
          printf ( "\n *** not ready for contact calculations "
                   "(program error), quit\n");
        break;

      case RESULT_SurfaceCalcsFailed  :
          printf ( "\n *** surface calculations failed, quit\n" );
        break;


      case RESULT_badCoordinateFile  :
          if (query)  {
            query->MMDB->GetInputBuffer ( MMDBErrLine,MMDBErrLineNo );
            printf ( " *** problem reading the coordinate file:\n"
                     "     '%s'\n"
                     " at line number %i:\n"
                     "     %s\n",
                     query->coorFile,MMDBErrLineNo,MMDBErrLine );
          } else
            printf ( "\n *** bad coordinate file\n");
        break;

      case RESULT_SRSFailed        :
          printf ( "\n *** SRS malfunction, quit\n");
        break;

      case RESULT_noSymOps           :
          printf ( "\n *** cannot load space symmetry group data, "
                   "quit\n");
        break;

      case RESULT_noSymOp            :
          printf ( "\n *** a symmetry operation not found in "
                   "dictionaries, quit\n" );
        break;

      case RESULT_noMacromolecules   :
          printf ( "\n *** no macromolecules found in input file, "
                   "quit\n");
        break;

      case RESULT_InterfaceCalcFailure :
          printf ( "\n *** interface calculations error\n" );
        break;

      case RESULT_noRCSBData         :
          printf ( "\n +++ RCSB translation data not found\n" );
        break;

      case RESULT_noRCSBassignment   :
          printf ( "\n +++ RCSB assignment failed (program error)\n" );
        break;

      case RESULT_cantWriteResults   :
          printf ( "\n *** cannot write results into session "
                   "directory '%s'\n",sessionID );
        break;

      case RESULT_cantReadResults   :
          printf ( "\n *** cannot read results from session "
                   "directory %s\n",sessionID );
        break;

      case RESULT_noAssemblyParameters :
          printf ( "\n *** cannot load assembly analysis parameters\n" );
        break;

      case RESULT_SessionDoesntExist :
          printf ( "\n *** session '%s' does not exist.\n",
                   sessionID );
        break;

      case RESULT_noSessionResults :
          printf ( "\n *** no calculation results found in "
                   "session '%s'\n",sessionID );
        break;

      case RESULT_cantWriteFile    :
          printf ( "\n *** a file cannot be opened for writing\n" );
        break;

      case RESULT_SessionErased    :
          printf ( "\n +++ session '%s' erased\n\n",sessionID );
        break;

      case RESULT_monomerNoOutOfRange :
          printf ( "\n +++ monomer serial number out of range\n\n" );
        break;

      case RESULT_interfaceNoOutOfRange :
          printf ( "\n +++ interface serial number out of range\n\n" );
        break;

      case RESULT_assemblyNoOutOfRange :
          printf ( "\n +++ assembly serial number out of range\n\n" );
        break;

      case RESULT_cantWriteRasmolData  :
          printf ( "\n *** cannot write data for rasmol "
                   "(disk full?)\n\n" );
        break;

      case RESULT_cantWriteMGData  :
          printf ( "\n *** cannot write data for CCP4-MG "
                   "(disk full?)\n\n" );
        break;

      case RESULT_noDimerFound  :
          printf ( "\n *** no (stable) dimers found\n\n" );
        break;

      default : printf ( "\n *** unknown result code %i, quit.\n",
                         resCode );
    }

  }


  void Data::makeNoAssembliesPage ( mmdb::io::RFile f )  {

    if (query->asmStatus==ASSMB_noSymOps)  {
      f.WriteLine (
       " *** Assembly infomation is not available because\n"
       " *** crystallographic data is either absent or incomplete.\n" );
    } else  {
      f.Write (
       " *** Results are not available due to the following reason:\n"
       " ***\n"
       " ***   " );
      f.Write ( getAsmStatus(query->asmStatus) );
      f.Write (
       " ***\n" );
      if (query->asmStatus==ASSMB_Overlap)
        f.Write (
       " *** This means that crystallographic data is not correct\n"
       " *** or precise, and overlapping symmetry mates of molecules\n"
       " *** were found at crystal reconstruction. PISA does not\n"
       " *** attempt to perform an assembly analysis if there are\n"
       " *** problematic interfaces. Overlapping interfaces are\n"
       " *** marked by star in the list of interfaces.\n"
       " ***\n" );
      f.Write (
       " *** If you think that this message is issued in error,\n"
       " *** please report it to author Eugene Krissinel at\n"
       " *** " );
      f.Write ( auth_email );
      f.Write ( ", providing all details\n"
       " *** of your input.\n" );
    }

  }


  mmdb::xml::PXMLObject Data::makeNoAssembliesXML ( mmdb::cpstr sessionName )  {
  mmdb::xml::PXMLObject xml,xml1;
  mmdb::pstr            S = NULL;

    xml = new mmdb::xml::XMLObject ( xml_pisa_results );
    if (sessionName)
      xml->AddObject ( new mmdb::xml::XMLObject(xml_name,sessionName) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_status,"Ok") );

    if (query->asmStatus==ASSMB_noSymOps)  {
      xml1 = new mmdb::xml::XMLObject ( xml_status_desc,
           "Assembly infomation is not available because "
           "crystallographic data is either absent or incomplete." );
    } else  {
      xml1 = new mmdb::xml::XMLObject ( xml_status_desc,
       "Results are not available due to the following reason: " );
      xml1->AddData ( getAsmStatus(query->asmStatus) );
      if (query->asmStatus==ASSMB_Overlap)
        xml1->AddData (
          ". This means that crystallographic data is not correct "
          "or precise, and overlapping symmetry mates of molecules "
          "were found at crystal reconstruction. PISA does not "
          "attempt to perform an assembly analysis if there are "
          "problematic interfaces. Overlapping interfaces are "
          "marked by star in the list of interfaces." );
    }

    xml->AddObject ( xml1 );

    xml->AddObject ( new mmdb::xml::XMLObject ( xml_status_note,
       mmdb::CreateCopCat ( S,
         "If you think that this message is issued in error, "
         "please report it to author Eugene Krissinel at ",
         auth_email,", providing all details of your input." ) ) );

    if (S)  delete[] S;

    return xml;

  }


  void Data::makeNoInterfacesXML ( mmdb::cpstr sessionName,
                                        mmdb::cpstr fileName )  {
  mmdb::xml::PXMLObject xml;
  mmdb::pstr            S = NULL;

    xml = new mmdb::xml::XMLObject ( xml_pisa_results );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_name,sessionName) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_status,"NOk") );

    xml->AddObject ( new mmdb::xml::XMLObject ( xml_status_desc,
       "Interface information is not found in results. Make sure "
       "that analysis stage has been performed." ) );

    xml->AddObject ( new mmdb::xml::XMLObject ( xml_status_note,
       mmdb::CreateCopCat ( S,
         "If you think that this message is issued in error, "
         "please report it to author Eugene Krissinel at ",
         auth_email,", providing all details of your input." ) ) );

    xml->WriteObject ( fileName,0,2 );

    if (S)  delete[] S;

    delete xml;

  }


  // =================================================================


  int  removeDir ( mmdb::cpstr dirpath )  {
  mmdb::pstr S;
  int  rc;
    S = new char[strlen(dirpath)+100];
    strcpy ( S,"rm -rf " );
    strcat ( S,dirpath   );
    rc = system ( S );
    delete[] S;
    return rc;
  }


  mmdb::cpstr getConfigurationMessage ( int confCode, mmdb::pstr S )  {

    switch (confCode)  {

      case CFG_Configured   : strcpy ( S,"Ok" );  break;
      case CFG_Unconfigured : strcpy ( S,"Configuration not invoked"   );
                             break;
      case CFG_noConfFile   : strcpy ( S,"No configuration file"       );
                             break;
      case CFG_cantOpenFile : strcpy ( S,"Can't open configuration "
                                         "file" );
                             break;
      case CFG_noDataRoot      : strcpy ( S,"No data root"      ); break;
      case CFG_noSessionPrefix : strcpy ( S,"No session prefix" ); break;
      case CFG_noSRSDir        : strcpy ( S,"No SRS path"       ); break;
      case CFG_noMolRefDir     : strcpy ( S,"No MolRef path"    ); break;
      case CFG_noPIStoreDir    : strcpy ( S,"No PIStore path"   ); break;
      case CFG_noRasmolCommand : strcpy ( S,"No Rasmol command" ); break;
      case CFG_noJmolCommand   : strcpy ( S,"No JMol command"   ); break;
      case CFG_noCCP4MGCommand : strcpy ( S,"No CCP4MG command" ); break;

      default :  sprintf ( S,"Unknown misconfiguration code %i",
                             confCode );

    }

    return S;

  }


  mmdb::realtype round_figure ( mmdb::realtype V, int digit )  {
  mmdb::realtype f;
    f = mmdb::Pow ( 10.0,digit );
    return mmdb::mround ( V*f ) / f;
  }

}  // namespace pisa
