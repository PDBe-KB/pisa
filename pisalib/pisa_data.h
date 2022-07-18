// $Id: pisa_data.h $
// =================================================================
//
//    04.03.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_data <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Data
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2014
//
// =================================================================
//

#ifndef __PISA_Data__
#define __PISA_Data__

#include "pisa_query.h"
#include "pisa_datastats.h"
#include "pisa_types.h"
#include "ssm/ssm_graph.h"
#include "ccp4srs/ccp4srs_manager.h"
#include "ccp4srs/ccp4srs_defs.h"

namespace pisa  {

  // ========================  Data  ============================

  DefineClass(Data);

  class Data  {

    public :
      Data ( mmdb::cpstr confPath );
      virtual ~Data();

      void        printResultMessage ( int resCode             );
      RESULT_CODE eraseSession       ( mmdb::cpstr sessionName );

      inline mmdb::pstr get_data_root  ()  { return dataRoot;   }
      inline mmdb::pstr get_rasmol_path()  { return rasmol_com; }
      inline mmdb::pstr get_jmol_path  ()  { return jmol_com;   }
      inline mmdb::pstr get_ccp4mg_path()  { return ccp4mg_com; }

      void set_ccp4mg_path ( mmdb::cpstr mg_path );

      mmdb::pstr getConfigurationTemplate ( mmdb::pstr & S );

      RESULT_CODE queryIntfStats ( PInterface interface,
                                   int       refParamNo,
                                   mmdb::rvector  stats );

      void setSessionPrefix ( mmdb::cpstr prefix );

    protected :

      ccp4srs::PManager SRS;  //!< common SRS resource
      PMolRefIndex   molRef;  //!< common molecular reference
      mmdb::pstr  agentsRef;  //!< crystallization agents reference
      PProSurf      proSurf;  //!< surface calculation module
      StatDistSet  URIStats;  //!< unreferenced interface stats
      PPDataStats intfStats;  //!< interface statistics for radar plots
      int        nIntfStats;  //!< number of interface parameters

      ///  Program configuration
      CONFIGURATION_STATUS confStatus; //!< configuration status
      mmdb::pstr  cfgPath;        //!< path to the Configuration File
      mmdb::pstr  syminfo_lib;    //!< path to CCP4 symmetry data file
      mmdb::pstr  dataRoot;       //!< general data root
      mmdb::pstr  srsDir;         //!< path to SRS directory
      mmdb::pstr  MolRefDir;      //!< path to MolRef directory
      mmdb::pstr  PIStoreDir;     //!< path to PIStore directory
      mmdb::pstr  helpDir;        //!< path to Help directory
      mmdb::pstr  rasmol_com;     //!< rasmol launch command
      mmdb::pstr  jmol_com;       //!< jmol launch command
      mmdb::pstr  ccp4mg_com;     //!< CCP4-MG launch command

      ///  Session configuration
      mmdb::pstr  session_prefix; //!< prefix for session directories
      mmdb::pstr  sessionID;      //!< a unique session ID
      mmdb::pstr  crDir;          //!< session directory

      ///  Server configuration
      mmdb::pstr    phpUri;       //!< URI to directory with php scripts
      mmdb::pstr    helpUri;      //!< URI to directory with help files
      mmdb::pstr    dnlUrl;       //!< Url to download files
      mmdb::pstr    jsrviewUri;   //!< URI to jsrview directory
      mmdb::pstr    pdbDir;       //!< path to PDB directory
      PDBDIR_FORMAT pdbFormat;    //!< file format in PDB directory
      int           expiry_time;  //!< time for cleaning the storage secs
      int           erase_time;   //!< time for unconditional cleaning secs

      ///  Query structure data
      PQueryData query;
      int        nCellLayers; //!< number of cell layers to check
      int        nCellOut;    //!< num. of cell layers to output in symops

      char  MMDBErrLine[300];  //!< filled by MMDB-reading routines
      int   MMDBErrLineNo;

      ///  general result code (signal)
//      RESULT_CODE resultCode;

      ///  Interface calculation parameters
      int            sphCodeNo;     //!< spherical code number
      mmdb::realtype probe_radius;  //!< radius of probe sphere

      void  FreeMemory   ();
      void  InitData     ();
      void  FreeIntfStats();

      bool  getConfPath  ( mmdb::io::RFile f, mmdb::pstr S, int lenS,
                           mmdb::cpstr   tag, mmdb::pstr & path,
                           bool dir );
      bool  getConfLine  ( mmdb::io::RFile f, mmdb::pstr S, int lenS,
                           mmdb::cpstr   tag, mmdb::pstr & line );
      bool  getConfInt   ( mmdb::io::RFile f, mmdb::pstr S, int lenS,
                           mmdb::cpstr   tag, int & iv );
      bool  getConfBool  ( mmdb::io::RFile f, mmdb::pstr S, int lenS,
                           mmdb::cpstr   tag, bool & B );
      bool  getConfOnOff ( mmdb::io::RFile f, mmdb::pstr S, int lenS,
                           mmdb::cpstr   tag, bool & B );
      void  readConfiguration ( mmdb::cpstr confPath );
      CONFIGURATION_STATUS ConfStatus   ();

      DICT_RC openMolRef       ( bool indexOnly );
      DICT_RC openAgentsRef    ();
      DICT_RC openDictionaries ( bool read_intf_stats );
      ccp4srs::CCP4SRS_RC getEnergyTypes   ( mmdb::PManager MMDB );
      void    makeProSurf      ();
      DICT_RC readIntfStats    ();

      RESULT_CODE readAssemblyParameters  ();

      void    makeCrDirPath  ( mmdb::cpstr sessionName );
      int     makeCrDir      ( mmdb::cpstr sessionName );
      int     checkCrDir     ( mmdb::cpstr sessionName );
      void    rmCrDir        ();

      mmdb::cpstr  makeFName ( mmdb::pstr & S, mmdb::cpstr FName,
                               int mod=-1 );

      mmdb::cpstr  makePDBFName ( mmdb::pstr & S, mmdb::cpstr pdbCode );
      mmdb::cpstr  makePDBFPath ( mmdb::pstr & S, mmdb::cpstr pdbCode );

      void    writeParams    ();
      void    readParams     ();
      virtual void  writeParams ( mmdb::io::RFile f );
      virtual void  readParams  ( mmdb::io::RFile f );

      int     writePIData    ();
      int     writeStructure ();
      int     writeAssemblies();
      void    writePIData    ( mmdb::io::RFile f );
      void    writeStructure ( mmdb::io::RFile f );
      void    writeAssemblies( mmdb::io::RFile f );

      RESULT_CODE readPIData     ();
      RESULT_CODE readStructure  ();
      RESULT_CODE readAssemblies ();
      void    readPIData     ( mmdb::io::RFile f );
      void    readStructure  ( mmdb::io::RFile f );
      void    readAssemblies ( mmdb::io::RFile f );

      RESULT_CODE StartTextOutput  ( mmdb::io::RFile f,
                                     mmdb::cpstr fileName );
      void    makeNoAssembliesPage ( mmdb::io::RFile f );
      mmdb::xml::PXMLObject makeNoAssembliesXML ( mmdb::cpstr sessionName );
      void    makeNoInterfacesXML  ( mmdb::cpstr sessionName,
                                     mmdb::cpstr fileName );

  };

  extern int            removeDir    ( mmdb::cpstr dirpath );
  extern mmdb::cpstr    getConfigurationMessage ( int confCode,
                                                  mmdb::pstr S   );
  extern mmdb::realtype round_figure ( mmdb::realtype V, int digit );

}  // namespace pisa

#endif

