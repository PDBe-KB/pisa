// $Id: pisa_query.h $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// -----------------------------------------------------------------
//
//  **** Module  :  PISA_Query <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::QueryData
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2014
//
// =================================================================
//

#ifndef __PISA_Query__
#define __PISA_Query__


#include "pisa_domains.h"
#include "pisa_interface.h"
#include "pisa_assemblies.h"
#include "pisa_types.h"

#include "mmdb2/mmdb_utils.h"
#include "ccp4srs/ccp4srs_manager.h"


namespace pisa  {

  //  =========================  QueryData  =========================

  DefineClass(QueryData);

  class QueryData : public mmdb::io::Stream  {

    public :
      DATASTATUS       dataStatus;
      mmdb::pstr       coorFile;
      mmdb::ERROR_CODE readErr;

      mmdb::realtype resolution;
      mmdb::realtype cell_a,cell_b,cell_c; //!< cell parameters
      mmdb::realtype cell_alpha,cell_beta,cell_gamma;
      int            OrthCode;           //!< and orthogonalization code
      mmdb::SymGroup spaceGroup;
      int            nSymOps;
      int            nNCSOps;
      mmdb::pstr     title;

      mmdb::PManager MMDB;
      PDomains       domains;
      PInterfaces    PI;
      PAssemblies    A;             //!< crystal assemblies
      PAssembly      Complex;       //!< complex given by coordinate part
      ASSMB_RC       asmStatus,cmpStatus;
      mmdb::pstr     warnUNKRes;    //!< warinings: unknown residues
      mmdb::pstr     warnExcluded;  //!< warinings: excluded residues
      mmdb::pstr     warnOverlap;   //!< warinings: interface overlap
      mmdb::pstr     warnAssembly;  //!< warinings from assembly building

      QueryData ();
      QueryData ( mmdb::io::RPStream Object );
      ~QueryData();

      void  FreeMemory     ();
      void  DeleteWarnings ();

      void setLigExclude ( mmdb::cpstr ligList );
      void setLigKey     ( LIGAND_KEY  ligProc );
      void setAsIsKey    ( AS_IS_KEY asisProc );

    
      DATASTATUS analyseStructure ( mmdb::cpstr  fileName,
                                    PMolRefIndex molRef,
                                    mmdb::cpstr  agentsRef );
      DATASTATUS analyseStructure ( mmdb::io::RFile f,
                                    mmdb::cpstr     fileSpec,
                                    PMolRefIndex    molRef,
                                    mmdb::cpstr     agentsRef );

      mmdb::pstr getDomainName     ( mmdb::pstr & S, int domainNo );
      mmdb::pstr getDomainRange    ( mmdb::pstr   S, int domainNo,
                                     int fieldLen );
      mmdb::pstr getChainID        ( mmdb::pstr   S, int domainNo );
      mmdb::pstr getDomainCode     ( mmdb::pstr & S, int domainNo );
      int        getDomainClass    ( int domainNo );
      mmdb::pstr getInterfaceName  ( mmdb::pstr & S, int interfaceNo );
      mmdb::pstr getInterfaceName1 ( mmdb::pstr & S, int interfaceNo );

      void AlignSimilarDomains ();
      void SortDomains         ();
      void makeSSGraphs        ();
      void AlignAllDomains     ( int startNo, int endNo );
      PROSURF_RC  CalcSurfaces ( PProSurf     proSurf,
                                 PMolRefIndex molRef );
      RESULT_CODE CalcContacts ( int nCellLayers );

      void calcDimensions      ();
      void checkResidueBSA     ();
      void SortInterfaces      ();

      void makeOverlapWarnings();
      bool         areWarnings();

      void setLigandsAssemble ( bool On ); //!< exclude on/off
      void excludeLigands     (); //!< uses ligExcl

      ASSMB_RC    calcAssemblies   ();
      RESULT_CODE assignRCSBSymOps ( mmdb::pstr rcsbSymOpFile );
      RESULT_CODE assignRCSBSymOps ( PPRCSBData rcsbData, int nRCSBData );

      int  getSimilarInterface (
        PDomain DA, mmdb::mat44 & Ta, int domain1, mmdb::realtype rmsdA,
        PDomain DB, mmdb::mat44 & Tb, int domain2, mmdb::realtype rmsdB,
        mmdb::mat44 & Tm );


      int         getNofInterfaces ();
      int         getNofITypes     ();
      int         getNofDomains    ();
      void        getNofDomains    ( int & nAADomains, int & nNADomains,
                                     int & nLigands );
      int         getNofNCSParents ();
      int         getNofAssemblies ();
      PInterface  getInterface     ( int interfaceNo );
      PPInterface getInterfaces    ();
      mmdb::pstr  getStructName    ( mmdb::pstr & S );
      mmdb::cpstr getStructName    ();
      PDomain     getDomain        ( int domainNo );  // 0..nDomains-1
      PPDomain    getDomains       ();
      int         getNCSParent     ( int domainNo );

      void  calcDomainStats ( int   domainID,
                              mmdb::realtype & bsa, // buried asa in interfaces
                              mmdb::realtype & DeltaG,
                              int & nHB,
                              int & nSB,
                              int & nDS,
                              int & nAt, // no. of buried atoms
                              int & nRes // no. of buried residues
                            );

      void getExclLigandList ( mmdb::PResName & ligName,
                               int & nLigNames );

      void readPIData       ( mmdb::io::RFile f );
      void readStructure    ( mmdb::io::RFile f );
      void readAssemblies   ( mmdb::io::RFile f );

      void writePIData      ( mmdb::io::RFile f );
      void writeStructure   ( mmdb::io::RFile f );
      void writeAssemblies  ( mmdb::io::RFile f );

      void outputSIData     ( mmdb::io::RFile f, int nCellLayers );
      void outputAsmData    ( mmdb::io::RFile f );
      void outputAssemblies ( mmdb::pstr outDir, int maxNofAsmSets );

      mmdb::xml::PXMLObject getInterfacesXML(int asis_param);

      void read  ( mmdb::io::RFile f );
      void write ( mmdb::io::RFile f );

    protected :
      mmdb::pstr  ligExcl;  //!< list of ligands to exclude from analysis
      LIGAND_KEY  ligKey;   //!< ligand processing mode

      AS_IS_KEY   asisKey;
      int         nInterfaces,nDomains,nNCSParents,nAssemblies;

      void  InitQueryData();

      ///  Reads external coordinate file
      int  getStructure ( mmdb::io::RFile f, mmdb::cpstr fileSpec );

  };

  DefineStreamFunctions(QueryData);

  extern void SetSyminfoLib    ( mmdb::pstr syminfo_lib );
  extern void SetMMDBReadFlags ( mmdb::PManager MMDB    );
  extern bool CheckDataStatus  ( DATASTATUS ds, RESULT_CODE & rc,
                                 mmdb::pstr & msg );

}  // namespace pisa


#endif
