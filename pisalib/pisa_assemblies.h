// $Id: pisa_assemblies.h $
// =================================================================
//
//    06.12.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_assemblies <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Assemblies
//       ~~~~~~~~~
//
//  (C) E. Krissinel, 2004-2013
//
// =================================================================
//

#ifndef __PISA_Assemblies__
#define __PISA_Assemblies__

#include "pisa_crystsplit.h"
#include "pisa_asmstock.h"


namespace pisa  {


  //  ========================  Assemblies  =========================

  DefineClass(Assemblies);

  class Assemblies : public mmdb::io::Stream  {

    public :
      PPCrystSplit crystSplit; //!< assembly sets obtained
      int        nCrystSplits; //!< total number of assembly sets
      int              nCSRes; //!< number of splits scored less than 8
      PAsmStock      asmStock; //!< assembly stock (bulk solution)
      int         nInterfaces; //!< number of interfaces
      bool        rcsb_symops; //!< True if rcsb symops have been assigned
      bool        orig_chains; //!< True if all ASU chains are used

      Assemblies ();
      Assemblies ( mmdb::io::RPStream Object );
      ~Assemblies();

      void Sort  ( PAssembly Complex ); //!< Complex may be NULL

      int  getMaxAsmSize ();
      int  getMaxAsmSize ( int CrystSplitNo );

      int  getNofAssembliesInSplits();

      void calcInterfaceScores ( PPInterface interface );

      void makeChainMapping    ();
      void makeOrientations    ( PDomains D, mmdb::PManager MMDB );
      void checkOriginalOrientations ( PDomains D    );
      void Orth2Frac           ( mmdb::PManager MMDB );
      void assignRCSBSymOps    ( PRCSBData  rcsbData );

      void TrimNofSets         ( mmdb::pstr & warning );

      void calcAsmStock        ( PDomains D, mmdb::PManager MMDB );


      mmdb::xml::PXMLObject getAssembliesXML ( mmdb::cpstr name,
                                               PDomains       D,
                                               PInterfaces   PI,
                                               int     nCellOut );

      void write ( mmdb::io::RFile f );
      void read  ( mmdb::io::RFile f );

    protected :
      void InitAssemblies();
      void FreeMemory    ();

  };

  DefineStreamFunctions(Assemblies);


  // =================================================================

  extern mmdb::cpstr getAsmStatus  ( ASSMB_RC asmRC  );
  extern void SetMaxNofCrystSplits ( int maxNCrystSplits );

}  // namespace pisa


#endif
