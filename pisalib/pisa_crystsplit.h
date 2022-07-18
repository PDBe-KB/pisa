// $Id: pisa_crystsplit.h $
// =================================================================
//
//    06.12.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_crystsplit <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::CrystSplit
//       ~~~~~~~~~
//
//  (C) E. Krissinel, 2004-2013
//
// =================================================================
//

#ifndef __PISA_CrystSplit__
#define __PISA_CrystSplit__

#include "pisa_assembly.h"

namespace pisa  {

  //  ========================  CrystSplit  =========================

  DefineClass(CrystSplit);

  class CrystSplit  {

    public :
      mmdb::ivector intf;  //!< vector of engaged interfaces
      PPAssembly       A;  //!< list of assemblies obtained
      int    nAssemblies;  //!< number of assemblies obtained
      int          Score;  //!< score 0-7
      bool     equiv_all;  //!< True if all assemblies are equivalent
      bool    stable_all;  //!< Trye if all assemblies are stable
      bool   orig_chains;  //!< True if all ASU chains are used

      CrystSplit ();
      ~CrystSplit();

      void calcScore    ();
      int  getMaxAsmSize();

      void Copy ( PMultimerSet multSet, PPInterface interface,
                  int nInterfaces, int nInterfaces0,
                  PPDomain D, mmdb::mat44 & rom );

      void checkOriginalOrientations ( PDomains Domains,
                                       mmdb::ivector icnt );

      mmdb::xml::PXMLObject getCrystSplitXML ( int     serNo,
                                               PDomains    D,
                                               PInterfaces PI,
                                               int   nCellOut );

      void write ( mmdb::io::RFile f, int nInterfaces );
      void read  ( mmdb::io::RFile f, int nInterfaces );

    protected :
      void FreeMemory();

  };


}  // namespace pisa


#endif
