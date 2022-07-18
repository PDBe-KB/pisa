// $Id: pisa_rcsbdat.h $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// -----------------------------------------------------------------
//
//  **** Module  :  PISA_RCSBDat <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::RCSBData
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2014
//
// =================================================================
//

#ifndef __PISA_RCSBDat__
#define __PISA_RCSBDat__

#include "mmdb2/mmdb_io_stream.h"
#include "mmdb2/mmdb_defs.h"

namespace pisa  {

  //  ==========================  RCSBData  =========================

  DefineClass(RCSBData)
  DefineStreamFunctions(RCSBData)

  class RCSBData : public mmdb::io::Stream  {

    public:

      mmdb::SymGroup spaceGroup;   //!< space symmetry group
      int            nSymOps;      //!< number of symops
      mmdb::ivector  rcsb_symop;   //!< rcsb_symop[i] is rcsb's symop
                                   /// number for ith symop in CCP4

      RCSBData ();
      RCSBData ( mmdb::cpstr spGroup, int nSO );
      RCSBData ( mmdb::io::RPStream Object );

      ~RCSBData();

      void write ( mmdb::io::RFile f );
      void read  ( mmdb::io::RFile f );

    protected :

      void InitRCSBData ( mmdb::cpstr spGroup, int nSO );
      void FreeMemory   ();

  };

  extern int getRCSBData ( mmdb::pstr fileName, mmdb::pstr spaceGroup,
                           RPRCSBData rcsbData );

  extern int getRCSBData ( mmdb::pstr fileName, PPRCSBData & rcsbData,
                           int & nRCSBData );

  extern void deleteRCSBData ( PPRCSBData & rcsbData, int nRCSBData );


}  // namespace pisa


#endif
