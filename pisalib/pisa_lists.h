// $Id: pisa_lists.h $
// =================================================================
//
//    20.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_lists <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Lists
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2013
//
// =================================================================
//

#ifndef __PISA_Lists__
#define __PISA_Lists__

#include "pisa_data.h"

namespace pisa  {

  // ========================  Lists  ===========================

  DefineClass(Lists);

  class Lists : public Data  {

    public :

      Lists ( mmdb::cpstr confPath );
      ~Lists();

      RESULT_CODE ListInterfaces ( mmdb::cpstr sessionName,
                                   mmdb::cpstr fileName );
      RESULT_CODE ListMonomers   ( mmdb::cpstr sessionName,
                                   mmdb::cpstr fileName );
      RESULT_CODE ListAssemblies ( mmdb::cpstr sessionName,
                                   mmdb::cpstr fileName );
      RESULT_CODE ListStock      ( mmdb::cpstr sessionName,
                                   mmdb::cpstr fileName );

      RESULT_CODE ListInterfaces ( mmdb::io::RFile f );
      RESULT_CODE ListMonomers   ( mmdb::io::RFile f );
      RESULT_CODE ListAssemblies ( mmdb::io::RFile f );
      RESULT_CODE ListStock      ( mmdb::io::RFile f );


      RESULT_CODE ListInterfaces_csv ( mmdb::io::RFile f );
      RESULT_CODE ListMonomers_csv   ( mmdb::io::RFile f );
      RESULT_CODE ListAssemblies_csv ( mmdb::io::RFile f );


    protected :
      void InitLists();

      void makePQSTable          ( mmdb::io::RFile f,
                                   int Score1, int Score2 );
      void makeComplexAsIsPage   ( mmdb::io::RFile f );
      void makeAssembliesPage    ( mmdb::io::RFile f );
      void makeStockPage         ( mmdb::io::RFile f );
      void writePQSTableAssembly ( mmdb::io::RFile f, PAssembly A,
                                   int setNo, int serNo,
                                   mmdb::pstr & F );

      void makePQSTable_csv          ( mmdb::io::RFile f,
                                       int Score1, int Score2 );
      void makeComplexAsIsPage_csv   ( mmdb::io::RFile f );
      void makeAssembliesPage_csv    ( mmdb::io::RFile f );
      void writePQSTableAssembly_csv ( mmdb::io::RFile f, PAssembly A,
                                       int setNo, int serNo,
                                       mmdb::pstr & F );

  };

}  // namespace pisa

#endif
