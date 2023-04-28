// $Id: pisa_detail.h $
// =================================================================
//
//    04.03.16   <--  Date of Last Modification.
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
//  (C) E. Krissinel 2007-2016
//
// =================================================================
//

#ifndef __PISA_Detail__
#define __PISA_Detail__

#include "pisa_lists.h"

namespace pisa  {

  // ========================  Detail  ==========================

  DefineClass(Detail);

  class Detail : public Lists  {

    public :

      Detail ( mmdb::cpstr confPath );
      ~Detail();

      RESULT_CODE DetailInterface   ( mmdb::cpstr sessionName,
                                      mmdb::cpstr fileName,
                                      int  serialNo );
      RESULT_CODE DetailMonomer     ( mmdb::cpstr sessionName,
                                      mmdb::cpstr fileName,
                                      int  serialNo );
      RESULT_CODE DetailAssembly    ( mmdb::cpstr sessionName,
                                      mmdb::cpstr fileName,
                                      int  serialNo );
      RESULT_CODE Remark350         ( mmdb::cpstr sessionName,
                                      mmdb::cpstr fileName,
                                      int  serialNo );
      RESULT_CODE MakeInterfacesXML ( mmdb::cpstr sessionName,
                                      mmdb::cpstr fileName );
      RESULT_CODE MakeAssembliesXML ( mmdb::cpstr sessionName,
                                      mmdb::cpstr fileName );
      void   set_As_Is_Key_xml ( AS_IS_KEY as_is_Proc );

    protected :
    
      AS_IS_KEY   asisKey;

      void   InitDetail();

      void   printInterfaceSummary  ( mmdb::io::RFile f,
                                      PInterface interface );
      mmdb::PAtom GetAtom           ( int serNum );
      void   printBondTable         ( mmdb::io::RFile f,
                                      mmdb::cpstr     Title,
                                      mmdb::cpstr     noBondMsg,
                                      PInterface      interface,
                                      PIBond          Bond,
                                      int             nBonds );
      void   printIntResTable       ( mmdb::io::RFile f,
                                      mmdb::cpstr     Title,
                                      mmdb::cpstr     noResMsg,
                                      PInterface      interface,
                                      int             structNo );

      void   printMonomerSummary    ( mmdb::io::RFile f, int domainNo );
      void   printResidueData       ( mmdb::io::RFile f, int domainNo );

      void   printAssemblySummary   ( mmdb::io::RFile f, PAssembly A );
      void   printEngagedInterfaces ( mmdb::io::RFile f, PAssembly A );
      void   printAssembledMonomers ( mmdb::io::RFile f, PAssembly A );
      void   printTransformations   ( mmdb::io::RFile f, PAssembly A );

      void   printRemark350         ( mmdb::io::RFile f, PAssembly A,
                                      bool fractional );
      void   printRemark350M        ( mmdb::io::RFile f );
      void   printRemark350         ( mmdb::io::RFile f, PDomain D,
                                      int  domainNo );
    
  };

}  // namespace pisa

#endif

