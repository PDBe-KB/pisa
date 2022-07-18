// $Id: pisa_view.h $
// =================================================================
//
//    20.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_view <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::View
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2013
//
// =================================================================
//

#ifndef __PISA_View__
#define __PISA_View__

#include "pisa_data.h"

namespace pisa  {

  // =========================  View  ===========================

  DefineClass(View)

  class View : public Data  {

    public :

      View ( mmdb::cpstr confPath );
      ~View();

      RESULT_CODE ViewInput     ( mmdb::cpstr   sessionName,
                                  VIEWER_TARGET target,
                                  mmdb::cpstr   fileName );
      RESULT_CODE ViewMonomer   ( mmdb::cpstr   sessionName,
                                  int           serialNo,
                                  VIEWER_TARGET target,
                                  mmdb::cpstr   fileName );
      RESULT_CODE ViewInterface ( mmdb::cpstr   sessionName,
                                  int           serialNo,
                                  VIEWER_TARGET target,
                                  mmdb::cpstr   fileName );
      RESULT_CODE ViewAssembly  ( mmdb::cpstr   sessionName,
                                  int           serialNo,
                                  bool          dissociated,
                                  VIEWER_TARGET target,
                                  mmdb::cpstr   fileName );

    protected :

      void InitPISAView ();
      int  launchRasmol ( mmdb::pstr FN );
      int  launchMG     ( mmdb::pstr FN );

  };

}  // namespace pisa


#endif
