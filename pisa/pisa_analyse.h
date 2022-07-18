// $Id: pisa_analyse.h $
// =================================================================
//
//    20.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_analyse <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Analyse
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2013
//
// =================================================================
//

#ifndef __PISA_Analyse__
#define __PISA_Analyse__

#include "pisalib/pisa_data.h"

#undef _dimer_special

namespace pisa  {

  // =======================  Analyse  ==========================

  DefineClass(Analyse);

  class Analyse : public Data  {

    public :

      Analyse ( mmdb::cpstr confPath );
      ~Analyse();

      void setLigExclude ( mmdb::pstr ligList );
      void setLigKey     ( LIGAND_KEY ligProc );

      RESULT_CODE analyse ( mmdb::pstr sessionName, mmdb::pstr coorFile );

  #ifdef  _dimer_special
      RESULT_CODE dimerSpecial ( mmdb::pstr sessionName,
                                 mmdb::pstr pdbCode,
                                 mmdb::pstr pdbDir, int & intfNo );
  #endif

    protected :
      mmdb::pstr ligExcl;
      LIGAND_KEY ligKey;

      void InitAnalyse();

  };

}  // namespace pisa

#endif
