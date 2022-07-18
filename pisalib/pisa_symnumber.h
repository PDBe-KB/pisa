// $Id: pisa_symnumber.h $
// =================================================================
//
//    20.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_symnumber <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::SymNumber
//       ~~~~~~~~~
//
//  (C) E. Krissinel, 2007-2013
//
// =================================================================
//

#ifndef __PISA_SymNumber__
#define __PISA_SymNumber__

#include "pisa_domain.h"

namespace pisa  {

  // =========================  SymNumber  ==========================

  // forward definitions
  DefineClass(Multimer);
  DefineClass(Assembly);

  DefineClass(SymNumber);

  class SymNumber {

    public :

      SymNumber ();
      ~SymNumber();

      int getSymNumber ( PMultimer U, int n, PDomains D,
                         mmdb::mat44 * TMatrix );
      int getSymNumber ( PAssembly Asm, PDomains D );

      void getStorage ( mmdb::rmatrix & S1, mmdb::rmatrix & S2,
                        mmdb::rvector & S3, mmdb::rvector & S4 );

    protected :
      mmdb::imatrix  C;
      PTFrame        F;
      mmdb::rmatrix  A,U,V;
      mmdb::rvector  W,RV1;
      int            nUnits;

      void FreeMemory      ();
      int  superposeFrames ( mmdb::mat44 & T,
                             mmdb::ivector f1, mmdb::ivector f2,
                             int nFrames );
      int  calcSymNumber   ();

    private :
      int allocSize;

  };

}  // namespace pisa

#endif
