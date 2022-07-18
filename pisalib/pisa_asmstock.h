// $Id: pisa_asmstock.h $
// =================================================================
//
//    19.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_asmstock <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::AsmStock
//       ~~~~~~~~~
//
//  (C) E. Krissinel, 2004-2013
//
// =================================================================
//

#ifndef __PISA_AsmStock__
#define __PISA_AsmStock__

#include "pisa_assembly.h"

namespace pisa  {


  //  =========================  AsmStock  ==========================

  DefineClass(AsmStock);

  class AsmStock  {

    public :

      enum RETURN_CODE { Ok,UnknownType,WrongProductSize };

      PPAssembly          A;  //!< list of assemblies in stock
      int       nAssemblies;  //!< number of assemblies in stock
      mmdb::rvector asuConc;  //!< asu concentration profile
      mmdb::rmatrix asmConc;  //!< assembly concentration profiles
      int             nConc;  //!< concentration profile length

      AsmStock ( int n );
      ~AsmStock();

      RETURN_CODE addToStock   ( PAssembly       assembly,
                                 mmdb::ivector & count );
      RETURN_CODE processStock ( PDomains        D,
                                 mmdb::realtype  minConc,
                                 mmdb::realtype  maxConc );

      mmdb::realtype getRelativeConc     ( int n );
      mmdb::realtype getAggregationIndex ( int asmNo, int n );

      void write ( mmdb::io::RFile f, int nInterfaces );
      void read  ( mmdb::io::RFile f, int nInterfaces );

    protected :
      int  asuSize;  //!< number of macromolecules in ASU
      void InitAsmStock ( int n );
      void FreeMemory   ();

    private:
      int  nAsmAlloc;

  };


}  // namespace pisa


#endif
