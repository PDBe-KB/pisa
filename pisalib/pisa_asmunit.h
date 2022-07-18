// $Id: pisa_asmunit.h $
// =================================================================
//
//    05.10.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_asmunit <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::AsmUnit
//       ~~~~~~~~~
//
//  (C) E. Krissinel, 2004-2013
//
// =================================================================
//

#ifndef __PISA_AsmUnit__
#define __PISA_AsmUnit__

#include "pisa_engunits.h"

namespace pisa  {

  //  =========================  AsmUnit  ==========================

  DefineClass(AsmUnit);

  class AsmUnit  {

    public :
      int           id;      //!< domain id (serial number in ASU)
      int           type;    //!< domain type id
      int           dclass;  //!< domain class
      mmdb::ivector intfl[_max_n_int]; //!< lists of interfaces (>=0: molecule
                       ///  No in assembly; -1: link to outside;
                       ///                  -2: not present)
      int           dissAsm;    //!< dissociated assembly no. 1,2,...
      int           symOpNo;    //!< symop (ASU) serial number
      int           rcsb_symop; //!< rcsb symop serial number
      mmdb::pstr    symOp;      //!< symmetry operation
      int           ncsOpNo;    //!< NCS opeartion serial number
      int           cell_i,cell_j,cell_k;  //!< unit cell ( -n..0..n );
      bool          fixed;      //!< true if a fixed ligand
      mmdb::mat44   T;          //!< transformation matrix
      mmdb::mat44   TF;         //!< transformation matrix in fract. coordinates
      mmdb::ChainID visualID;   //!< chain ID for visualisation and download

      AsmUnit ();
      ~AsmUnit();

      mmdb::xml::PXMLObject getAsmUnitXML ( PPDomain D, int nCellOut );

      void Copy ( PAsmUnit M );

      void write ( mmdb::io::RFile f, int nInterfaces );
      void read  ( mmdb::io::RFile f, int nInterfaces );

      void makeUnit ( RMonRef R, PPInterface interface,
                      int nInterfaces, int nInterfaces0,
                      PPDomain D, mmdb::mat44 & rom );

    protected :
      void FreeMemory();

  };

}  // namespace pisa


#endif
