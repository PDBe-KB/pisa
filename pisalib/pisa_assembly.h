// $Id: pisa_assembly.h $
// =================================================================
//
//    24.03.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_asmsets <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Assembly
//       ~~~~~~~~~
//
//  (C) E. Krissinel, 2004-2014
//
// =================================================================
//

#ifndef __PISA_Assembly__
#define __PISA_Assembly__

#include "pisa_asmunit.h"

namespace pisa  {


  //  =========================  Assembly  ==========================

  DefineClass(Assembly);

  class Assembly  {

    public :
      PPAsmUnit       M;         //!< list of AsmUnits (their ids) in assembly
      int             serNo;       //!< assembly serial number
      int             asmSize;     //!< assembly size (number of AsmUnits)
      int             mmSize;      //!< macromolecular size (number of macromols)
      int             freeSize;    //!< number of free (not fixed) monomers
      mmdb::realtype  dissEn;      //!< maximal enthalpy of dissociation
      mmdb::realtype  entropy;     //!< entropy change at dissociation
      mmdb::realtype  freeEn;      //!< maximal free energy of dissociation
      mmdb::realtype  dissIntArea; //!< dissociation interface area
      mmdb::realtype  dissEn0;     //!< ground-level enthalpy of dissociation
      mmdb::realtype  entropy0;    //!< ground-level entropy change at dissociation
      mmdb::realtype  freeEn0;     //!< ground-level free energy of dissociation
      mmdb::realtype  asa;         //!< accessible surface area
      mmdb::realtype  bsa;         //!< buried surface area
      mmdb::realtype  seGain;      //!< solvation energy gain
      int             nUC;         //!< number of this assemblies in unit cells
      int             type;        //!< assembly type id in assembly set
      int             nDiss;       //!< number of dissociating parts
      int             symNumber;   //!< assembly symmetry number
      bool            stable;      //!< True if assembly is stable

      Assembly ();
      ~Assembly();

      mmdb::pstr getFormula     ( mmdb::pstr & S, PDomains D, int maxLen,
                                  mmdb::cpstr leadStr, bool html=false );
      mmdb::pstr getComposition ( mmdb::pstr & S, PDomains D, int maxLen,
                                  mmdb::cpstr leadStr, bool html=false );
      mmdb::pstr getDissPattern ( mmdb::pstr & S, PDomains D, int maxLen,
                                  mmdb::cpstr leadStr, bool html=false );
      void getEngagedInterfaces ( mmdb::ivector intf, int nInterfaces );
      void getDissInterfaces    ( mmdb::ivector intf, int nInterfaces );

      mmdb::pstr getFileSafeName ( mmdb::pstr & S, PDomains D );

      mmdb::PManager getAssemblyStructure ( mmdb::PManager MMDB,
                                            PDomains   D );

      void calcStats ( int              domainID, // domainID
                       int &            nOcc,     // occurence number
                       mmdb::realtype & dom_bsa,
                       mmdb::realtype & dom_DeltaG,
                       int &            dom_nHB,
                       int &            dom_nSB,
                       int &            dom_nDS,
                       int &            dom_nAt,  // no. atoms in interfaces
                       int &            dom_nRes, // no. res-s in interfaces
                       PPInterface      interface,
                       PPDomain         D,
                       int              nInterfaces );

      /// Returns type of first macromolecular monomeric unit, or first
      /// ligand unit if mmSize==0.
      int getFirstMonType();

      /// Returns the occurence number of specified monomer type in the
      /// assembly.
      int getMonTypeOccurence ( int monType );

      void Copy  ( PMultimer U, PPInterface interface,
                   int nInterfaces, int nInterfaces0,
                   PPDomain D, mmdb::mat44 & rom );

      /// Convenience function to average assembly data (not including
      /// the array of monomers). nSample=0 initializes this assembply
      /// with data from A, future calls make the averaging and
      /// autpincrement of nSample.
      void averageData ( PAssembly A, int & nSample );

      ///  makeChainMapping() fills AsmUnit::visualID fields, which
      ///  are used as chain IDs when assembly is visualised or
      ///  downloaded
      void makeChainMapping ();

      ///  Assemblies that are coming from the assembling engine
      /// (CAssembler), may be placed in any unit cell and in
      /// any part of a unit cell. makeOrientation() brings the assembly
      /// into original ASU and unit cell, if that is necessary and
      /// possible.
      void makeOrientation ( PDomains D, mmdb::PManager MMDB );

      void countOriginalOrientations ( PPDomain D, mmdb::ivector icnt );

      void Orth2Frac        ( mmdb::PManager MMDB );

      void assignRCSBSymOps ( PRCSBData rcsbData );
      
      mmdb::pstr sizeName ( mmdb::pstr S );  // returns S filled

      mmdb::xml::PXMLObject getasisAssemblyXML ( PDomains    D,
                                             PInterfaces PI,
                                             int   nCellOut);
      
      mmdb::xml::PXMLObject getAssemblyXML ( PDomains    D,
                                             PInterfaces PI,
                                             int   nCellOut,
                                             int      Score);

      void write ( mmdb::io::RFile f, int nInterfaces );
      void read  ( mmdb::io::RFile f, int nInterfaces );

    protected :
      void FreeMemory  ();

  };


}  // namespace pisa


#endif
