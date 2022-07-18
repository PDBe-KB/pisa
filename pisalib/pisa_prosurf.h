// $Id: pisa_prosurf.h $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  ProSurf <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::ProSurf
//       ~~~~~~~~~
//
//  (C) E. Krissinel, 2003-2014
//
// =================================================================
//

#ifndef  __PISA_ProSurf__
#define  __PISA_ProSurf__

#include "mmdb2/mmdb_manager.h"

#include "pisa_molref.h"

namespace pisa  {

  // =========================  ProSurf  ============================

  enum PROSURF_RC  {
    PROSURF_Stop2     = -2,
    PROSURF_Stop1     = -1,
    PROSURF_Ok        = 0,
    PROSURF_noAtoms   = 1,
    PROSURF_noManager = 2
  };

  enum SURFSEL_KEY  {
    SURFF_Select     = 0,
    SURFF_dontSelect = 1,
    SURFF_Selected   = 2
  };

  DefineClass(ProSurf);

  class ProSurf : public mmdb::io::Stream  {

    public :
      // Public fields:
      //   Selection handles for surface atoms
      int  selHndSurf1;    // 1st or only structure
      int  selHndSurf2;    // 2nd structure
      //   Selection handles for interface atoms
      int  selHndInt1;     // 1st structure
      int  selHndInt2;     // 2nd structure
      //   Solvent-Accessible Surface areas
      mmdb::realtype surfArea1;  // SAS area of 1st or only structure
      mmdb::realtype surfArea2;  // SAS area 2nd structure
      mmdb::realtype surfArea12; // SAS area of both structures in contact
      //   Interface areas
      mmdb::realtype intArea1;   // area of the interface of 1st structure
      mmdb::realtype intArea2;   // area of the interface of 2nd structure

      ProSurf ();
      ProSurf ( mmdb::io::RPStream Object );
      ~ProSurf();

      void  setCodeNo         ( int            code_no     );
      void  setProbeRadius    ( mmdb::realtype probeRadius );
      void  setExcludeSolvent ( bool           exclSolvent );

      //   In the following function, Atom[0..nAtoms-1] should be part of
      // a structure in an MMDB Manager. The Manager is used for
      // bricking, which is disposed on exit, and selection of the
      // surface atoms (from those given in Atom).
      //   If atomSAS is not NULL, it must point on a real vector
      // [0..nAtoms-1], where calcSurface(..) returns individual solvent
      // accessible surface areas for atoms.
      //   Returns:
      //     PROSURF_Ok        success
      //     PROSURF_noAtoms   nAtoms <= 0
      //     PROSURF_noManager Atom do not belong to coordinate
      //                       hierarchy.
      // NOTE 1: if Atom contains only NULL pointers, the function still
      //         returns PROSURF_Ok and a valid selection handle,
      //         however that selection is empty and Area is equal to 0.0
      // NOTE 2: if Atom does not represent all atoms in MMDB Manager M,
      //         other atoms are merely ignored. Thus, selecting only
      //         internal atoms in Atom will yield calculation of surface
      //         around them, even though it is not solvent-accessible
      //         in the reality.
      //   Public fields engaged on output:
      //     selHndSurf1  - selection handle for surface atoms
      //     surfArea1    - SAS area, A^2
      PROSURF_RC calcSurface ( mmdb::PPAtom  atom,
                               int           nAtoms,
                               PMolRefIndex  molRef,
                               mmdb::rvector atomSAS=NULL );

      //   This function is a wrapper over the previous one. All the
      // structure is supplied in MMDB manager M, however any part
      // of it may be selected using CID.
      PROSURF_RC calcSurface ( mmdb::PManager M,
                               PMolRefIndex   molRef,
                               mmdb::cpstr    CID=NULL );


      //  calcInterface(..) will calculate the solvent accessible surface
      // (SAS) of structures Atom1 of length nAtoms1 and Atom2 of length
      // nAtoms2, as well as an interface (defined as SAS buried because
      // of structure contact) between them. Both Atom1 and Atom2 must
      // belong to the same coordinate hierarchy.
      //   surfFlag allows for the following options:
      //    SURFF_Select      will select the surface atoms
      //    SURFF_dontSelect  will not select the surface atoms thus
      //                      saving time and selection space in M.
      //                      The surface selection handles in PID will
      //                      be set to zero
      //    SURFF_Selected    assumes that surface atoms of Atom1 and
      //                      Atom2 are already selected (perhaps by
      //                      a previous call to calcInterface or
      //                      calcSurface), and the selection handles
      //                      are given in the corresponding public
      //                      fields (selHndSurf1 and selHndSurf2).
      //                      This will save time. If selection handles
      //                      are set zero, then the function works as
      //                      if SURFF_Select flag were given.
      //   Returns:
      //     PROSURF_Ok       success
      //     PROSURF_NoAtoms  nAtoms1<=0 or nAtoms2<=0 .
      //   All public fields are engaged on output.
      PROSURF_RC  calcInterface (
                           mmdb::PPAtom   atom1,
                           int            nAtoms1,
                           mmdb::PPAtom   atom2,
                           int            nAtoms2,
                           PMolRefIndex   molRef,
                           SURFSEL_KEY    surfFlag    = SURFF_Select,
                           mmdb::rvector  atomSAS1    = NULL,
                           mmdb::rvector  intAtomSAS1 = NULL,
                           mmdb::rvector  atomSAS2    = NULL,
                           mmdb::rvector  intAtomSAS2 = NULL );

      //   This function is a wrapper over the previous one. All the
      // structure is supplied in MMDB manager M, however any part
      // of it may be selected using CID.
      PROSURF_RC  calcInterface (
                           mmdb::PManager M,
                           mmdb::cpstr    CID1,
                           mmdb::cpstr    CID2,
                           PMolRefIndex   molRef,
                           SURFSEL_KEY    surfFlag=SURFF_Select );

      mmdb::realtype getMaxGap();  // 2*probe_radius

    protected :

      //  Spherical Code
      mmdb::rvector scx,scy,scz,sca; // code points and areas
      int     codeNo,scn;      // code number and number of points
      //  Bricking
      mmdb::realtype max_radius;     // estimate of maximal atom radius
      //  Options
      mmdb::realtype probe_radius;   // probe (solvent) sphere radius
      bool  excludeSolvent; // True by default

      void ResetResults     ();
      void InitProSurf      ();
      void FreeSCMemory     ();
      void calcSphericalCode();
      void GetAtomSet       ( mmdb::PManager M, mmdb::cpstr CID,
                              int & selHnd,     mmdb::PPAtom & Atom,
                              int & nAtoms );

  };


  // ========================  SolvEnergy  ==========================

  enum ASP_TYPE  {
    ASP_Special   = -2,
    ASP_Unknown   = -1,
    ASP_Other     =  0,
    ASP_Carbon    =  1,
    ASP_Sulphur   =  2,
    ASP_NeutralNO =  3,
    ASP_ChargedN  =  4,
    ASP_ChargedO  =  5,
    ASP_C_D_S     =  6,
    ASP_C_D_B     =  7,
    ASP_C_ARO     =  8,
    ASP_C_ALI     =  9,
    ASP_N_ARO     = 10,
    ASP_N_AMI     = 11,
    ASP_O_D_P     = 12,
    ASP_O_D_B     = 13,
    ASP_O_D_S     = 14,
    ASP_O_BYL     = 15,
    ASP_P_DNA     = 16,
    nASPs         = 17   //!< number of ASP types
  };

  extern mmdb::realtype ASP[nASPs];

  DefineClass(SolvEnergy);

  class SolvEnergy : public mmdb::io::Stream  {

    public :
      mmdb::realtype solvEnergy; // solvation energy of folding
      mmdb::realtype DeltaG;     // solvation energy gain of interface
      mmdb::realtype PValue;     // P-value of the interface solv-n energy gain
      mmdb::realtype aveDeltaG;  // average Delta G for the interface
      mmdb::realtype probDeltaG; // probable Delta G for the interface
      mmdb::realtype aspArea[nASPs]; // partial atom surface areas

      SolvEnergy ();
      SolvEnergy ( mmdb::io::RPStream Object );
      ~SolvEnergy();

      mmdb::realtype calcSurfSolvEnergy (
                         mmdb::PPAtom  atom, mmdb::rvector atomSAS,
                         int   nAtoms, mmdb::rvector atomSE,
                         PMolRefIndex  molRef );
      mmdb::realtype calcIntSolvEnergy  (
                         mmdb::PPAtom  atom, mmdb::rvector atomSAS,
                         mmdb::rvector intSAS, int nAtoms,
                         mmdb::rvector atomSE,
                         PMolRefIndex  molRef );
      mmdb::realtype calcIntPValue (
                         mmdb::PPAtom  atom, mmdb::rvector atomSAS,
                         int   nAtoms, int nIntAtoms,
                         mmdb::realtype intDeltaG,
                         PMolRefIndex  molRef );

    protected :
      void ResetResults  ();
      void InitSolvEnergy();

      void calcPValue ( mmdb::rvector surfEn, int nSurfAtoms,
                        int  nIntAtoms, mmdb::realtype intDeltaG );

  };


  extern void   writeASPs ( mmdb::io::RFile f );
  extern void   readASPs  ( mmdb::io::RFile f );

  extern int           aspType ( mmdb::PAtom a );
  extern mmdb::realtype getASP ( mmdb::PAtom a );

  extern void makeASPCIF     ( mmdb::mmcif::PLoop mmCIFLoop );
  extern void writeStructCIF ( mmdb::pstr FName );

}  // namespace pisa

#endif
