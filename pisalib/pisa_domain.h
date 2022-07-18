// $Id: pisa_domain.h $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_domain <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Domain
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2014
//
// =================================================================
//

#ifndef __PISA_Domain__
#define __PISA_Domain__

#include "mmdb2/mmdb_manager.h"
#include "ssm/ssm_graph.h"

#include "pisa_prosurf.h"
#include "pisa_types.h"

namespace pisa  {

  //  ========================  Domain  =============================

  DefineClass(Domains);

  DefineClass(Domain);
  DefineStreamFunctions(Domain);

  const int frameLen = 4;

  typedef mmdb::vect3 TFrame[frameLen];
  typedef TFrame & RTFrame;
  typedef TFrame * PTFrame;

  typedef mmdb::realtype DFrame[frameLen][frameLen];
  typedef DFrame & RDFrame;
  typedef DFrame * PDFrame;


  class Domain : public mmdb::io::Stream  {

    public :
      mmdb::pstr     range;      //!< domain range (selection of residues)
      mmdb::pstr     typeId;     //!< domain type Id for output
      int            ncsOpNo;    //!< NCS operation 0..nNCSOps, 0<->identity
      int            ncsParent;  //!< parential domain for NCS mates
      int            selHndSurf; //!< selection handle for the domain surface
      int            nAtoms;     //!< total number of atoms in the domain
      int            nRes;       //!< total number of residues in the domain
      int            nRes0;      //!< total number of AA/NT res-s in the domain
      int            nSurfAtoms; //!< number of surafce atoms in the domain
      int            nSurfRes;   //!< number of surface residues in the domain
      int            type;       //!< domain type 1,2,...
      DOMAIN_CLASS   dclass;     //!< domain class protein/dna/ligand DCLASS_XXXX
      int            symNumber;  //!< symmetry number
      bool           assemble;   //!< True if to be inculded into assembling
      bool           agent;      //!< True if precipitation agent
      mmdb::mat44    ncs_m;      //!< NCS matrix
      mmdb::mat33    itm;        //!< inertia tensor matrix in mass center
      mmdb::vect3    im;         //!< principal moments of inertia in mass center
      mmdb::mat33    itb;        //!< canonical coordinate reference
      TFrame         frame;      //!< reference frame
      mmdb::realtype surfArea;   //!< surface area
      mmdb::realtype DeltaG;     //!< relative solvation energy
      mmdb::realtype mx,my,mz;   //!< center of mass
      mmdb::realtype weight;     //!< molecular weight
      mmdb::realtype R;          //!< minimal radius of outer sphere
      mmdb::rvector  asa;        //!< vector[0..nRes-1] of accessible surf. area
                                 ///   of individual residues
      mmdb::rvector  solvEn;     //!< vector[0..nRes-1] of sol-n energy changes
                                 ///   of individual residues

      Domain ();
      Domain ( mmdb::io::RPStream Object );
      ~Domain();

      void SetChainRange     ( mmdb::cpstr    chainID );
      void SetLigandRange    ( mmdb::PResidue res     );
      void SetLigandResidues ( mmdb::cpstr    ligRes  );

      int SelectDomain ( mmdb::PManager       MMDB,
                         mmdb::SELECTION_TYPE selType,
                         int   modelNo,
                         bool  removeHydrogens );
      int SelectDomain ( int   selHnd,
                         mmdb::PManager       MMDB,
                         mmdb::SELECTION_TYPE selType,
                         mmdb::SELECTION_KEY  selKey,
                         int   modelNo,
                         bool  removeHydrogens );

      void  makeReferenceFrame ();
      void  getReferenceFrame  ( RTFrame F, mmdb::mat44 & TM );
      void  setReferenceFrame  ( RTFrame F  );

      void  calcDimensions    ( mmdb::PManager MMDB );
      void  getInertiaTensor  ( mmdb::mat33 & itn,
                                mmdb::mat44 & TMatrix,
                                mmdb::realtype x,
                                mmdb::realtype y,
                                mmdb::realtype z );

      bool  isContact      ( PDomain D );
      bool  isContact      ( PDomain D, mmdb::mat44 & TM );
      bool  isContact      ( mmdb::mat44 & T, PDomain D,
                             mmdb::mat44 & TM );

      void  MakeContactBricks ( mmdb::PManager MMDB,
                                mmdb::mat44 & TMatrix,
                                int & selHnd );
      bool  isContact         ( mmdb::mat44 & TM, PDomain D,
                                mmdb::PManager MMDB, int & selHnd );


      PROSURF_RC calcSurface ( mmdb::PManager MMDB,
                               PProSurf       ProSurf,
                               PMolRefIndex   MolRef );

      mmdb::pstr  getDomainRange      ( mmdb::pstr S, int fieldLen );
      mmdb::pstr  getFileSafeName     ( mmdb::pstr S, int fieldLen );
      mmdb::pstr  getDomainClass      ( mmdb::pstr S );
      mmdb::pstr  getLigandRange      ( mmdb::pstr S );
      mmdb::pstr  getChainID          ( mmdb::pstr S );
      mmdb::pstr  getDomainID         ( mmdb::pstr S );
      mmdb::pstr  getDomainRange_html ( mmdb::pstr S,
                                        bool separate=false );

      bool checkLigandRange ( PDomain D );

      void CopyParentData ( PPDomain Domain );

      void write ( mmdb::io::RFile f );
      void read  ( mmdb::io::RFile f );

    protected :
      mmdb::pstr ligandRes;  //!< for protein/rna/dna, list of residues
                             /// that were recognized as ligands

      void InitDomain();
      void FreeMemory();
      void getCID    ( mmdb::pstr chID, mmdb::pstr resName, int & resNo,
                       mmdb::pstr insCode );

  };


}  // namespace pisa

#endif
