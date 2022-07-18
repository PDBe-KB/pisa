// $Id: pisa_asmunit.cpp $
// =================================================================
//
//    20.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_asmunit <implementation>
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

#include <string.h>

#include "pisa_asmunit.h"
#include "pisa_defs.h"

namespace pisa  {

  // =========================  AsmUnit  ===========================

  AsmUnit::AsmUnit()  {
  int i;

    id     = 0;            // domain id (serial number in asu)
    type   = 0;            // domain type id
    dclass = DCLASS_None;  // domain class id

    for (i=0;i<_max_n_int;i++)
      intfl[i] = NULL;  // lists of interfaces (>=0: AsmUnit No
                        //   in assembly; -1: link to outside;
                        //                -2: not present)

    dissAsm    = 0;     // dissociated assembly no.
    symOpNo    = 0;     // symop (ASU) serial number
    rcsb_symop = 0;     // rcsb symop serial number
    symOp      = NULL;  // symmetry operation
    ncsOpNo    = 0;     // NCS operation serial number
    cell_i     = 0;
    cell_j     = 0;
    cell_k     = 0;     // unit cell ( -n..0..n );

    fixed      = false; // true if a fixed ligand

    visualID[0] = char(0);

  }

  AsmUnit::~AsmUnit()  {
    FreeMemory();
  }

  void AsmUnit::FreeMemory()  {
  int i;
    for (i=0;i<_max_n_int;i++)
      mmdb::FreeVectorMemory ( intfl[i],0 );
    if (symOp)  {
      delete[] symOp;
      symOp = NULL;
    }
  }

  void AsmUnit::makeUnit ( RMonRef R, PPInterface interface,
                           int nInterfaces, int nInterfaces0,
                           PPDomain D, mmdb::mat44 & rom )  {
  PMonomer  M1;
  PMultimer U;
  int       i,j,k,n, dic,djc,dkc, intNo;

    FreeMemory();

    id     = R.M->id;
    type   = D[R.M->ncsParent]->type;
    dclass = D[R.M->ncsParent]->dclass;

    if (nInterfaces>0)  {
      U = R.M->U;
      for (i=0;i<nInterfaces;i++)  {
        intNo = interface[i]->id-1;
        n = 0;
        while (n<_max_n_int)
          if (R.M->L[n])  {
            if (!intfl[n])  {
              mmdb::GetVectorMemory ( intfl[n],nInterfaces0,0 );
              for (j=0;j<nInterfaces0;j++)
                intfl[n][j] = -2;
            }
            M1 = R.M->L[n][i].M;
            if (M1)  {
              if (M1->id2>0)  {
                // check that L[n][i] links to the right unit cell
                k = M1->id2 - 1;
                dic = R.i + R.M->L[n][i].i - U->R[k].i;
                djc = R.j + R.M->L[n][i].j - U->R[k].j;
                dkc = R.k + R.M->L[n][i].k - U->R[k].k;
                if (dic || djc || dkc)
                      intfl[n][intNo] = -1;  // link to outside
                else  intfl[n][intNo] = k;
              } else // link to outside the assembly
                intfl[n][intNo] = -1;
            }
            n++;
          } else
            n = _max_n_int;

      }

    }

    symOpNo    = R.M->symOpNo; // symop (ASU) serial number
    rcsb_symop = symOpNo;      // just for the moment
    ncsOpNo    = 0;            // temporary value for NCS operation number

    cell_i = R.i;
    cell_j = R.j;
    cell_k = R.k;  // unit cell ( -n..0..n );

    for (i=0;i<4;i++)
      for (j=0;j<4;j++)
        T[i][j] = R.M->uct[i][j];
    T[0][3] += R.i*rom[0][0] + R.j*rom[0][1] + R.k*rom[0][2];
    T[1][3] += R.i*rom[1][0] + R.j*rom[1][1] + R.k*rom[1][2];
    T[2][3] += R.i*rom[2][0] + R.j*rom[2][1] + R.k*rom[2][2];

  }

  mmdb::xml::PXMLObject AsmUnit::getAsmUnitXML ( PPDomain D,
                                                 int nCellOut )  {
  mmdb::xml::PXMLObject xml;
  char        S[500];

    xml = new mmdb::xml::XMLObject ( xml_asmunit );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asmunit_name,
                                       D[id]->getDomainID(S)) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asmunit_visual_id,visualID) );

    //addXML0 ( xml,xml_asmunit_rxx,T[0][0] );
    //addXML0 ( xml,xml_asmunit_rxy,T[0][1] );
    //addXML0 ( xml,xml_asmunit_rxz,T[0][2] );
    //addXML0 ( xml,xml_asmunit_tx ,T[0][3] );
    //addXML0 ( xml,xml_asmunit_ryx,T[1][0] );
    //addXML0 ( xml,xml_asmunit_ryy,T[1][1] );
    //addXML0 ( xml,xml_asmunit_ryz,T[1][2] );
    //addXML0 ( xml,xml_asmunit_ty ,T[1][3] );
    //addXML0 ( xml,xml_asmunit_rzx,T[2][0] );
    //addXML0 ( xml,xml_asmunit_rzy,T[2][1] );
    //addXML0 ( xml,xml_asmunit_rzz,T[2][2] );
    //addXML0 ( xml,xml_asmunit_tz ,T[2][3] );

    //addXML0 ( xml,xml_asmunit_frxx,TF[0][0] );
    //addXML0 ( xml,xml_asmunit_frxy,TF[0][1] );
    //addXML0 ( xml,xml_asmunit_frxz,TF[0][2] );
    //addXML0 ( xml,xml_asmunit_ftx ,TF[0][3] );
    //addXML0 ( xml,xml_asmunit_fryx,TF[1][0] );
    //addXML0 ( xml,xml_asmunit_fryy,TF[1][1] );
    //addXML0 ( xml,xml_asmunit_fryz,TF[1][2] );
    //addXML0 ( xml,xml_asmunit_fty ,TF[1][3] );
    //addXML0 ( xml,xml_asmunit_frzx,TF[2][0] );
    //addXML0 ( xml,xml_asmunit_frzy,TF[2][1] );
    //addXML0 ( xml,xml_asmunit_frzz,TF[2][2] );
    //addXML0 ( xml,xml_asmunit_ftz ,TF[2][3] );

    sprintf ( S,"%i_%1i%1i%1i",rcsb_symop,
              cell_i+nCellOut,cell_j+nCellOut,cell_k+nCellOut );

    xml->AddObject ( new mmdb::xml::XMLObject(xml_asmunit_sym_id,S) );

    return xml;

  }

  void AsmUnit::Copy ( PAsmUnit M )  {

    FreeMemory();

    id     = M->id;      // domain id (serial number in asu)
    type   = M->type;    // domain type id
    dclass = M->dclass;  // domain class id
    /*
      mmdb::ivector intfl[_max_n_int]; //!< lists of interfaces (>=0: AsmUnit
                       //!<   No in assembly; -1: link to outside;
                       //!<                   -2: not present)
    */
    dissAsm    = M->dissAsm;     // dissociated assembly no. 1,2,...
    symOpNo    = M->symOpNo;     // symop (ASU) serial number
    rcsb_symop = M->rcsb_symop;  // rcsb symop serial number
    mmdb::CreateCopy ( symOp,M->symOp );  // symmetry operation
    ncsOpNo    = M->ncsOpNo;     // NCS opeartion serial number
    cell_i     = M->cell_i;
    cell_j     = M->cell_j;
    cell_k     = M->cell_k;      // unit cell ( -n..0..n );
    mmdb::Mat4Copy ( M->T ,T  ); // transformation matrix
    mmdb::Mat4Copy ( M->TF,TF ); // transf. matrix in fract. coordinates
                                 // chain ID for visual-n and download
    strncpy ( visualID,M->visualID,sizeof(mmdb::ChainID) );

    fixed = M->fixed;

  }

  void AsmUnit::write ( mmdb::io::RFile f, int nInterfaces )  {
  int  i,j;
  mmdb::byte Version,t;

    Version = 6;
    f.WriteByte ( &Version );

    f.WriteInt ( &id     );
    f.WriteInt ( &type   );
    f.WriteInt ( &dclass );

    t = 0x01;
    i = 0;
    while (i<_max_n_int)
      if (intfl[i])  {
        f.WriteByte ( &t );
        for (j=0;j<nInterfaces;j++)
          f.WriteInt ( &(intfl[i][j]) );
        i++;
      } else
        i = _max_n_int;
    t = 0x00;
    f.WriteByte   ( &t );

    f.WriteInt    ( &dissAsm    );
    f.WriteInt    ( &symOpNo    );
    f.WriteInt    ( &rcsb_symop );
    f.CreateWrite ( symOp       );
    f.WriteInt    ( &ncsOpNo    );
    f.WriteInt    ( &cell_i     );
    f.WriteInt    ( &cell_j     );
    f.WriteInt    ( &cell_k     );
    for (i=0;i<3;i++)
      for (j=0;j<4;j++)
        f.WriteReal ( &(T[i][j]) );

    for (i=0;i<3;i++)
      for (j=0;j<4;j++)
        f.WriteReal ( &(TF[i][j]) );

    f.WriteTerLine ( visualID );
    f.WriteBool    ( &fixed   );

  }

  void AsmUnit::read  ( mmdb::io::RFile f, int nInterfaces )  {
  int        i,j;
  mmdb::byte Version,t;

    FreeMemory();

    f.ReadByte ( &Version );

    f.ReadInt ( &id );
    if (Version>5)  {
      f.ReadInt ( &type   );
      f.ReadInt ( &dclass );
    }

    if (Version>1)  {
      i = 0;
      while (i<_max_n_int)  {
        f.ReadByte ( &t );
        if (t)  {
          mmdb::GetVectorMemory ( intfl[i],nInterfaces,0 );
          for (j=0;j<nInterfaces;j++)
            f.ReadInt ( &(intfl[i][j]) );
          i++;
        } else
          i = _max_n_int;
      }
    } else if (nInterfaces>0)  {
      mmdb::GetVectorMemory ( intfl[0],nInterfaces,0 );
      mmdb::GetVectorMemory ( intfl[1],nInterfaces,0 );
      for (i=0;i<nInterfaces;i++)  {
        f.ReadInt ( &(intfl[0][i]) );
        f.ReadInt ( &(intfl[1][i]) );
      }
    }

    f.ReadInt    ( &dissAsm );
    f.ReadInt    ( &symOpNo );
    if (Version>3)
         f.ReadInt ( &rcsb_symop );
    else rcsb_symop = symOpNo;
    f.CreateRead ( symOp    );
    f.ReadInt    ( &ncsOpNo );
    f.ReadInt    ( &cell_i  );
    f.ReadInt    ( &cell_j  );
    f.ReadInt    ( &cell_k  );
    for (i=0;i<3;i++)
      for (j=0;j<4;j++)
        f.ReadReal ( &(T[i][j]) );
    T[3][0] = 0.0;
    T[3][1] = 0.0;
    T[3][2] = 0.0;
    T[3][3] = 1.0;

    if (Version>2) {
      for (i=0;i<3;i++)
        for (j=0;j<4;j++)
          f.ReadReal ( &(TF[i][j]) );
      TF[3][0] = 0.0;
      TF[3][1] = 0.0;
      TF[3][2] = 0.0;
      TF[3][3] = 1.0;
    }

    if (Version>4)  f.ReadTerLine ( visualID );
              else  visualID[0] = char(0);

    if (Version>5)
      f.ReadBool ( &fixed );

  }


}  // namespace pisa
