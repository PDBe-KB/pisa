// $Id: pisa_rcsbdat.cpp $
// =================================================================
//
//    19.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// -----------------------------------------------------------------
//
//  **** Module  :  PISA_RCSBDat <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::RCSBData
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2013
//
// =================================================================
//

#include <string.h>

#include "pisa_rcsbdat.h"

namespace pisa  {

  //  ==========================  RCSBData  =========================

  RCSBData::RCSBData() : mmdb::io::Stream()  {
    InitRCSBData ( "",0 );
  }

  RCSBData::RCSBData ( mmdb::cpstr spGroup, int nSO )
          : mmdb::io::Stream()  {
    InitRCSBData ( spGroup,nSO );
  }

  RCSBData::RCSBData ( mmdb::io::RPStream Object )
          : mmdb::io::Stream ( Object )  {
    InitRCSBData ( "",0 );
  }

  RCSBData::~RCSBData()  {
    FreeMemory();
  }

  void RCSBData::InitRCSBData ( mmdb::cpstr spGroup, int nSO )  {
  int i;
    strcpy ( spaceGroup,spGroup );
    nSymOps = nSO;
    if (nSymOps>0)  {
      mmdb::GetVectorMemory ( rcsb_symop,nSO,0 );
      for (i=0;i<nSymOps;i++)
        rcsb_symop[i] = -1;
    } else
      rcsb_symop = NULL;

  }

  void RCSBData::FreeMemory()  {
    mmdb::FreeVectorMemory ( rcsb_symop,0 );
    spaceGroup[0] = char(0);
    nSymOps       = 0;
  }

  void RCSBData::write ( mmdb::io::RFile f )  {
  int i;
    f.WriteTerLine ( spaceGroup,false );
    f.WriteInt     ( &nSymOps         );
    for (i=0;i<nSymOps;i++)
      f.WriteInt ( &(rcsb_symop[i]) );
  }

  void RCSBData::read ( mmdb::io::RFile f )  {
  int i,l;
    l = nSymOps;
    f.ReadTerLine ( spaceGroup,false );
    f.ReadInt     ( &nSymOps         );
    if (nSymOps>l)  {
      mmdb::FreeVectorMemory ( rcsb_symop,0 );
      mmdb::GetVectorMemory  ( rcsb_symop,nSymOps,0 );
    }
    for (i=0;i<nSymOps;i++)
      f.ReadInt ( &(rcsb_symop[i]) );
  }



  int getRCSBData ( mmdb::pstr fileName, mmdb::pstr spaceGroup,
                    RPRCSBData rcsbData )  {
  mmdb::io::File  f;
  int             nRCSBData,i;
  bool            done;

    if (rcsbData)  {
      delete rcsbData;
      rcsbData = NULL;
    }

    if (!spaceGroup)  return 2;

    if (!spaceGroup[0])  return 3;

    f.assign ( fileName,false,true );
    if (f.reset(true))  {

      f.ReadInt ( &nRCSBData );
      done = false;
      for (i=0;(i<nRCSBData) && (!done);i++)  {
        StreamRead ( f,rcsbData );
        if (!strcmp(spaceGroup,rcsbData->spaceGroup))
          done = true;
      }
      f.shut();

      if (!done)  {
        if (rcsbData)  {
          delete rcsbData;
          rcsbData = NULL;
        }
        return 1;
      } else
        return 0;

    } else
      return -1;

  }

  int getRCSBData ( mmdb::pstr fileName, PPRCSBData & rcsbData,
                    int & nRCSBData )  {
   mmdb::io::File f;
  int   i;

    deleteRCSBData ( rcsbData,nRCSBData );

    f.assign ( fileName,false,true );
    if (f.reset(true))  {

      f.ReadInt ( &nRCSBData );
      rcsbData = new PRCSBData[nRCSBData];
      for (i=0;i<nRCSBData;i++)  {
        rcsbData[i] = NULL;
        StreamRead ( f,rcsbData[i] );
      }
      f.shut();

      return 0;

    } else
      return -1;

  }

  void deleteRCSBData  ( PPRCSBData & rcsbData, int nRCSBData )  {
  int i;

    if (rcsbData)  {

      for (i=0;i<nRCSBData;i++)
        if (rcsbData[i])  delete rcsbData[i];
      delete[] rcsbData;
      rcsbData = NULL;

    }

    nRCSBData = 0;

  }

  MakeStreamFunctions(RCSBData)

}  // namespace pisa
