// $Id: pisa_asmstock.cpp $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_asmstock <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::AsmStock
//       ~~~~~~~~~
//
//  (C) E. Krissinel, 2004-2014
//
// =================================================================
//

#include "pisa_asmstock.h"
#include "chem_equilibrium.h"
#include "mmdb2/mmdb_math_.h"

namespace pisa  {

  //  =========================  AsmStock  ==========================

  AsmStock::AsmStock ( int n )  {
    InitAsmStock ( n );
  }

  AsmStock::~AsmStock()  {
    FreeMemory();
  }

  void AsmStock::InitAsmStock ( int n )  {
  mmdb::realtype lgC,dLgC;
  int            i,j;

    if (n>0)  {

      A = new PAssembly[n];  // list of assemblies in stock
      for (i=0;i<n;i++)
        A[i] = NULL;

      nConc = 1201;
      mmdb::GetVectorMemory ( asuConc,nConc,0 );
      mmdb::GetMatrixMemory ( asmConc,n,nConc,0,0 );

      lgC  = -12.0;
      dLgC = -lgC/(nConc-1);
      for (i=0;i<nConc;i++)  {
        asuConc[i] = mmdb::math::exp10 ( lgC );
        for (j=0;j<n;j++)
          asmConc[j][i] = 0.0;
        lgC += dLgC;
      }

    } else  {

      A       = NULL;
      asuConc = NULL;
      asmConc = NULL;
      nConc   = 0;

    }

    nAssemblies = n;        // number of assemblies in stock
    nAsmAlloc   = n;
    asuSize     = 1;

  }

  void AsmStock::FreeMemory()  {
  int i;

    if (A)  {
      for (i=0;i<nAsmAlloc;i++)
        if (A[i])  delete A[i];
      delete[] A;
      A = NULL;
    }

    mmdb::FreeVectorMemory ( asuConc,0 );
    mmdb::FreeMatrixMemory ( asmConc,nAsmAlloc,0,0 );

    nAssemblies = 0;
    nAsmAlloc   = 0;
    nConc       = 0;

  }

  AsmStock::RETURN_CODE AsmStock::addToStock ( PAssembly assembly,
                                               mmdb::ivector & count ) {
  PPAssembly    A1;
  mmdb::ivector c1;
  int           i,n1,t;

    t = assembly->type;
    if (t<0)
      return UnknownType;

    if (t>nAssemblies)  {
      n1 = t + 10;
      A1 = new PAssembly[n1];
      mmdb::GetVectorMemory ( c1,n1,0 );
      for (i=0;i<nAssemblies;i++)  {
        A1[i] = A[i];
        c1[i] = count[i];
      }
      for (i=nAssemblies;i<n1;i++)  {
        A1[i] = NULL;
        c1[i] = 0;
      }
      if (A)  delete[] A;
      A = A1;
      mmdb::FreeVectorMemory ( count,0 );
      count = c1;
      nAssemblies = n1;
    }

    if (!A[t])  {
      A[t] = new Assembly();
      count[t] = 0;
    }

    A[t]->averageData ( assembly,count[t] );

    return Ok;

  }

  AsmStock::RETURN_CODE AsmStock::processStock ( PDomains D,
                                           mmdb::realtype minConc,
                                           mmdb::realtype maxConc )  {
  mmdb::ivector     species;
  mmdb::ivector     products;
  chem::Equilibrium equil;
  mmdb::realtype    lgC,dLgC;
  int               i,j,k, nSpecies,nProducts, speciesType;
  bool              B;


    D->getNofDomains ( i,j,k );
    asuSize = i + j;

    lgC  = log10(minConc);
    dLgC = (log10(maxConc)-lgC)/(nConc-1);
    for (i=0;i<nConc;i++)  {
      asuConc[i] = mmdb::math::exp10 ( lgC );
      for (j=0;j<nAssemblies;j++)
        asmConc[j][i] = 0.0;
      lgC += dLgC;
    }

    mmdb::GetVectorMemory ( species ,nAssemblies,0 );
    mmdb::GetVectorMemory ( products,nAssemblies,0 );

    nSpecies  = 0;
    nProducts = 0;
    k         = 0;
    for (i=0;i<nAssemblies;i++)
      if (A[i])  {
        if (A[i]->mmSize>0)  {
          if (A[i]->freeSize<=1)  species [nSpecies++]  = k;
                            else  products[nProducts++] = k;
          if (k<i)  {
            A[k] = A[i];
            A[i] = NULL;
          }
          k++;
        } else  {
          delete A[i];
          A[i] = NULL;
        }

      }

//    for (i=k;i<nAssemblies;i++)
//      A[i] = NULL;

    nAssemblies = k;  // == nSpecies + nProducts

    equil.allocate ( nSpecies,nProducts );

    for (i=0;i<nSpecies;i++)  {
      speciesType = A[species[i]]->getFirstMonType();
      for (j=0;j<nProducts;j++)
        equil.addSpecies ( j,i,
                     A[products[j]]->getMonTypeOccurence(speciesType) );
    }

    B = true;
    for (i=0;i<nProducts;i++)  {
      equil.setLnKdiss ( i,-A[products[i]]->freeEn0/TRconst );
      if (equil.productSize(i)!=A[products[i]]->mmSize)
        B = false;
    }

    if (!B)  {
      mmdb::FreeVectorMemory ( species ,0 );
      mmdb::FreeVectorMemory ( products,0 );
      return WrongProductSize;
    }

    for (i=0;i<nSpecies;i++)
      equil.setStoichiometry ( i,
              D->getMonTypeOccurence(A[species[i]]->getFirstMonType()) );

    equil.prepareSolution();

    for (i=0;i<nConc;i++)  {
      equil.equilibrate ( asuConc[i],1.0e-12 );
      for (j=0;j<nSpecies;j++)
        asmConc[species[j]][i] = equil.speciesConc(j);
      for (j=0;j<nProducts;j++)
        asmConc[products[j]][i] = equil.productConc(j);
    }

    mmdb::FreeVectorMemory ( species ,0 );
    mmdb::FreeVectorMemory ( products,0 );

    return Ok;

  }

  mmdb::realtype AsmStock::getRelativeConc ( int n )  {
    if ((!asuConc) || (n>=nConc))  return 0.0;
    return asuConc[n]/asuConc[nConc-1];
  }

  mmdb::realtype AsmStock::getAggregationIndex ( int asmNo, int n )  {

    if ((!asuConc) || (n>=nConc) || (asmNo>=nAssemblies))
      return 0.0;

    return A[asmNo]->mmSize*asmConc[asmNo][n]/asuSize/asuConc[n];

  }

  void AsmStock::write ( mmdb::io::RFile f, int nInterfaces )  {
  int        i,j;
  mmdb::byte Version;

    Version = 1;
    f.WriteByte ( &Version );

    f.WriteInt ( &nAssemblies );
    for (i=0;i<nAssemblies;i++)
      A[i]->write ( f,nInterfaces );

    f.WriteInt ( &nConc );
    for (i=0;i<nConc;i++)  {
      f.WriteReal ( &(asuConc[i]) );
      for (j=0;j<nAssemblies;j++)
        f.WriteReal ( &(asmConc[j][i]) );
    }

    f.WriteInt ( &asuSize );

  }

  void AsmStock::read ( mmdb::io::RFile f, int nInterfaces )  {
  int  i,j;
  mmdb::byte Version;

    FreeMemory();

    f.ReadByte ( &Version );

    f.ReadInt ( &nAssemblies );
    if (nAssemblies>0)  {
      A = new PAssembly[nAssemblies];
      for (i=0;i<nAssemblies;i++)  {
        A[i] = new Assembly();
        A[i]->read ( f,nInterfaces );
      }
    }
    nAsmAlloc = nAssemblies;

    f.ReadInt ( &nConc );
    if (nConc>0)  {
      mmdb::GetVectorMemory ( asuConc,nConc,0 );
      mmdb::GetMatrixMemory ( asmConc,nAsmAlloc,nConc,0,0 );
      for (i=0;i<nConc;i++)  {
        f.ReadReal ( &(asuConc[i]) );
        for (j=0;j<nAssemblies;j++)
          f.ReadReal ( &(asmConc[j][i]) );
      }
    }

    f.ReadInt ( &asuSize );

  }


}  // namespace pisa
