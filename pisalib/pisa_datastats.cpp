// $Id: pisa_datastats.h $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_datastats <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  DataStats
//       ~~~~~~~~~
//
//  (C) E. Krissinel, 2013-2014
//
// =================================================================
//

#include <stdlib.h>
#include <string.h>

#include "pisa_datastats.h"
#include "mmdb2/mmdb_utils.h"


//  ======================  StatDistribution  ======================

pisa::StatDistribution::StatDistribution()  {
  parts       = NULL;
  prob        = NULL;
  in_assembly = NULL;
  nParts      = 0;
  nPartsAlloc = 0;
}

pisa::StatDistribution::~StatDistribution()  {
  clear();
}

void pisa::StatDistribution::clear()  {
  if (parts)  delete[] parts;
  if (prob)   delete[] prob;
  mmdb::FreeVectorMemory ( in_assembly,0 );
  parts       = NULL;
  prob        = NULL;
  nParts      = 0;
  nPartsAlloc = 0;
}

void pisa::StatDistribution::print()  {

  printf (
    "\n"
    " ------------------------------------------------------------\n"
    " nParts = %i\n",
    nParts
    );

  for (int i=0;i<=nParts;i++)
    printf ( " %5i  %12.7g %12.7g\n",i,parts[i],prob[i] );

}


void pisa::StatDistribution::addData ( mmdb::realtype v,
                                       bool  inAssembly )  {
mmdb::shortreal *p1;
mmdb::ovector    a1;

  if (nParts>=nPartsAlloc)  {
    nPartsAlloc += 5000;
    p1 = new mmdb::shortreal[nPartsAlloc+1];
    mmdb::GetVectorMemory ( a1,nPartsAlloc+1,0 );
    for (int i=0;i<nParts;i++)  {
      p1[i] = parts[i];
      a1[i] = in_assembly[i];
    }
    if (parts)  delete[] parts;
    mmdb::FreeVectorMemory ( in_assembly,0 );
    parts       = p1;
    in_assembly = a1;
  }

  parts      [nParts] = v;
  in_assembly[nParts] = inAssembly;
  nParts++;

}


class SortVector : public mmdb::QuickSort  {
  public :
    SortVector () : QuickSort() { data2 = NULL; }
    ~SortVector() {}
    virtual int  Compare ( int i, int j );
    virtual void Swap    ( int i, int j );
    void Sort ( mmdb::shortreal *sortdata, mmdb::ovector auxdata,
                int data_len );
  protected:
    mmdb::ovector data2;
};


int SortVector::Compare ( int i, int j )  {
// sort by increasing data[i]
  if (((mmdb::shortreal*)data)[i]<((mmdb::shortreal*)data)[j]) return -1;
  if (((mmdb::shortreal*)data)[i]>((mmdb::shortreal*)data)[j]) return  1;
  return 0;
}

void SortVector::Swap ( int i, int j )  {
mmdb::shortreal b;
bool            a;
  b = ((mmdb::shortreal*)data)[i];
  ((mmdb::shortreal*)data)[i] = ((mmdb::shortreal*)data)[j];
  ((mmdb::shortreal*)data)[j] = b;
  a = data2[i];
  data2[i] = data2[j];
  data2[j] = a;
}

void SortVector::Sort ( mmdb::shortreal *sortdata,
                        mmdb::ovector    auxdata, int data_len )  {
  data2 = auxdata;
  mmdb::QuickSort::Sort ( sortdata,data_len );
}


void pisa::StatDistribution::makeParts ( int maxNParts )  {
SortVector       SV;
mmdb::shortreal *p;
mmdb::ivector    count;
mmdb::realtype   d,d0;
int              n, partSize, np,i,k,kp;

  n = nParts;
  if (n<=0)
    return;

  SV.Sort ( parts,in_assembly,n );

  np = mmdb::IMin ( maxNParts,n );
  partSize = mmdb::ifloor ( mmdb::realtype(n)/np );
  if (partSize*np<n)  partSize++;

  if (prob)  delete[] prob;
  prob = new mmdb::shortreal[np];
  mmdb::GetVectorMemory ( count,np,0 );
  p    = new mmdb::shortreal[np+1];

  nParts = 0;
  p[nParts] = parts[0];  // lower boundary belongs to the part
  if (in_assembly[0])  kp = 1;
                 else  kp = 0;
  k = 1;
  for (i=1;i<n;i++)  {
    if ((i==n-1) || ((k>=partSize) && (parts[i]>parts[i-1])))  {
      //  prob[i] contains probability for interface to be internal
      //  to an assembly, if reference parameter is between p[i]
      //  and p[i+1]. count[i] is number of counts used to calculate
      //  that probability.
      prob [nParts]  = kp;
      prob [nParts] /= k;
      count[nParts]  = k;
      nParts++;
      p[nParts] = parts[i];
      k  = 0;
      kp = 0;
    }
    if (in_assembly[i])  kp++;
    k++;
  }

  // In the following loop we optimise partitions such that empty
  // partitions are merged with their neighbours
  do  {
    n  = nParts-1;
    d0 = 0.01;
    k  = -1;
    for (i=0;i<n;i++)  {
      d = fabs(prob[i+1]-prob[i]);
      if (d<d0)  {
        d0 = d;
        k = i;
      }
    }
    if (k>=0)  {
      prob [k] = (count[k]*prob[k] + count[k+1]*prob[k+1]) /
                 (count[k] + count[k+1]);
      count[k] = count[k] + count[k+1];
      p    [k] = p[k+1];
      for (i=k+1;i<n;i++)  {
        prob [i] = prob [i+1];
        count[i] = count[i+1];
        p    [i] = p    [i+1];
      }
      p[n] = p[nParts];
      nParts = n;
    }
  } while (k>=0);

/*
  do  {
    n      = nParts;
    nParts = 0;
    i      = 0;
    while (i<n)  {
      k = i+1;
      if (fabs(prob[k]-prob[i])<0.01)  {
        prob [nParts] = (count[i]*prob[i] + count[k]*prob[k]) /
                        (count[i]+count[k]);
        count[nParts] = count[i] + count[k];
        p    [nParts] = p[k];
        i++;
      } else  {
        prob [nParts] = prob [i];
        count[nParts] = count[i];
        p    [nParts] = p    [i];
      }
      i++;
      nParts++;
    }
    p[nParts] = p[n];
  } while (nParts<n);
*/

/*
  for (i=1;i<nParts;i++)
    prob[i] += prob[i-1];

  for (i=0;i<nParts;i++)
    prob[i] /= n;
*/

  if (parts)  delete[] parts;
  mmdb::FreeVectorMemory ( count,0 );
  mmdb::FreeVectorMemory ( in_assembly,0 );
  parts = p;

}


mmdb::realtype pisa::StatDistribution::getStats ( mmdb::realtype v )  {
int  i,i1,i2;

  if (nParts==0)        return 0.5;
  if (nParts==1)        return prob[0];
  if (v<parts[0])       return prob[0];
  if (v>=parts[nParts]) return prob[nParts-1];

  i1 = 0;
  i2 = nParts;
  while (i2-i1>1)  {
    i = (i2+i1)/2;
    if (parts[i]>v)      i2 = i;
    else if (parts[i]<v) i1 = i;
    else {
      i1 = i;
      i2 = i;
    }
  }

  return prob[i1];

/*
  if (i1==0)         return prob[0]/2.0;
  if (i1==nParts-1)  return prob[i1-1] + (1.0-prob[i1-1])/2.0;
  return prob[i1-1] + (prob[i1+1]-prob[i1])/2.0;
*/

/*
  mmdb::realtype r;
  if (i1==0)         r = prob[0]/2.0;
  else if (i1==nParts-1)  r = prob[i1-1] + (1.0-prob[i1-1])/2.0;
  else r = prob[i1-1] + (prob[i1+1]-prob[i1])/2.0;

  if (r<0.0)  {
    printf ( " ------------------\n" );
    printf ( "negative 2: nParts=%i, i1=%i\n",nParts,i1 );
    if (i1>0)  printf ( "     -1: %10.4g",prob[i1-1] );
    printf ( ";  0: %10.4g",prob[i1] );
    if (i1<nParts)  printf ( ";  +1: %10.4g",prob[i1+1] );
    printf ( "\n" );
  }

  return r;
*/

}

mmdb::realtype pisa::StatDistribution::getMinValue()  {
  if (nParts<0)  return -mmdb::MaxReal;
  return parts[0];
}

mmdb::realtype pisa::StatDistribution::getMaxValue()  {
  if (nParts<0)  return mmdb::MaxReal;
  return parts[nParts];
}


void pisa::StatDistribution::write ( mmdb::io::File & f )  {
mmdb::realtype b;

  f.WriteInt ( &nParts );
  for (int i=0;i<=nParts;i++)  {
    b = parts[i];
    f.WriteFloat ( &b );
  }
  for (int i=0;i<nParts;i++)  {
    b = prob[i];
    f.WriteFloat ( &b );
  }

}

void pisa::StatDistribution::read ( mmdb::io::File & f )  {
mmdb::realtype b;

  clear();

  f.ReadInt ( &nParts       );
  parts = new mmdb::shortreal[nParts+1];
  prob  = new mmdb::shortreal[nParts];
  for (int i=0;i<=nParts;i++)  {
    f.ReadFloat ( &b );
    parts[i] = b;
  }
  for (int i=0;i<nParts;i++)  {
    f.ReadFloat ( &b );
    prob[i] = b;
  }

}


//  ========================  StatDistSet  =========================

pisa::StatDistSet::StatDistSet()  {
  SD      = NULL;
  nParams = 0;
}

pisa::StatDistSet::~StatDistSet()  {
  clear();
}

void pisa::StatDistSet::clear() {
  if (SD)  {
    for (int i=0;i<nParams;i++)
      if (SD[i])  delete SD[i];
    delete[] SD;
    SD = NULL;
  }
  nParams = 0;
}


void pisa::StatDistSet::initSet ( int n_params, int param_no )  {

  clear();

  nParams = n_params;
  SD      = new PStatDistribution[nParams];
  for (int i=0;i<nParams;i++)
    if (i!=param_no)  SD[i] = new StatDistribution();
                else  SD[i] = NULL;

}

void pisa::StatDistSet::addData ( mmdb::rvector p, bool inAssembly )  {
  for (int i=0;i<nParams;i++)
    if (SD[i])
      SD[i]->addData ( p[i],inAssembly );
}

void pisa::StatDistSet::makeParts ( int maxNParts )  {
  for (int i=0;i<nParams;i++)
    if (SD[i])
      SD[i]->makeParts ( maxNParts );
}

void pisa::StatDistSet::getStats ( mmdb::rvector stats )  {
  for (int i=0;i<nParams;i++)
    if (SD[i])
      stats[i] = SD[i]->getStats ( stats[i] );
}

void pisa::StatDistSet::write ( mmdb::io::File & f )  {
int  i;
bool b;

  f.WriteInt ( &nParams );
  for (i=0;i<nParams;i++)  {
    b = (SD[i]!=NULL);
    f.WriteBool ( &b );
    if (SD[i])
      SD[i]->write ( f );
  }

}

void pisa::StatDistSet::read ( mmdb::io::File & f )  {
int  i;
bool b;

  clear();

  f.ReadInt ( &nParams );
  SD = new PStatDistribution[nParams];
  for (i=0;i<nParams;i++)  {
    f.ReadBool ( &b );
    if (b)  {
      SD[i] = new StatDistribution();
      SD[i]->read ( f );
    } else
      SD[i] = NULL;
  }

}


//  ==========================  DataScan  ==========================

pisa::DataScan::DataScan()  {
  distSet   = NULL;
  partition = NULL;
  nParts    = 0;
}

pisa::DataScan::~DataScan()  {
  clear();
}

void pisa::DataScan::clear()  {
  if (partition)  {
    delete[] partition;
    partition = NULL;
  }
  if (distSet)  {
    for (int i=0;i<nParts;i++)
      if (distSet[i])  delete distSet[i];
    delete[] distSet;
    distSet = NULL;
  }
  nParts = 0;
}

void pisa::DataScan::initScan ( int n_parts, mmdb::rvector parts,
                                int nParams, int paramNo )  {
int i;

  clear();

  nParts = n_parts;
  partition = new mmdb::shortreal[nParts+1];
  for (i=0;i<=nParts;i++)
    partition[i] = parts[i];

  distSet = new PStatDistSet[nParts];
  for (i=0;i<nParts;i++)  {
    distSet[i] = new StatDistSet();
    distSet[i]->initSet ( nParams,paramNo );
  }

  scanParam = paramNo;

}

int pisa::DataScan::getPartNo ( mmdb::realtype value )  {
int i1,i2,i;

  i1 = 0;
  i2 = nParts;
  while (i2-i1>1)  {
    i = (i1+i2)/2;
    if (value<partition[i])      i2 = i;
    else if (value>partition[i]) i1 = i;
    else  {
      i1 = i;
      i2 = i;
    }
  }

  if (i1==nParts)  i1--;
  return i1;

}

void pisa::DataScan::addData ( mmdb::rvector p, bool inAssembly )  {
int n = getPartNo ( p[scanParam] );
  distSet[n]->addData ( p,inAssembly );
}

void pisa::DataScan::makeParts ( int maxNParts )  {
  for (int i=0;i<nParts;i++)
    distSet[i]->makeParts ( maxNParts );
}

void pisa::DataScan::printPartition()  {
  printf ( "\n ---- partition\n");
  for (int i=0;i<=nParts;i++)
    printf ( " %5i %12.7g\n",i,partition[i] );
}

void pisa::DataScan::write ( mmdb::io::File & f )  {
mmdb::realtype b;
int            i;

  f.WriteInt ( &nParts    );
  f.WriteInt ( &scanParam );

  for (i=0;i<=nParts;i++)  {
    b = partition[i];
    f.WriteFloat ( &b );
  }

  for (i=0;i<nParts;i++)
    distSet[i]->write ( f );

}

void pisa::DataScan::read ( mmdb::io::File & f )  {
mmdb::realtype b;
int            i;

  clear();

  f.ReadInt ( &nParts    );
  f.ReadInt ( &scanParam );

  partition = new mmdb::shortreal[nParts+1];

  for (i=0;i<=nParts;i++) {
    f.ReadFloat ( &b );
    partition[i] = b;
  }

  distSet = new PStatDistSet[nParts];
  for (i=0;i<nParts;i++)  {
    distSet[i] = new StatDistSet();
    distSet[i]->read ( f );
  }

}


//  =========================  DataStats  ==========================

pisa::DataStats::DataStats()  {
  dataName  = NULL;
  scan      = NULL;
  paramID   = NULL;
  paramName = NULL;
  procKey   = NULL;
  buf       = NULL;
  nParams   = 0;
}

pisa::DataStats::~DataStats()  {
  clear();
}

void pisa::DataStats::clear()  {
  SDS.clear();
  if (scan)  {
    for (int i=0;i<nParams;i++)  {
      if (scan[i])      delete   scan[i];
      if (paramID[i])   delete[] paramID[i];
      if (paramName[i]) delete[] paramName[i];
    }
    delete[] scan;
    scan  = NULL;
    mmdb::FreeVectorMemory ( paramID  ,0 );
    mmdb::FreeVectorMemory ( paramName,0 );
    mmdb::FreeVectorMemory ( procKey  ,0 );
    mmdb::FreeVectorMemory ( buf      ,0 );
  }
  nParams = 0;
  if (dataName)  {
    delete[] dataName;
    dataName =NULL;
  }
}

void pisa::DataStats::initStats ( int n_parameters )  {

  clear();

  nParams = n_parameters;
  scan    = new PDataScan[nParams];
  mmdb::GetVectorMemory ( paramID  ,nParams,0 );
  mmdb::GetVectorMemory ( paramName,nParams,0 );
  mmdb::GetVectorMemory ( procKey  ,nParams,0 );
  mmdb::GetVectorMemory ( buf      ,nParams,0 );
  for (int i=0;i<nParams;i++)  {
    scan     [i] = new DataScan();
    paramID  [i] = NULL;
    paramName[i] = NULL;
    procKey  [i] = DataStats::Raw;
    buf      [i] = 0.0;
  }

}

void pisa::DataStats::setDataName ( mmdb::cpstr name )  {
  mmdb::CreateCopy ( dataName,name );
}


void pisa::DataStats::setParamID ( int paramNo, mmdb::cpstr ID )  {
  if ((paramNo>=0) && (paramNo<nParams))
    mmdb::CreateCopy ( paramID[paramNo],ID );
}

void pisa::DataStats::setParamName ( int paramNo, mmdb::cpstr name )  {
  if ((paramNo>=0) && (paramNo<nParams))
    mmdb::CreateCopy ( paramName[paramNo],name );
}

void pisa::DataStats::setProcKey ( int paramNo, PROC_KEY key )  {
  if ((paramNo>=0) && (paramNo<nParams))
    procKey[paramNo] = key;
}

void pisa::DataStats::initScan ( int paramNo, mmdb::rvector partition,
                                 int nPartitions )  {

  if (paramNo==0)
    SDS.initSet ( nParams,-1 );

  if ((paramNo<0) || (paramNo>=nParams))
    return;

  if (!scan[paramNo])
    scan[paramNo] = new DataScan();

  scan[paramNo]->initScan ( nPartitions,partition,nParams,paramNo );

}

void pisa::DataStats::convertParams ( mmdb::rvector p )  {
  for (int i=0;i<nParams;i++)
    switch (procKey[i])  {
      default:
      case DataStats::Raw     :  buf[i] =  p[i];       break;
      case DataStats::Minus   :  buf[i] = -p[i];       break;
      case DataStats::MinusLog:  buf[i] = -log(p[i]);
    }
}

void pisa::DataStats::addData ( mmdb::rvector p, bool inAssembly )  {
  convertParams ( p );
  SDS.addData ( buf,inAssembly );
  for (int i=0;i<nParams;i++)  {
    mmdb::RSwap ( p[i],buf[i] );
    scan[i]->addData ( buf,inAssembly );
    mmdb::RSwap ( p[i],buf[i] );
  }
}

void pisa::DataStats::addData ( mmdb::rvector p, bool inAssembly,
                                RStatDistSet D )  {
  convertParams ( p );
  D.addData ( buf,inAssembly );
}

void pisa::DataStats::makeParts ( int maxNParts )  {
  SDS.makeParts ( maxNParts );
  for (int i=0;i<nParams;i++)
    scan[i]->makeParts ( maxNParts );
}

pisa::PDataScan pisa::DataStats::getDataScan ( int paramNo )  {
  if ((0<=paramNo) && (paramNo<nParams))
    return scan[paramNo];
  return NULL;
}

int pisa::DataStats::getStats ( int scanParamNo, mmdb::rvector stats,
                                PStatDistSet D ) {
PStatDistSet distSet;
int          partNo;

  if (scanParamNo==-1)  {
    distSet = &SDS;
  } else if (scanParamNo==-2)  {
    distSet = D;
  } else   {

    if ((scanParamNo<0) || (scanParamNo>=nParams))  {
      for (int i=0;i<nParams;i++)
        stats[i] = 0.0;
      return -1;
    }

    partNo  = scan[scanParamNo]->getPartNo  ( stats[scanParamNo] );
    distSet = scan[scanParamNo]->getDistSet ( partNo );

  }

  if (!distSet)  {
    for (int i=0;i<nParams;i++)
      stats[i] = 0.0;
    return -2;
  }

  convertParams ( stats );
  distSet->getStats ( buf );
  for (int i=0;i<nParams;i++)
    stats[i] = buf[i];

  return 0;

}

bool pisa::DataStats::write ( mmdb::cpstr fileName )  {
mmdb::io::File f;
bool           ok = false;

  f.assign ( fileName,false,true );
  if (f.rewrite())  {
    write ( f );
    f.shut();
    ok = true;
  }

  return ok;

}

bool pisa::DataStats::read ( mmdb::cpstr fileName )  {
mmdb::io::File f;
bool           ok = false;

  f.assign ( fileName,false,true );
  if (f.rewrite())  {
    ok = read ( f );
    f.shut();
  }

  return ok;

}


#define fileLabel  "##DATA:STATS"

void pisa::DataStats::write ( mmdb::io::File & f )  {

  f.WriteFile ( fileLabel,strlen(fileLabel) );

  f.WriteInt    ( &nParams );
  f.CreateWrite ( dataName );
  SDS.write ( f );
  for (int i=0;i<nParams;i++)  {
    f.CreateWrite ( paramID  [i]  );
    f.CreateWrite ( paramName[i]  );
    f.WriteInt    ( &(procKey[i]) );
    scan[i]->write ( f );
  }

}

bool pisa::DataStats::read ( mmdb::io::File & f )  {
char label[20];
int  n_parameters;

  clear();

  f.ReadFile ( label,strlen(fileLabel) );
  if (strncmp(label,fileLabel,strlen(fileLabel)))
    return false;

  f.ReadInt ( &n_parameters );
  initStats ( n_parameters  );

  f.CreateRead ( dataName );
  SDS.read ( f );
  for (int i=0;i<nParams;i++)  {
    f.CreateRead ( paramID  [i] );
    f.CreateRead ( paramName[i] );
    f.ReadInt    ( &(procKey[i]) );
    scan[i]->read ( f );
  }

  return true;

}
