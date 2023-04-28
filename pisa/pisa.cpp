// $Id: pisa.cpp $
// =================================================================
//
//    15.02.19   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// -----------------------------------------------------------------
//
//  **** Module  :  PISA standalone implementation
//       ~~~~~~~~~
//  **** Project :  Protein Interfaces, Surfaces and Assembles PISA
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2019
//
// =================================================================
//


#include <string.h>
#include <stdlib.h>

#include "pisa.h"
#include "pisalib/pisa_lists.h"
#include "pisalib/pisa_detail.h"
#include "pisalib/pisa_view.h"

#include "pisa_analyse.h"

#include<iostream>
using namespace std;

void printInstructions ( mmdb::cpstr pisaName )  {
  printf (
    "\n"
    " Protein Interactions, Surfaces and Assemblies (PISA)\n"
    " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    " v%i.%i.%i built on %s\n"
    " with SSM v%i.%i.%i, SRS v%i.%i.%i, MMDB v%i.%i.%i\n"
    "\n"
    " (C) Eugene Krissinel 2007-2019\n"
    "\n"
    " .============================================================.\n"
    " !   Use of this program is subject to citing the following   !\n"
    " !   publication:                                             !\n"
    " !      E. Krissinel and K. Henrick (2007). Inference of      !\n"
    " !      macromolecular assemblies from crystalline state.     !\n"
    " !      J. Mol. Biol. 372, 774-797.                           !\n"
    " `============================================================'\n"
    "\n"
    " USAGE:\n"
    " ~~~~~~\n"
    "\n"
    " 1. Perform PISA analysis of a structure given in PDB or\n"
    "    mmCIF file 'coorfile':\n"
    "\n"
    " # %s name " job_analyse_tag " coorfile [options] [cfg]\n"
    "\n"
    " where 'name' is a mandatory session name, [cfg] stands\n"
    " for optional configuration file. Configuration file must\n"
    " be specified unless pointed out by environmental variable\n"
    " %s. [options] may be used to specify optional\n"
    " processing keys:\n"
    "\n"
    " " job_lig_exclude_tag "'list'"
                            " - exclude ligands in the given list\n"
    "               when calculating the oligomeric state. The\n"
    "               list should be a comma-separated, no spaces,\n"
    "               case-sensitive set of ligand names. Two\n"
    "               entries have special meaning:\n"
    "                 %s    :"
                                            " exclude all ligands\n"
    "                 %s :"
                                            " exclude ligands"
                                      " mentioned in\n"
    "                           'agents.dat' file (specified in\n"
    "                           configuration file [cfg])\n"
    " " job_lig_auto_tag "  - auto processing for ligands (default);\n"
    "               ligands may be forcefully fixed on protein\n"
    "               surface subject to task complexity and data\n"
    "               available\n"
    " " job_lig_fixed_tag " - fix ligands on protein surface\n"
    " " job_lig_free_tag "  - do not fix ligands on protein surface.\n"
    "\n"
    " Example:\n"
    " %s 3gcb -analyse 3gcb.pdb " job_lig_exclude_tag "'%s,ATP' pisa.cfg\n"
    "\n"
    " After finishing, the results will be stored\n"
    " in session directory identified by session name 'name'.\n"
    " Use the following commands for retrieving the results.\n"
    "\n"
    " 2. Retrieve lists of interfaces, monomers, assemblies and\n"
    "    assembly stock:\n"
    "\n"
   " # %s name " job_list_tag " {" job_interfaces_tag "|" job_monomers_tag
               "|" job_assemblies_tag "|" job_stock_tag "} [cfg]\n"
    "\n"
    " 3. View an interface, monomer, assembly or a dissociate\n"
    "    in Rasmol:\n"
    "\n"
    " # %s name " job_view_tag " spec serial_no [cfg]\n"
    "\n"
    " where spec={" job_interface_tag "|" job_monomer_tag "|" job_assembly_tag
    "|" job_dissociate_tag "}, and\n"
    " serial_no is serial number shown in the corresponding list.\n"
    "\n"
    " 4. View an interface, monomer, assembly or a dissociate\n"
    "    in CCP4-MG:\n"
    "\n"
    " # %s name " job_mg_tag " spec serial_no [cfg]\n"
    "\n"
    " 5. Download an interface, monomer, assembly or a dissociate:\n"
    "\n"
    " # %s name " job_download_tag " spec serial_no [cfg] > output_file\n"
    " # %s name " job_pdb_tag " spec serial_no [cfg] > output_file\n"
    " # %s name " job_cif_tag " spec serial_no [cfg] > output_file\n"
    "\n"
    " Keys " job_download_tag " and " job_pdb_tag " are equivalent. They"
    " produce output in PDB format.\n Key " job_cif_tag " produces output"
    " in mmCIF format.\n"
    "\n"
    " 6. Get details of an interface, monomer or assembly:\n"
    "\n"
    " # %s name " job_detail_tag " spec serial_no [cfg]\n"
    "\n"
    " 7. Generate Remark 350:\n"
    "\n"
    " # %s name " job_350n_tag " assembly_serial_no [cfg]\n"
    "\n"
    " Special use:\n"
    "  assembly_serial_no=0: generate remark for complexe made by\n"
    "                        the content of asu (i.e. the bare content\n"
    "                        of input file), without symmetry analysis)\n"
    "  assembly_serial_no<0: generate remarks 300 and 350 for the whole\n"
    "                        crystal split number -assembly_serial_no.\n"
    "                        E.g., in order to output remarks 300 and\n"
    "                        350 for the most probable split, use\n"
    "                        assembly_serial_no=-1.\n"
    "\n"
    " 8. Remark 350 for automatic annotation:\n"
    "\n"
    " # %s name " job_350_tag " [cfg]\n"
    "\n"
    " If stable assemblies were identified, then this is identical to\n"
    "\n"
    " # %s name " job_350n_tag " -1 [cfg]\n"
    "\n"
    " However if no stable assemblies were found, then the full set of\n"
    " remarks 350 for individual chains is generated.\n"
    "\n"
    " 9. Generate XML output:\n"
    "\n"
    " # %s name " job_xml_tag " {interfaces|assemblies} [cfg] "
    "> outputfile.xml\n"
    "\n"
    " 10. Erase session data:\n"
    "\n"
    " # %s name " job_erase_tag " [cfg]\n"
    "\n"
    " 11. Generate template configuration file:\n"
    "\n"
    " # %s " job_cfg_template_tag " [cfg]\n"
    "\n"
    " 12. Analyse interfaces for input assembly as is, omit calculating assemblies:\n"
    "\n"
    " # %s 3gcb -analyse 3gcb.pdb " job_as_is_on_tag " [cfg]\n"
    "\n"
    "\n",
    pisa::MAJOR_VERSION,pisa::MINOR_VERSION,pisa::MICRO_VERSION,
    pisa::BuildDate,
    ssm::MAJOR_VERSION,ssm::MINOR_VERSION,ssm::MICRO_VERSION,
    ccp4srs::MAJOR_VERSION,ccp4srs::MINOR_VERSION,ccp4srs::MICRO_VERSION,
    mmdb::MAJOR_VERSION,mmdb::MINOR_VERSION,mmdb::MICRO_VERSION,
    pisaName,
      pisa::conf_file,
        pisa::lig_excl_all_key,
        pisa::lig_excl_agents_key,
      pisaName,
        pisa::lig_excl_agents_key,
    pisaName,pisaName,pisaName,pisaName,pisaName,pisaName,pisaName,
    pisaName,pisaName,pisaName,pisaName,pisaName,pisaName,pisaName);

#ifdef _dimer_special
  printf (
    " 10. Special dimer processing:\n"
    "\n"
    " # %s name "job_dimer_special_tag" pdbCode pdbDir outDir [cfg]\n"
    "\n"
    " Here, PISA checks whether 'pdbCode' is dimer with and without\n"
    " ligands. If it is, a line of dimer parameters is output:\n"
    "     'ASA BSA DGi DG DG0 Nhb Nsb Nds PDB'\n"
    " (DG0 is dissociation barrier with ligands), and two files:\n"
    " 1xyz_dimer.pdb (the dimer) and 1xyz_i.pdb (a different\n"
    " interface) are written into outDir.\n"
    "\n"
    " 11. Special dimer loop processing:\n"
    "\n"
    " # %s name "job_dimer_loop_tag" pdbList pdbDir outDir [cfg]\n"
    "\n"
    " Here, PISA reads PDB codes from 'pdbList' and then applies\n"
    " 'special dimer' processing (cf. above) to all of them.\n"
    "\n",
    pisaName,pisaName );
#endif

}


void printLogo()  {
  printf (
  "\n"
  " .==============================================================.\n"
  " !                                                              !\n"
  " !      Protein Interfaces, Surfaces and Assemblies (PISA)      !\n"
  " !                    v%i.%i.%i from %s                    !\n"
  " !           uses SSM v%i.%i.%i, SRS v%i.%i.%i, MMDB v%i.%i.%i"
                                                       "           !\n"
  " !                                                              !\n"
  " +--------------------------------------------------------------+\n"
  " !          Written by Eugene B. Krissinel 2007-2019.           !\n"
  " !               CCP4, Research Complex at Harwell              !\n"
  " !         Rutherford Appleton Laboratory, Didcot, OX11 0FA     !\n"
  " !                        UNITED KINGDOM                        !\n"
  " !              Tel/Fax: ++44 (1235) 567725/7720                !\n"
  " !             E-mail: eugene.krissinel@stfc.ac.uk              !\n"
  " +--------------------------------------------------------------+\n"
  " !  Use of this program is subject to citing the following      !\n"
  " !  publication:                                                !\n"
  " !      E. Krissinel and K. Henrick (2007). Inference of        !\n"
  " !      macromolecular assemblies from crystalline state.       !\n"
  " !      J. Mol. Biol. 372, 774-797.                             !\n"
  " !  Other relevant citations may be found at                    !\n"
  " !       http://www.ebi.ac.uk/msd-srv/prot_int/picite.html      !\n"
  " `=============================================================='\n"
  "\n"
  "\n",
  pisa::MAJOR_VERSION,pisa::MINOR_VERSION,pisa::MICRO_VERSION,
  pisa::BuildDate,
  ssm::MAJOR_VERSION,ssm::MINOR_VERSION,ssm::MICRO_VERSION,
  ccp4srs::MAJOR_VERSION,ccp4srs::MINOR_VERSION,ccp4srs::MICRO_VERSION,
  mmdb::MAJOR_VERSION,mmdb::MINOR_VERSION,mmdb::MICRO_VERSION
   );
}

class InputData {

  public:
    mmdb::pstr       confFile;  // is not allocated, just pointer copied
    mmdb::pstr       ligExcl;
    int              serNo;
    pisa::LIGAND_KEY ligKey;
    pisa::AS_IS_KEY  asisKey;

    InputData ();
    ~InputData();

    void init   ();
    void dispose();
    void print  ();

};

InputData::InputData()  {
  init();
}

InputData::~InputData()  {
  dispose();
}

void InputData::init()  {
  confFile = NULL;
  ligKey   = pisa::LIGANDS_Auto;
  asisKey  = pisa::AS_IS_off;
  ligExcl  = NULL;
  serNo    = 0;
}

void InputData::dispose()  {
  if (ligExcl)   delete[] ligExcl;
  confFile = NULL;
  ligExcl  = NULL;
}

void InputData::print()  {
  if (confFile)  printf ( " confFile = '%s'\n",confFile );
           else  printf ( " confFile = NULL\n" );
  if (ligExcl)   printf ( " ligExcl  = '%s'\n",ligExcl  );
           else  printf ( " ligExcl  = NULL\n" );
  printf ( " serNo   = %i\n",serNo  );
  printf ( " ligKey  = %i\n",ligKey );
  printf ( " asisKey  = %i\n",asisKey );
}


JOB_CONTROL_ID checkArguments ( int argc, char **argv,
                                      InputData & inpData )  {
mmdb::pstr endptr;
int        k;

  inpData.dispose();
  inpData.init   ();

  if (argc<2)
    return JOB_Help;

  if (argc<=3)  {

    if (!strcasecmp(argv[1],job_cfg_template_tag))  {
      if ((argc<2) || (argc>3))  return JOB_Help;
      if (argc==3)  inpData.confFile = argv[2];
      return JOB_Cfg_Template;
    }

  }

  if (argc<3)  {

    if (argc==2)  {
      if (!strcasecmp(argv[1],job_logo_tag))
        return JOB_Logo;
      if (!strcasecmp(argv[1],job_help_tag))
        return JOB_Help;
    }

    return JOB_Help;

  }

  if (!strcasecmp(argv[2],job_analyse_tag))  {
    if (argc<4)  return JOB_Help;
    for (k=4;k<argc;k++)  {
      if (!strcmp(argv[k],job_lig_auto_tag))
        inpData.ligKey = pisa::LIGANDS_Auto;
      else if (!strcmp(argv[k],job_lig_fixed_tag))
        inpData.ligKey = pisa::LIGANDS_FixAll;
      else if (!strcmp(argv[k],job_lig_free_tag))
        inpData.ligKey = pisa::LIGANDS_FreeAll;
      else if (!strncmp(argv[k],job_lig_exclude_tag,
                        strlen(job_lig_exclude_tag)))  {
        mmdb::CreateCopy ( inpData.ligExcl,
                           argv[k]+strlen(job_lig_exclude_tag) );
        mmdb::DelSpaces ( inpData.ligExcl,'\'' );
        mmdb::DelSpaces ( inpData.ligExcl,' '  );
      }
      else if (!strcmp(argv[k],job_as_is_on_tag))
        inpData.asisKey = pisa::AS_IS_on;
      else if (k==argc-1)
        inpData.confFile = argv[k];
      else
        return JOB_Help;
    }
    return JOB_Analyse;
  }

#ifdef _dimer_special
  if (!strcasecmp(argv[2],job_dimer_special_tag))  {
    if ((argc<6) || (argc>7))  return JOB_Help;
    if (argc==7)  inpData.confFile = argv[6];
    return JOB_DimerSpecial;
  }
  if (!strcasecmp(argv[2],job_dimer_loop_tag))  {
    if ((argc<6) || (argc>7))  return JOB_Help;
    if (argc==7)  inpData.confFile = argv[6];
    return JOB_DimerLoop;
  }
#endif

  if (!strcasecmp(argv[2],job_list_tag))  {
    if ((argc<4) || (argc>5))  return JOB_Help;
    if (argc==5)  inpData.confFile = argv[4];
    if (!strncmp(argv[3],job_interfaces_tag,4))
      return JOB_List_Interfaces;
    if (!strncmp(argv[3],job_monomers_tag,4))
      return JOB_List_Monomers;
    if (!strncmp(argv[3],job_assemblies_tag,4))
      return JOB_List_Assemblies;
    if (!strncmp(argv[3],job_stock_tag,4))
      return JOB_List_Stock;
    return JOB_Help;
  }

  if (!strcasecmp(argv[2],job_view_tag))  {
    if ((argc<5) || (argc>6))  return JOB_Help;
    if (argc==6)  inpData.confFile = argv[5];
    inpData.serNo = mmdb::mround ( strtod(argv[4],&endptr) );
    if ((!inpData.serNo) && (endptr==argv[4]))
      return JOB_Help;
    if (!strncmp(argv[3],job_interface_tag,4))
      return JOB_View_Interface;
    if (!strncmp(argv[3],job_monomer_tag,4))
      return JOB_View_Monomer;
    if (!strncmp(argv[3],job_assembly_tag,4))
      return JOB_View_Assembly;
    if (!strncmp(argv[3],job_dissociate_tag,4))
      return JOB_View_Dissociate;
    return JOB_Help;
  }

  if (!strcasecmp(argv[2],job_mg_tag))  {
    if ((argc<5) || (argc>6))  return JOB_Help;
    if (argc==6)  inpData.confFile = argv[5];
    inpData.serNo = mmdb::mround ( strtod(argv[4],&endptr) );
    if ((!inpData.serNo) && (endptr==argv[4]))
      return JOB_Help;
    if (!strncmp(argv[3],job_interface_tag,4))
      return JOB_MG_Interface;
    if (!strncmp(argv[3],job_monomer_tag,4))
      return JOB_MG_Monomer;
    if (!strncmp(argv[3],job_assembly_tag,4))
      return JOB_MG_Assembly;
    if (!strncmp(argv[3],job_dissociate_tag,4))
      return JOB_MG_Dissociate;
    return JOB_Help;
  }

  if ((!strcasecmp(argv[2],job_download_tag)) ||
      (!strcasecmp(argv[2],job_pdb_tag))      ||
      (!strcasecmp(argv[2],job_cif_tag)))  {
    if ((argc<5) || (argc>6))  return JOB_Help;
    if (argc==6)  inpData.confFile = argv[5];
    inpData.serNo = mmdb::mround ( strtod(argv[4],&endptr) );
    if ((!inpData.serNo) && (endptr==argv[4]))
      return JOB_Help;
    if (strcasecmp(argv[2],job_cif_tag))  {
      if (!strncmp(argv[3],job_interface_tag,4))
        return JOB_Download_Interface;
      if (!strncmp(argv[3],job_monomer_tag,4))
        return JOB_Download_Monomer;
      if (!strncmp(argv[3],job_assembly_tag,4))
        return JOB_Download_Assembly;
      if (!strncmp(argv[3],job_dissociate_tag,4))
        return JOB_Download_Dissociate;
    } else  {
      if (!strncmp(argv[3],job_interface_tag,4))
        return JOB_CIF_Interface;
      if (!strncmp(argv[3],job_monomer_tag,4))
        return JOB_CIF_Monomer;
      if (!strncmp(argv[3],job_assembly_tag,4))
        return JOB_CIF_Assembly;
      if (!strncmp(argv[3],job_dissociate_tag,4))
        return JOB_CIF_Dissociate;
    }
    return JOB_Help;
  }

  if (!strcasecmp(argv[2],job_detail_tag))  {
    if ((argc<5) || (argc>6))  return JOB_Help;
    if (argc==6)  inpData.confFile = argv[5];
    inpData.serNo = mmdb::mround ( strtod(argv[4],&endptr) );
    if ((!inpData.serNo) && (endptr==argv[4]))
      return JOB_Help;
    if (!strncmp(argv[3],job_interface_tag,4))
      return JOB_Detail_Interface;
    if (!strncmp(argv[3],job_monomer_tag,4))
      return JOB_Detail_Monomer;
    if (!strncmp(argv[3],job_assembly_tag,4))
      return JOB_Detail_Assembly;
    return JOB_Help;
  }

  if (!strcasecmp(argv[2],job_350n_tag))  {
    if ((argc<4) || (argc>5))  return JOB_Help;
    if (argc==5)  inpData.confFile = argv[4];
    inpData.serNo = mmdb::mround ( strtod(argv[3],&endptr) );
    if ((!inpData.serNo) && (endptr==argv[3]))
      return JOB_Help;
    return JOB_Remark_350n;
  }

  if (!strcasecmp(argv[2],job_350_tag))  {
    if ((argc<3) || (argc>4))  return JOB_Help;
    if (argc==4)  inpData.confFile = argv[3];
    inpData.serNo = mmdb::MaxInt4;
    return JOB_Remark_350;
  }

  if (!strcasecmp(argv[2],job_xml_tag))  {
    if ((argc<4))  return JOB_Help;
    //if ((argc<4) || (argc>5))  return JOB_Help;
    for (k=4;k<argc;k++)  {
      if (!strcmp(argv[k],job_as_is_on_tag))
        inpData.asisKey = pisa::AS_IS_on;
      else if (k==argc-1)
        inpData.confFile = argv[k];
    }
    //if (argc==5)  inpData.confFile = argv[4];
    if (!strncmp(argv[3],job_interfaces_tag,4))
      return JOB_Make_XML_Interfaces;
    if (!strncmp(argv[3],job_assemblies_tag,4))
      return JOB_Make_XML_Assemblies;
  }

  if (!strcasecmp(argv[2],job_erase_tag))  {
    if ((argc<3) || (argc>4))  return JOB_Help;
    if (argc==4)  inpData.confFile = argv[3];
    return JOB_Erase;
  }

  return JOB_Help;

}

#ifdef _dimer_special
void  GetPDBList ( mmdb::cpstr FName, PIDCode & pdbCode, int & nCodes,
                   int & iStart )  {
CFile f;
char  S[1000];
pstr  p;

  pdbCode = NULL;
  nCodes  = 0;

  f.assign ( FName,true );
  if (!f.reset(true))  return;

  f.ReadLine ( S,sizeof(S) );
  iStart = atoi(S);

  pdbCode = new IDCode[40000];
  while (!f.FileEnd())  {
    f.ReadLine ( S,sizeof(S) );
    p = S;
    while (*p==' ')  p++;
    while (*p && (*p!=' '))  p++;
    while (*p==' ')  p++;
    strncpy ( pdbCode[nCodes],p,4 );
    pdbCode[nCodes][4] = char(0);
    nCodes++;
  }

}
#endif

int main ( int argc, char ** argv, char ** env )  {
UNUSED_ARGUMENT(env);
pisa::PData         pisaData;
pisa::PAnalyse      pisaAnalyse;
pisa::PLists        pisaLists;
pisa::PView         pisaView;
pisa::PDetail       pisaDetail;
InputData           inpData;
mmdb::pstr          S = NULL;
pisa::VIEWER_TARGET target;
pisa::RESULT_CODE   rc;

  mmdb::InitMatType();
  ssm::InitGraph();

  target = pisa::TARGET_NONE;

  rc = pisa::RESULT_Ok;

  switch (checkArguments(argc,argv,inpData)) {

    case JOB_Analyse :
//        inpData.print();
        pisaAnalyse = new pisa::Analyse  ( inpData.confFile );
        pisaAnalyse->setLigKey           ( inpData.ligKey   );
	pisaAnalyse->setAsIsKey          ( inpData.asisKey  );
        pisaAnalyse->setLigExclude       ( inpData.ligExcl  );
        rc = pisaAnalyse->analyse        ( argv[1],argv[3]  );
        pisaAnalyse->printResultMessage  ( rc );
        delete pisaAnalyse;
	pisaDetail = new pisa::Detail       ( inpData.confFile );
        pisaDetail->set_As_Is_Key_xml       ( inpData.asisKey  );
	delete pisaDetail;
      break;

#ifdef _dimer_special
    case JOB_DimerSpecial :
        pisaAnalyse = new pisa::Analyse ( confFile );
        rc = pisaAnalyse->dimerSpecial  ( argv[1],argv[3],argv[4],
                                          serNo );
        pisaAnalyse->printResultMessage ( rc );
        delete pisaAnalyse;
        if (rc==RESULT_Ok)  {
          pisaView = new pisa::View ( confFile );
          confFile = NULL;
          if (argv[5][strlen(argv[5])-1]=='/')  {
            target = 1;
            CreateCopCat ( confFile,argv[5],argv[3],"_dimer.pdb" );
          } else  {
            target = 0;
            CreateCopCat ( confFile,argv[5],"/",argv[3],"_dimer.pdb" );
          }
          rc = pisaView->ViewAssembly ( argv[1],1,false,
                                        pisa::TARGET_Download,confFile );
          if (rc==RESULT_Ok)  {
            if (target)
                 CreateCopCat ( confFile,argv[5],argv[3],"_i.pdb" );
            else CreateCopCat ( confFile,argv[5],"/",argv[3],"_i.pdb");
            rc = pisaView->ViewInterface ( argv[1],serNo,
                                           pisa::TARGET_Download,confFile );
          }
          if (confFile)  delete[] confFile;
          pisaView->printResultMessage ( rc );
          delete pisaView;
        }
      break;

    case pisa::JOB_DimerLoop :
        PIDCode  pdbCode;
        pstr     FName;
        int      nCodes,i,iStart;
        GetPDBList ( argv[3],pdbCode,nCodes,iStart );
        for (i=0;i<nCodes;i++)  {
          pisaAnalyse = new pisa::Analyse  ( inpData.confFile );
          inpData.serNo = -iStart;
          printf ( "\n ==== %4s\n",pdbCode[i] );
          rc = pisaAnalyse->DimerSpecial  ( argv[1],pdbCode[i],argv[4],
                                            inpData.serNo );
          pisaAnalyse->printResultMessage ( rc );
          delete pisaAnalyse;
          if (rc==RESULT_Ok)  {
            iStart++;
            pisaView = new pisa::View ( confFile );
            FName    = NULL;
            if (argv[5][strlen(argv[5])-1]=='/')  {
              target = 1;
              CreateCopCat ( FName,argv[5],pdbCode[i],"_dimer.pdb.gz" );
            } else  {
              target = 0;
              CreateCopCat ( FName,argv[5],"/",pdbCode[i],
                                                      "_dimer.pdb.gz" );
            }
            rc = pisaView->ViewAssembly ( argv[1],1,false,
                                          pisa::TARGET_Download,FName );
            if (rc==RESULT_Ok)  {
              if (target)
                   CreateCopCat ( FName,argv[5],pdbCode[i],"_i.pdb.gz" );
              else CreateCopCat ( FName,argv[5],"/",pdbCode[i],
                                  "_i.pdb.gz");
              rc = pisaView->ViewInterface ( argv[1],serNo,
                                             pisa::TARGET_Download,FName );
            }
            if (FName)  delete[] FName;
            pisaView->printResultMessage ( rc );
            delete pisaView;
          }
        }
        if (pdbCode)  delete[] pdbCode;
      break;
#endif

    case JOB_List_Interfaces :
        pisaLists = new pisa::Lists ( inpData.confFile );
        rc = pisaLists->ListInterfaces ( argv[1],"stdout" );
        pisaLists->printResultMessage ( rc );
        delete pisaLists;
      break;

    case JOB_List_Monomers :
        pisaLists = new pisa::Lists ( inpData.confFile );
        rc = pisaLists->ListMonomers ( argv[1],"stdout" );
        pisaLists->printResultMessage ( rc );
        delete pisaLists;
      break;

    case JOB_List_Assemblies :
        pisaLists = new pisa::Lists ( inpData.confFile );
        rc = pisaLists->ListAssemblies ( argv[1],"stdout" );
        pisaLists->printResultMessage ( rc );
        delete pisaLists;
      break;

    case JOB_List_Stock :
        pisaLists = new pisa::Lists ( inpData.confFile );
        rc = pisaLists->ListStock ( argv[1],"stdout" );
        pisaLists->printResultMessage ( rc );
        delete pisaLists;
      break;

    case JOB_Download_Interface :
        target = pisa::TARGET_Download;
    case JOB_CIF_Interface      :
        if (target<0) target = pisa::TARGET_CIF;
    case JOB_MG_Interface       :
        if (target<0) target = pisa::TARGET_CCP4MG;
    case JOB_View_Interface     :
        if (target<0) target = pisa::TARGET_Rasmol;
        pisaView = new pisa::View ( inpData.confFile );
        rc = pisaView->ViewInterface ( argv[1],inpData.serNo,
                                       target,"stdout" );
        pisaView->printResultMessage ( rc );
        delete pisaView;
      break;

    case JOB_Download_Monomer :
        target = pisa::TARGET_Download;
    case JOB_CIF_Monomer       :
        if (target<0) target = pisa::TARGET_CIF;
    case JOB_MG_Monomer       :
        if (target<0) target = pisa::TARGET_CCP4MG;
    case JOB_View_Monomer     :
        if (target<0) target = pisa::TARGET_Rasmol;
        pisaView = new pisa::View ( inpData.confFile );
        rc = pisaView->ViewMonomer ( argv[1],inpData.serNo,
                                     target,"stdout" );
        pisaView->printResultMessage ( rc );
        delete pisaView;
      break;

    case JOB_Download_Assembly :
        target = pisa::TARGET_Download;
    case JOB_CIF_Assembly       :
        if (target<0) target = pisa::TARGET_CIF;
    case JOB_MG_Assembly       :
        if (target<0) target = pisa::TARGET_CCP4MG;
    case JOB_View_Assembly     :
        if (target<0) target = pisa::TARGET_Rasmol;
        pisaView = new pisa::View ( inpData.confFile );
        rc = pisaView->ViewAssembly  ( argv[1],inpData.serNo,false,
                                       target,"stdout" );
        pisaView->printResultMessage ( rc );
        delete pisaView;
      break;

    case JOB_Download_Dissociate :
        target = pisa::TARGET_Download;
    case JOB_CIF_Dissociate       :
        if (target<0) target = pisa::TARGET_CIF;
    case JOB_MG_Dissociate       :
        if (target<0) target = pisa::TARGET_CCP4MG;
    case JOB_View_Dissociate     :
        if (target<0) target = pisa::TARGET_Rasmol;
        pisaView = new pisa::View ( inpData.confFile );
        rc = pisaView->ViewAssembly  ( argv[1],inpData.serNo,true,
                                       target,"stdout" );
        pisaView->printResultMessage ( rc );
        delete pisaView;
      break;

    case JOB_Detail_Interface    :
        pisaDetail = new pisa::Detail ( inpData.confFile );
        rc = pisaDetail->DetailInterface ( argv[1],"stdout",
                                           inpData.serNo );
        pisaDetail->printResultMessage   ( rc );
        delete pisaDetail;
      break;

    case JOB_Detail_Monomer      :
        pisaDetail = new pisa::Detail ( inpData.confFile );
        rc = pisaDetail->DetailMonomer ( argv[1],"stdout",
                                         inpData.serNo );
        pisaDetail->printResultMessage ( rc );
        delete pisaDetail;
      break;

    case JOB_Detail_Assembly     :
        pisaDetail = new pisa::Detail ( inpData.confFile );
        rc = pisaDetail->DetailAssembly ( argv[1],"stdout",
                                          inpData.serNo );
        pisaDetail->printResultMessage  ( rc );
        delete pisaDetail;
      break;

    case JOB_Remark_350n         :
    case JOB_Remark_350          :
        pisaDetail = new pisa::Detail ( inpData.confFile );
        rc = pisaDetail->Remark350 ( argv[1],"stdout",
                                     inpData.serNo );
        pisaDetail->printResultMessage  ( rc );
        delete pisaDetail;
      break;

    case JOB_Make_XML_Interfaces :
        pisaDetail = new pisa::Detail       ( inpData.confFile );
	pisaDetail->set_As_Is_Key_xml       ( inpData.asisKey  );
        rc = pisaDetail->MakeInterfacesXML ( argv[1],"stdout" );
        pisaDetail->printResultMessage     ( rc );
        delete pisaDetail;
      break;

    case JOB_Make_XML_Assemblies :
        pisaDetail = new pisa::Detail       ( inpData.confFile );
	pisaDetail->set_As_Is_Key_xml       ( inpData.asisKey  );
        rc = pisaDetail->MakeAssembliesXML ( argv[1],"stdout" );
        pisaDetail->printResultMessage     ( rc );
        delete pisaDetail;
      break;

    case JOB_Logo :
        printLogo();
        rc = pisa::RESULT_Logo;
      break;

    case JOB_Erase :
        pisaData = new pisa::Data ( inpData.confFile );
        rc = pisaData->eraseSession ( argv[1] );
        pisaData->printResultMessage ( rc );
        rc = pisa::RESULT_Ok;
        delete pisaData;
      break;

    case JOB_Cfg_Template :
        pisaData = new pisa::Data ( inpData.confFile );
        pisaData->printResultMessage ( rc );
        S = NULL;
        printf ( "%s",pisaData->getConfigurationTemplate(S) );
        if (S)  delete[] S;
        rc = pisa::RESULT_Ok;
        delete pisaData;
      break;

    default : printInstructions ( argv[0] );
              rc = pisa::RESULT_Instructions;

  }

  return rc;

}
