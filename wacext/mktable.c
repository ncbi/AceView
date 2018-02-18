#include "dna.h"

typedef struct enzmoins_struct {
  char *name;
  char *seq;
  int length;
  int overhang;
}ENZMOINS;

#include "../wacext/enzlist.h"

static  ENZMOINS enzListe[]={
	{"AatI","AGGCCT",6,6},
	{"AatII","GACGTC",6,2},
	{"Acc65I","GGTACC",6,1},
	{"AccI","GTmkAC",6,6},
	{"AccII","CGCG",4,6},
	{"AccIII","TCCGGA",6,6},
	{"AciI","CCGC",4,6},
	{"AcyI","GrCGyC",6,6},
	{"AfaI","GTAC",4,6},
	{"AflII","CTTAAG",6,6},
	{"AflIII","ACryGT",6,6},
	{"AgeI","ACCGGT",6,1},
	{"AhaII","GrCGyC",6,6},
	{"AhaIII","TTTAAA",6,6},
	{"AluI","AGCT",4,6},
	{"Alw21I","GwGCwC",6,6},
	{"Alw26I","GTCTCnnnnn",10,6},
	{"Alw44I","GTGCAC",6,6},
	{"AlwI","GGATCnnnnn",10,6},
	{"AlwNI","CAGnnnCTG",9,6},
	{"AocI","CCTnAGG",7,6},
	{"AosI","TGCGCA",6,6},
	/* {"ApaBI","GCAnnnnnTGC",11,6}, */
	{"ApaI","GGGCCC",6,2},
	{"ApaLI","GTGCAC",6,6},
	{"ApyI","CCwGG",5,6},
	{"AscI","GGCGCGCC",8,1},
	{"AseI","ATTAAT",6,6},
	{"AsnI","ATTAAT",6,6},
	{"Asp700I","GAAnnnnTTC",10,6}, 
	{"Asp718I","GGTACC",6,6},
	{"AspHI","GwGCwC",6,6},
	{"AspI","GACnnnGTC",9,6},
	{"AsuI","GGnCC",5,6},
	{"AsuII","TTCGAA",6,6},
	{"AvaI","CyCGrG",6,6},
	{"AvaII","GGwCC",5,1},
	{"AvaIII","ATGCAT",6,6},
	{"AviII","TGCGCA",6,6},
	{"AvrII","CCTAGG",6,6},
	{"AxyI","CCTnAGG",7,6},
	/* {"BaeI","ACnnnnGTAyC",11,6}, */
	{"BalI","TGGCCA",6,6},
	{"BamHI","GGATCC",6,1},
	{"BanI","GGyrCC",6,6},
	{"BanII","GrGCyC",6,6},
	{"BanIII","ATCGAT",6,6},
	{"BbeI","GGCGCC",6,6},
	{"BbiII","GrCGyC",6,6},
	{"BbrPI","CACGTG",6,6},
	{"BbsI","GAAGAC",6,6},
	{"BbuI","GCATGC",6,6},
	{"BbvI","GCAGCnnnnnnnnnnnn",17,6},
	{"BbvII","GAAGACnnnnnn",12,6},
	{"BccI","CCATC",5,6},
	{"BcefI","ACGGCnnnnnnnnnnnnn",18,6},
	/*	{"BcgI","GCAnnnnnnTCGnnnnnnnnnnnn",24,6}, */
	{"BclI","TGATCA",6,6},
	{"BcnI","CCsGG",5,6},
	{"BetI","wCCGGw",6,6},
	{"BfrI","CTTAAG",6,6},
	/* {"BglI","GCCnnnnnGGC",11,6}, */
	{"BglII","AGATCT",6,3},
	{"BinI","GGATCnnnnn",10,6},
	{"BlnI","CCTAGG",6,6},
	{"Bme18I","GGwCC",5,6},
	{"BmyI","GdGChC",6,6},
	{"Bpu10I","CCTnAGC",7,6},
	{"Bpu14I","TTCGAA",6,6},
	{"Bpu1102I","GCTnAGC",7,6},
	{"BsaAI","yACGTr",6,6},
	{"BsaBI","GATnnnnATC",10,6},
	{"BsaHI","GrCGyC",6,6},
	{"BsaI","GGTCTCnnnnn",11,6},
	{"BsaJI","CCnnGG",6,6},
	{"BscI","ATCGAT",6,6},
	{"Bse21I","CCTnAGG",7,6},
	{"BsePI","GCGCGC",6,6},
	{"BsgI","GTGCAGnnnnnnnnnnnnnnnn",22,6},
	{"BsiI","CTCGTG",6,6},
	{"BsiWI","CGTACG",6,2},
	/*	{"BsiYI","CCnnnnnnnGG",11,6}, */
	{"BslI","CCnnnnnnnGG",11,6},
	{"BsmAI","GTCTCnnnnn",10,6},
	{"BsmI","GAATGCn",7,6},
	{"Bsp50I","CGCG",4,6},
	{"Bsp68I","TCGCGA",6,6},
	{"Bsp106I","ATCGAT",6,6},
	{"Bsp119I","TTCGAA",6,6},
	{"Bsp120I","GGGCCC",6,6},
	{"Bsp1286I","GdGChC",6,6},
	{"BspAI","GATC",4,6},
	{"BspCI","CGATCG",6,6},
	{"BspDI","ATCGAT",6,6},
	{"BspEI","TCCGGA",6,2},
	{"BspGI","CTGGAC",6,6},
	{"BspHI","TCATGA",6,6},
	{"BspMI","ACCTGCnnnnnnnn",14,6},
	{"BspMII","TCCGGA",6,6},
	{"BsrI","ACTGGn",6,6},
	{"BssCI","GGCC",4,6},
	{"BssHII","GCGCGC",6,2},
	{"BssT1I","CCwwGG",6,6},
	{"BstBI","TTCGAA",6,6},
	{"BstEII","GGTnACC",7,6},
	{"BstI","GGATCC",6,6},
	{"BstNI","CCwGG",5,6},
	{"BstPI","GGTnACC",7,6},
	{"BstUI","CGCG",4,6},
	{"BstVI","CTCGAG",6,6},
	/* {"BstXI","CCAnnnnnnTGG",12,6}, */
	{"BstYI","rGATCy",6,6},
	{"Bsu15I","ATCGAT",6,6},
	{"Bsu36I","CCTnAGG",7,6},
	{"BsuRI","GGCC",4,6},
	{"CauII","CCsGG",5,6},
	{"CcrI","CTCGAG",6,6},
	{"CelII","GCTnAGC",7,6},
	{"CfoI","GCGC",4,6},
	{"Cfr9I","CCCGGG",6,6},
	{"Cfr10I","rCCGGy",6,6},
	{"Cfr13I","GGnCC",5,6},
	{"Cfr42I","CCGCGG",6,6},
	{"CfrI","yGGCCr",6,6},
	{"ClaI","ATCGAT",6,6},
	{"CpoI","CGGwCCG",7,6},
	{"Csp6I","GTAC",4,6},
	{"Csp45I","TTCGAA",6,6},
	{"CspI","CGGwCCG",7,6},
	{"CviJI","rGCy",4,6},
	{"CviRI","TGCA",4,6},
	{"CvnI","CCTnAGG",7,6},
	{"DdeI","CTnAG",5,6},
	{"DpnI","GATC",4,6},
	{"DpnII","GATC",4,6},
	{"DraI","TTTAAA",6,6},
	{"DraII","rGGnCCy",7,6},
	{"DraIII","CACnnnGTG",9,6},
	/* {"DrdI","GACnnnnnnGTC",12,6}, */
	{"DrdII","GAACCA",6,6},
	{"DsaI","CCryGG",6,6},
	{"DsaV","CCnGG",5,6},
	{"EaeI","yGGCCr",6,6},
	{"EagI","CGGCCG",6,2},
	{"Eam1104I","CTCTTCnnnn",10,6},
	/* 	{"Eam1105I","GACnnnnnGTC",11,6}, */
	{"EarI","CTCTTCnnnn",10,6},
	{"EciI","TCCGCC",6,6},
	{"Ecl136II","GAGCTC",6,6},
	{"EclXI","CGGCCG",6,6},
	{"Eco24I","GrGCyC",6,6},
	{"Eco31I","GGTCTCnnnnn",11,6},
	{"Eco32I","GATATC",6,6},
	{"Eco47I","GGwCC",5,6},
	{"Eco47III","AGCGCT",6,6},
	{"Eco52I","CGGCCG",6,6},
	{"Eco57I","CTGAAGnnnnnnnnnnnnnnnn",22,6},
	{"Eco64I","GGyrCC",6,6},
	{"Eco72I","CACGTG",6,6},
	{"Eco81I","CCTnAGG",7,6},
	{"Eco88I","CyCGrG",6,6},
	{"Eco91I","GGTnACC",7,6},
	{"Eco105I","TACGTA",6,6},
	{"Eco130I","CCwwGG",6,6},
	{"Eco147I","AGGCCT",6,6},
	/* {"EcoNI","CCTnnnnnAGG",11,6}, */
	{"EcoO65I","GGTnACC",7,6},
	{"EcoO109I","rGGnCCy",7,6},
	{"EcoRI","GAATTC",6,1},
	{"EcoRII","CCwGG",5,6},
	{"EcoRV","GATATC",6,1},
	{"EcoT14I","CCwwGG",6,6},
	{"EcoT22I","ATGCAT",6,6},
	{"EheI","GGCGCC",6,6},
	{"Esp3I","CGTCTCnnnnn",11,6},
	{"EspI","GCTnAGC",7,6},
	{"FauI","CCCGCnnnnnn",11,6},
	{"FbaI","TGATCA",6,6},
	{"FdiII","TGCGCA",6,6},
	{"FinI","GTCCC",5,6},
	{"Fnu4HI","GCnGC",5,6},
	{"FnuDII","CGCG",4,6},
	{"FokI","GGATGnnnnnnnnnnnnn",18,6},
	{"FseI","GGCCGGCC",8,6},
	{"FsiI","rAATTy",6,6},
	{"FspI","TGCGCA",6,6},
	{"GdiII","yGGCCG",6,6},
	{"GsuI","CTGGAGnnnnnnnnnnnnnnnn",22,6},
	{"HaeI","wGGCCw",6,6},
	{"HaeII","rGCGCy",6,6},
	{"HaeIII","GGCC",4,6},
	{"HapII","CCGG",4,6},
	{"HgaI","GACGCnnnnnnnnnn",15,6},
	{"HgiAI","GwGCwC",6,6},
	{"HgiCI","GGyrCC",6,6},
	/* {"HgiEII","ACCnnnnnnGGT",12,6}, */
	{"HgiJII","GrGCyC",6,6},
	{"HhaI","GCGC",4,6},
	{"Hin1I","GrCGyC",6,6},
	{"Hin6I","GCGC",4,6},
	{"HincII","GTyrAC",6,6},
	{"HindII","GTyrAC",6,6},
	{"HindIII","AAGCTT",6,3},
	{"HinfI","GAnTC",5,6},
	{"HinP1I","GCGC",4,6},
	{"HpaI","GTTAAC",6,6},
	{"HpaII","CCGG",4,6},
	{"HphI","GGTGAnnnnnnnn",13,6},
	{"KasI","GGCGCC",6,2},
	{"Kpn2I","TCCGGA",6,6},
	{"Kpn378I","CCGCGG",6,6},
	{"KpnI","GGTACC",6,2},
	{"Ksp632I","CTCTTCnnnn",10,6},
	{"KspI","CCGCGG",6,6},
	{"Kzo9I","GATC",4,6},
	{"LspI","TTCGAA",6,6},
	{"MaeI","CTAG",4,6},
	{"MaeII","ACGT",4,6},
	{"MaeIII","GTnAC",5,6},
	{"MamI","GATnnnnATC",10,6}, 
	{"MboI","GATC",4,6},
	{"MboII","GAAGAnnnnnnnn",13,6},
	{"McrI","CGryCG",6,6},
	{"MfeI","CAATTG",6,6},
	{"MflI","rGATCy",6,6},
	{"MluI","ACGCGT",6,2},
	{"MlyI","GACTCnnnnn",10,6},
	/* {"MmeI","TCCrACnnnnnnnnnnnnnnnnnnnn",26,6}, */
	{"MnlI","CCTCnnnnnnn",11,6},
	{"MraI","CCGCGG",6,6},
	{"MroI","TCCGGA",6,6},
	{"MscI","TGGCCA",6,6},
	{"MseI","TTAA",4,6},
	{"MspI","CCGG",4,6},
	{"MstI","TGCGCA",6,6},
	{"MstII","CCTnAGG",7,6},
	{"MvaI","CCwGG",5,6},
	{"MvnI","CGCG",4,2},
	/*	{"MwoI","GCnnnnnnnGC",11,6}, */
	{"NaeI","GCCGGC",6,6},
	{"NarI","GGCGCC",6,6},
	{"NciI","CCsGG",5,6},
	{"NcoI","CCATGG",6,2},
	{"NdeI","CATATG",6,6},
	{"NdeII","GATC",4,6},
	{"NheI","GCTAGC",6,1},
	{"NlaIII","CATG",4,6},
	{"NlaIV","GGnnCC",6,6},
	{"NotI","GCGGCCGC",8,4},
	{"NruI","TCGCGA",6,6},
	{"NsiI","ATGCAT",6,3},
	{"NspBII","CmGCkG",6,6},
	{"NspHI","rCATGy",6,6},
	{"NspHII","GGwCC",5,6},
	{"NspI","rCATGy",6,6},
	{"NspII","GdGChC",6,6},
	{"NspIII","CyCGrG",6,6},
	{"NspIV","GGnCC",5,6},
	{"NspV","TTCGAA",6,6},
	{"NunII","GGCGCC",6,6},
	{"PacI","TTAATTAA",8,1},
	{"PaeI","GCATGC",6,6},
	{"PaeR7I","CTCGAG",6,6},
	{"PalI","GGCC",4,6},
	/* {"PflMI","CCAnnnnnTGG",11,6}, */
	{"PleI","GAGTCnnnnn",10,6},
	{"PmaCI","CACGTG",6,6},
	{"Pme55I","AGGCCT",6,6},
	{"PmlI","CACGTG",6,6},
	{"PpuMI","rGGwCCy",7,6},
	/*	{"PshAI","GACnnnnGTC",10,6}, */
	{"PssI","rGGnCCy",7,6},
	{"PstI","CTGCAG",6,3},
	{"PvuI","CGATCG",6,6},
	{"PvuII","CAGCTG",6,6},
	{"RleAI","CCCACAnnnnnnnnnnnn",18,6},
	{"RmaI","CTAG",4,6},
	{"RsaI","GTAC",4,6},
	{"RspXI","TCATGA",6,6},
	{"RsrII","CGGwCCG",7,6},
	{"SacI","GAGCTC",6,3},
	{"SacII","CCGCGG",6,6},
	{"SalI","GTCGAC",6,3},
	{"SapI","GCTCTTCnnnn",11,6},
	{"Sau3AI","GATC",4,6},
	{"Sau96I","GGnCC",5,6},
	{"SauI","CCTnAGG",7,6},
	{"ScaI","AGTACT",6,6},
	{"ScrFI","CCnGG",5,6},
	{"SduI","GdGChC",6,6},
	{"SecI","CCnnGG",6,6},
	{"SfaNI","GCATCnnnnnnnnn",14,6},
	{"SfcI","CTyrAG",6,6},
	{"SfeI","CTryAG",6,6},
	/*	{"SfiI","GGCCnnnnnGGCC",13,6}, */
	{"Sfr303I","CCGCGG",6,6},
	{"SfuI","TTCGAA",6,6},
	{"SgrAI","CrCCGGyG",8,6},
	{"SinI","GGwCC",5,6},
	{"SlaI","CTCGAG",6,6},
	{"SmaI","CCCGGG",6,6},
	{"SnaBI","TACGTA",6,6},
	{"SnaI","GTATAC",6,6},
	{"SnoI","GTGCAC",6,6},
	{"SpeI","ACTAGT",6,2},
	{"SphI","GCATGC",6,2},
	{"SplI","CGTACG",6,6},
	{"SpoI","TCGCGA",6,6},
	{"Sse8387I","CCTGCAGG",8,6},
	{"SspI","AATATT",6,6},
	{"SstI","GAGCTC",6,6},
	{"SstII","CCGCGG",6,6},
	{"StuI","AGGCCT",6,6},
	{"StyI","CCwwGG",6,6},
	{"SwaI","ATTTAAAT",8,6},
	{"TaqI","TCGA",4,6},
	{"1","GACCGAnnnnnnnnnnn",17,6},
	{"2","CACCCAnnnnnnnnnnn",17,6},
	{"TfiI","GAwTC",5,6},
	{"ThaI","CGCG",4,6},
	{"Tsp45I","GTsAC",5,6},
	{"TspAI","CCwGG",5,6},
	{"TspEI","AATT",4,6},
	{"Tth111I","GACnnnGTC",9,6},
	{"Tth111II","CAArCAnnnnnnnnnnn",17,6},
	{"TthHB8I","TCGA",4,6},
	{"Uba1108I","TCGTAG",6,6},
	/* {"Van91I","CCAnnnnnTGG",11,6},  */
	{"VneI","GTGCAC",6,6},
	{"VspI","ATTAAT",6,6},
	{"XbaI","TCTAGA",6,1},
	/* {"XcmI","CCAnnnnnnnnnTGG",15,6}, */
	{"XhoI","CTCGAG",6,1},
	{"XhoII","rGATCy",6,6},
	{"XmaI","CCCGGG",6,2},
	{"XmaIII","CGGCCG",6,6},
	{0,0,0,0}
  };


/*********************************************************/

static int enzOrder (const void *a, const void *b)
{
  const ENZ *za = (const ENZ*) a, *zb = (const ENZ*) b ;
  int n, bza=0, bzb=0 ;
  char *cp ;

  /* bizare en queue */
  cp = za->seq -  1 ;
  while (!bza && *++cp)
    switch (*cp)
      {
      case 'A': 
      case 'T': 
      case 'G': 
      case 'C': 
	break ;
      default:
	bza = 1000 ;
	break ;
      }
  cp = zb->seq -  1 ;
  while (!bzb && *++cp)
    switch (*cp)
      {
      case 'A': 
      case 'T': 
      case 'G': 
      case 'C': 
	break ;
      default:
	bzb = 1000 ;
	break ;
      }

  n = bza - bzb ;
  if (n) return n ;

  /* court en tete */
  n = za->length - zb->length ;
  if (n) return n ;

  /* ordre alpha des sites */
  n = strcmp(za->seq, zb->seq) ;
   if (n) return n ;

   /* ordre alpha des moms */
   n = lexstrcmp(za->name, zb->name) ;

   return n ;
}

/*********************************************************/


Array enzInit (void) 
{
  char *dna;
  int j = 0 ;
  static Array enzArray = 0 ;

  ENZ *enz ;
  ENZMOINS *current ;

 if (!enzArray)
   {
     enzArray = arrayCreate (255, ENZ) ;  
     for (current=enzListe; current->name; current++)
       {
	 dna = strnew(current->seq,0);
	 dnaEncodeString(dna);
	 enz = arrayp (enzArray, j++, ENZ) ;
	 
	 enz->name = current->name ;
	 enz->seq = current->seq ;
	 enz->dna = dna ;
	 enz->length = strlen (enz->seq) ;
	 enz->overhang = current -> overhang ;
       }
     
     arraySort (enzArray, enzOrder) ;
   }
  return enzArray ;
}

/*********************************************************/
/*********************************************************/
/*

void main (void) 
{
  Array enzArray = enzInit () ;
  int j  ;
  ENZ *enz ;

  for (j=0 ; j < arrayMax (enzArray) ; j++)
    {
      enz = arrayp (enzArray, j, ENZ) ;
      printf ("{%s,%s,%d,6}\n",enz->name,enz->seq,enz->length);
    }
  printf ("done\n") ;  
}
*/
/*********************************************************/
/*********************************************************/
