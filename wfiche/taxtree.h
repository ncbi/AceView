#ifndef TAXTREE_DEF
#define  TAXTREE_DEF
#if 1
/* do not edit randomly
 * the code depends on the initial AND terminal dots etc 
 * the complete tree is below
 *
 * If you add a line this affects
 *    kantormegaparse
 *      which will produce an ace file for that modified line
 * hence you MUST edit the acedb schema in 4 places (sorry)
 *   *Tax_tree AND Tax_count ) * (Class Product AND Class kantor)
 *
 *  If you suppress a line, this will create parse errors
 *  everywhere in the human pipeline, because of the huge
 *  hysteresis of the system, so beware.
 *
 * fortunately there is no need to touch the aceview code
 * or the titlelizer code (makemrna.c wfiche...c etc)
 */

static char *myTaxTree=
"Archaea ....................................  0 \n"
"Viruses ....................................  0 \n"
"Bacteria ...................................  0 \n"
"..Escherichia coli ......................... -1 \n"
"Eukaryota ..................................  0 \n"
"..Mycetozoa ................................  0 \n"
"....Dictyostelium discoideum ............... -1 \n"
"..Viridiplantae ............................  0 \n"
"....Chlamydomonas reinhardtii .............. -1 \n" 
"....Arabidopsis thaliana ................... -1 \n" 
"....Zea mays ............................... -1 \n" 
"..Fungi Metazoa group ......................  0 \n"
"....Fungi ..................................  0 \n"
"......Saccharomyces cerevisiae ............. -1 \n"
"......Schizosaccharomyces pombe ............ -1 \n"
"......Neurospora crassa .................... -1 \n"
"....Metazoa ................................  0 \n"
"......Bilateria ............................  0 \n"
"........Pseudocoelomata ....................  0 \n"
"..........Caenorhabditis elegans ........... -1 \n"
"........Protostomia ........................  0 \n"
"..........Mollusca .........................  0 \n"
"..........Panarthropoda ....................  0 \n"
"............Anopheles gambiae str. PEST .... -1 \n"
"............Drosophila melanogaster ........ -1 \n"
"........Deuterostomia ......................  0 \n" /* division in embryogenesis*/
"..........Echinodermata ....................  0 \n"
"..........Chordata .........................  0 \n" /* nerve */
"............Vertebrata .....................  0 \n" /* vertebra */
"..............Teleostomi ...................  0 \n" /* bones */
"................Teleostei...................  0 \n" /* Sushi */
"..................Takifugu rubripes ........ -1 \n" /* */
"..................Danio rerio .............. -1 \n" /* */
"................Amphibia ...................  0 \n" /* amphibians - lay eggs and get rid of them */
"..................Xenopus laevis ........... -1 \n"
"................Amniota ....................  0 \n" /* baby is developped by mothers resources */
"..................Mammalia .................  0 \n" /* wife of papalia */
"....................Eutheria ...............  0 \n" /* standard mammals */
"......................Mus musculus ......... -1 \n"
"......................Rattus norvegicus .... -1 \n"
"......................Homo sapiens ......... -1 \n"
 ; 
#else
/* do not delete, this is the complete taxonomy tree */
static char *myTaxTree=
"Archaea .....................................................  0 \n"
"Viruses .....................................................  0 \n"
"Bacteria ....................................................  0 \n"
"..Escherichia coli .......................................... -1 \n"
"Eukaryota ...................................................  0 \n"
"..Plasmodium falciparum ..................................... -1 \n" 
"..Mycetozoa .................................................  0 \n"
"....Dictyostelium discoideum ................................ -1 \n"
"..Viridiplantae .............................................  0 \n"
"....Brassicales .............................................  0 \n"  /* up limit of ara */
"......Arabidopsis thaliana .................................. -1 \n" 
"..Fungi Metazoa group .......................................  0 \n"
"....Fungi ...................................................  0 \n"
"......Saccharomyces cerevisiae .............................. -1 \n"
"......Schizosaccharomyces pombe ............................. -1 \n"
"....Metazoa .................................................  0 \n"
"......Eumetazoa .............................................  0 \n"
"........Cnidaria ............................................  0 \n"
"........Ctenophora ..........................................  0 \n"
"........Bilateria ...........................................  0 \n"
"..........Coelomata .........................................  0 \n" /* cavity in the body */
"............Deuterostomia ...................................  0 \n" /* division in embrygenesys*/
"..............Chordata ......................................  0 \n" /* nerve */
"................Craniata ....................................  0 \n" /* head */
"..................Vertebrata ................................  0 \n" /* vertebra */
"....................Gnathostomata ...........................  0 \n" /* jaws */
"......................Chondrichthyes ........................  0 \n" /* cartilages */
"......................Teleostomi ............................  0 \n" /* bones */
"........................Euteleostomi ........................  0 \n"
"..........................Takifugu rubripes ................. -1 \n" /* */
"..........................Sarcopterygii .....................  0 \n" /* arms and lges */
"............................Tetrapoda .......................  0 \n" /* four of them */
"..............................Amphibia ......................  0 \n" /* amphibians - lay eggs and get rid of them */
"..............................Amniota .......................  0 \n" /* baby is developped by mothers resources */
"................................Sauropsida ..................  0 \n" /* dynasaurus, lizard  branch */
"................................Mammalia ....................  0 \n" /* wife of papalia */
"..................................Prototheria ...............  0 \n" /* egg laying */
"..................................Theria ....................  0 \n" /* no eggs */
"....................................Metatheria ..............  0 \n" /*  cangoroo type babycare */
"....................................Eutheria ................  0 \n" /* standard mammals */
"......................................Carnivora .............  0 \n" /* dogs .. cats */ 
"......................................Cetartiodactyla .......  0 \n" /* whales, pigs */
"......................................Perissodactyla ........  0 \n" /* horses ... */
"......................................Proboscidea ...........  0 \n" /* elephand , longnoses */
"......................................Rodentia ..............  0 \n" /* rodents */ /* up limit of mouse/rat */
"........................................Mus musculus ........ -1 \n"
"........................................Rattus norvegicus.... -1 \n"
"......................................Insectivora ...........  0 \n" /* hedgehogs, insect eaters */
"......................................Lagomorpha ............  0 \n" /* rabit, hare */
"......................................Primates ..............  0 \n" /* primates */ /* up limit of human/monkeys */
"........................................Catarrhini ..........  0 \n" /* */
"..........................................Cercopithecidae ...  0 \n" /* old world monkeys */
"............................................Macaca mulatta .. -1 \n" /* my cousin */
"..........................................Hominidae .........  0 \n" /* apes */
"............................................Gorilla ......... -1 \n" /* my neighbour */
"............................................Pan ............. -1 \n" /* chimps */
"............................................Pongo ........... -1 \n" /* oran gutans */
"............................................Homo ............  0 \n" /* you */
"..............................................Homo sapiens .. -1 \n" /* me */
"..............Echinodermata .................................  0 \n"
"............Protostomia .....................................  0 \n"
"..............Mollusca ......................................  0 \n"
"..............Panarthropoda .................................  0 \n" /* up limit of droso */
"................Drosophila melanogaster ..................... -1 \n"
"..........Pseudocoelomata ...................................  0 \n" /* up limit of worms */
"............Nematoda ........................................  0 \n"
"..............Caenorhabditis elegans ........................ -1 \n" ; 

#endif 

#endif
