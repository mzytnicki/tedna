#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "globals.hpp"


unsigned short Globals::KMER                     = 61;
int            Globals::NB_THREADS               = 2;
long           Globals::SIZE_THREAD              = 10000000;
unsigned long  Globals::NB_READS                 = 0;
KmerNb         Globals::MIN_COUNT                = 3;
float          Globals::NB_REPETITIONS           = 2;
float          Globals::FREQUENCY_DIFFERENCE     = 2.5;
int            Globals::MIN_NB_NODES             = 500;
unsigned int   Globals::MAX_NB_NODES             = 100000;
unsigned int   Globals::MAX_PATHS                = 100;
unsigned int   Globals::MIN_TE_SIZE              = 500;
unsigned int   Globals::MAX_TE_SIZE              = 30000;
unsigned int   Globals::NB_SMALL_GRAPHS          = 10000;
short          Globals::EROSION_STRENGTH         = 100;
unsigned int   Globals::BUBBLE_SIZE              = 1000;
int            Globals::MIN_LTR_SIZE             = 50;
int            Globals::MAX_LTR_SIZE             = 5000;
unsigned short Globals::SHORT_KMER_SIZE          = 15;
Penalty        Globals::MAX_MERGE_SIZE           = 500;
Penalty        Globals::MIN_MERGE_SIZE           = 20;
Penalty        Globals::PENALTY_SIZE             = 1;
Penalty        Globals::PENALTY_MISMATCH         = 10;
Penalty        Globals::PENALTY_INDEL            = 30;
Penalty        Globals::MAX_PENALTY              = 500;
float          Globals::MIN_IDENTITY             = 0.8;
float          Globals::MAX_IDENTITY             = 0.3;
unsigned int   Globals::MERGE_MAX_NB             = 10000;
unsigned int   Globals::MERGE_MAX_NODES          = 10;
unsigned long  Globals::MIN_SCAFFOLD_KMERS       = 10;
unsigned long  Globals::MAX_SCAFFOLD_COUNTS      = 10000000;
unsigned int   Globals::SCAFFOLD_MAX_EV          = 5;
bool           Globals::FASTA_INPUT              = false;
string         Globals::CHECK;

const string Globals::VERSION = "1.2.4";
