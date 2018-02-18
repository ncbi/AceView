#ifndef _seqIOABI_h_
#define _seqIOABI_h_

#include <sys/types.h> /* off_t */
#include "os.h"

/*
 * The ABI magic number - "ABIF"
 */
#define ABI_MAGIC	((int_4) ((((('A'<<8)+'B')<<8)+'I')<<8)+'F')

/*
 * The index is located towards the end of the ABI trace file.
 * It's location is given by a longword at a fixed place.
 */
#define IndexPO ((off_t)26)

#define IndexEntryLength 28

/*
 * Here are some labels we will be looking for, four chars packed
 * into an int_4
 */
#define DataEntryLabel    ((int_4) ((((('D'<<8)+'A')<<8)+'T')<<8)+'A')
#define BaseEntryLabel    ((int_4) ((((('P'<<8)+'B')<<8)+'A')<<8)+'S')
#define BasePosEntryLabel ((int_4) ((((('P'<<8)+'L')<<8)+'O')<<8)+'C')
#define SpacingEntryLabel ((int_4) ((((('S'<<8)+'P')<<8)+'A')<<8)+'C')
#define SignalEntryLabel  ((int_4) ((((('S'<<8)+'/')<<8)+'N')<<8)+'%')
#define FWO_Label         ((int_4) ((((('F'<<8)+'W')<<8)+'O')<<8)+'_')
#define MCHNLabel         ((int_4) ((((('M'<<8)+'C')<<8)+'H')<<8)+'N')
#define PDMFLabel         ((int_4) ((((('P'<<8)+'D')<<8)+'M')<<8)+'F')
#define SMPLLabel         ((int_4) ((((('S'<<8)+'M')<<8)+'P')<<8)+'L')
#define PPOSLabel         ((int_4) ((((('P'<<8)+'P')<<8)+'O')<<8)+'S')
#define CMNTLabel         ((int_4) ((((('C'<<8)+'M')<<8)+'N')<<8)+'T')

/*
 * From the ABI results file connected to `fp' whose index starts
 * at byte offset `indexO', return in `val' the `lw'th long word
 * from the `count'th entry labelled `label'.
 * The result is 0 for failure, or index offset for success.
 */
int getABIIndexEntryLW(FILE *fp, off_t indexO,
		       uint_4 label, uint_4 count, int lw,
		       uint_4 *val);

/*
 * Gets the offset of the ABI index.
 * Returns -1 for failure, 0 for success.
 */
int getABIIndexOffset(FILE *fp, uint_4 *indexO);

/*
 * Get an "ABI String". These strings are either pointed to by the index
 * offset, or held in the offset itself when the string is <= 4 characters.
 * The first byte of the string determines its length.
 * 'string' is a buffer 256 characters long.
 *
 * Returns -1 for failure, string length for success.
 */
int getABIString(FILE *fp, off_t indexO, uint_4 label, uint_4 count,
		 char *string);

int dump_labels(FILE *fp, off_t indexO);

#endif /* _seqIOABI_h_ */
