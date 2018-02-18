/*  File: systags.h
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov  1 11:47 1995 (rd)
 * Created: Tue Aug 30 19:17:15 1994 (mieg)
 *-------------------------------------------------------------------
 */

/*
   The systags are defined here are initialised in
   whooks/sysclass.c

   ALLWAYS edit the 2 files in a correlated way.
*/

/* $Id: systags.h,v 1.5 2008/01/20 01:42:17 mieg Exp $ */
#ifndef SYSTAGS_DEF_H
#define SYSTAGS_DEF_H

#define _Text  1
#define _AddressText  2
#define _Greek  3
#define _Russian  4
#define _Text1  5
#define _Text2  6
#define _Text3  7
#define __DNA1  8
#define __DNA2  9
#define __DNA3  10
#define __DNA4  11
#define __DNA5  12
#define __DNA6  13
#define __DNA7  14
#define __RNA1  15
#define __RNA2  16
#define __RNA3  17
#define __RNA4  18
#define __Protein1  19
#define __Protein2  20
#define __Protein3  21
#define _NextC  22
#define _LastC  23
		/*
		* all before _LastC are TEXT data fields.  If you say
		* TEXT in your model, it uses _Text.  The rest are
		* essentially unused.
		*/

#define _Int  24
#define _Unsigned  25
#define _Long  26		/* not supported */
#define _Long_Unsigned  27	/* not supported */
#define _Float  28
#define _DateType  29
#define _continuationKey  30
#define _LastN  31
		/*
		* all from _LastC to _LastN are numeric data items
		*/

#define _bsHere  32  /* Private to the kernel, used by queryexe */
#define _bsRight  33
#define _bsDown  34
#define ___sys35  35	/* available for reassignment */
#define ___sys36  36	/* available for reassignment */
#define ___sys37  37	/* available for reassignment */
#define ___sys38  38	/* available for reassignment */
#define ___sys39  39	/* available for reassignment */
#define _UNIQUE  40
#define _XREF  41
#define _ANY  42
#define _FREE  43
#define _REPEAT 44
#define _COORD 45  /* Int or Float subject to bsCoordShift() */
#define ___sys46  46	/* available for reassignment */
#define ___sys47  47	/* available for reassignment */
#define ___sys48  48	/* available for reassignment */
#define ___sys49  49	/* available for reassignment */
		/*
		* above here are "model modifiers"
		*/

		/*
		* below here are real tag names
		*/
#define _Date  50
#define _User  51
#define _Session 52
#define _BatPlus  53
#define _BatMinus  54
#define _GlobalBat 55
#define _Session_Title  56
#define _SessionLex 57
#define _VocLex  58
#define _GlobalLex  59
#define _Quoted_in  60
#define _CodeRelease 61
#define _DataRelease 62
#define _Created_from 63
#define _Image        64
#define _Pick_me_to_call  65	
#define _File 66
#define _Non_graphic 67
#define _Centre 68
#define _Destroyed_by_session 69
#define _Up_linked_to 70
#define _Permanent_session 71
#define _Secret  72
#define _This_session  73
#define _Related_tags 74
#define _Used_for 75
#define _Appears_in_source_code 76
#define _Uses_tags 77
#define _SourceFile 78
#define _Includes 79
#define _Included_by 80
#define _Parent_tag 81
#define _Appears_in_class 82
#define _Compiled_as 83
#define _Is_a_subclass_of 84
#define _Is_a_superclass_of 85
#define _Uses_class 86
#define _Visible 87
#define _Hidden 88
#define _Mask 89
#define _Belongs_to_class 90
#define _Filter 91
#define _Visibility 92
#define _Start 93
#define _Finish 94
#define _PROTECTED 95
#define _Constraints 96
#define _File_name 97
#define _IndexVersion 98

  /* These colors match in the same order those declared 
   * whooks/systags.h
   * whooks/sysclass.c twice
   * wh/colours.h
   * wh/colorRGB.h
   * w2/graphcolour.c  
   * w4/palette.c
   * waceview/swfc.c
   */	
#define _WHITE 100
#define _BLACK 101
#define _LIGHTGRAY 102
#define _DARKGRAY 103
#define _RED 104
#define _GREEN 105
#define _BLUE 106 
#define _YELLOW 107
#define _CYAN 108
#define _MAGENTA 109
#define _LIGHTRED 110
#define _LIGHTGREEN 111
#define _LIGHTBLUE 112
#define _DARKRED 113
#define _DARKGREEN 114
#define _DARKBLUE 115
#define _PALERED 116
#define _PALEGREEN 117
#define _PALEBLUE 118
#define _PALEYELLOW 119
#define _PALECYAN 120
#define _PALEMAGENTA 121
#define _BROWN 122
#define _ORANGE 123
#define _PALEORANGE 124
#define _PURPLE 125
#define _VIOLET 126
#define _PALEVIOLET 127
#define _GRAY 128
#define _PALEGRAY 129
#define _CERISE 130
#define _MIDBLUE 131
#define _LIGHTMAGENTA 132
#define _LIGHTCYAN 133
#define _DARKVIOLET 134
#define _LAVANDER 135
#define _BLUE1 136
#define _BLUE2 137
#define _BLUE3 138
#define _BLUE4 139
#define _BLUE5 140
#define _BLUE6 141
#define _BLUE7 142
#define _BLUE8 143
#define _GREEN1 144
#define _GREEN2 145
#define _GREEN3 146
#define _GREEN4 147
#define _GREEN5 148
#define _GREEN6 149
#define _GREEN7 150
#define _GREEN8 151
#define _RED1 152
#define _RED2 153
#define _RED3 154
#define _RED4 155
#define _RED5 156
#define _RED6 157
#define _RED7 158
#define _RED8 159
#define _LIGHTVIOLET 160
#define _DARKCYAN 161
#define _LIGHTORANGE 162

  /* tags 160 to 164 are reserved for 64 possible colors */

#define _LastSystemTag 499

	/*
	* special keys in class 1 - used for bootstrapping
	*/
/*GLOBAL_LEXIQUES 1<<24 ==  16777216  */

#define __Global_Bat  16777217
#define __lexi1   16777218
#define __lexa1   16777219
#define __voc1    16777220
#define __lexi2   16777221
#define __lexa2   16777222
#define __voc2    16777223
#define __lexi3   16777224
#define __lexa3   16777225
#define __voc3   16777226
#define __batPlus 16777227
#define __batMinus 16777228
#define __oldPlus 16777229
#define __oldMinus 16777230
#define __lexh1   16777231
#define __lexh2   16777232
#define __lexh3   16777233
/* the order of these number is not random but crucial
 * for back compatibility with prerelease use of the code at sanger
 */
#define __lext2   16777234
#define __lext3   16777235
#define __lext1   16777236
#define __superKey   16777237

/***************** end of file *******************/
#endif
