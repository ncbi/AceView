/*

generate_pseudoword
	Returns a string that looks like a pronouncable word.  n
	determines which word to return.

	People remember words more easily than numbers.  If you look at
	a number, any sequence of digits is just as good as any other.
	Words are chosen from a known set that follow certain rules.
	It will be easier to remember words that look "real", because
	people will treat them as new words that they just don't know
	yet.

	In fact, we are just making a number base that uses phonemes
	as digits, but nobody knows this.  (Actually, I know this
	and it still works.)

	This code supports different languages.  Presently, it can
	generate words that resemble English or Japanese.  You pass
	a parameter stating which you want.

English:

	Depending who you ask, English has somewhere in the range of 42
	to 45 phonemes.  I started with a fairly common set that has
	44 phonemes.

	I then reduced the set by removing sounds that I did not think
	I could unambiguously represent using the English alphabet.
	For example, the "u" sound in "put" has two standard spellings:
	"u" (put) and "oo" (wood).  But each of those spellings is also
	used for another sound ("ooze", "puke").

	I then further reduced the set by removing sounds that I anticipate
	would be difficult for speakers of English as a second language.
	For example, I am working for French scientists who inform me
	that the "th" sound is awful for French-speakers.  

	The coding is represented by using consonant-vowel pairs as
	digits in the numbering system.  The number base is number of
	vowels multiplied by the number of consonants.  This results in
	words that go consonant-vowel-consonant-vowel, which works quite
	nicely in English.

Japanese:

	Two of the three Japanese alphabets currently in use are
	strictly phonetic.  The list of letters in the alphabet is the
	list of valid phonemes.  The alphabets are named Katakana and
	Hirigana.  As far as I can tell, they are really two identical
	alphabets with different glyphs (pictures) for the same
	letters.  For example, in Katakana the sound "yo" looks like a
	mirror image of an English capital E.  In Hirigana, the same
	sound looks like the currency symbol for the British pound, a
	cursive L with a horizontal line.  (The third alphabet is named
	Kanji, and consists of a large number of a large number of
	Chinese symbols.)

	They say Japanese has "46 phonemes".  The are made of 5 vowels
	and 10 consonants.  A phoneme can be a vowel or a consonant
	followed by a vowel, but certain consonants cannot be followed
	by certain vowels.

	Because some consonants cannot have certain vowels after them,
	I am handling pseudo-japanese with a totally different
	algorithm.  There is just a table of sounds to choose from.

	I am leaving out the vowel-only sounds because that implies
	that some of the words we generate would be "aa", "aaa",
	"aaaa".  

	I am leaving out "wo", which is valid according to the hirigana
	table I am looking at, because I have a second hand report that
	it is not good.

	I am leaving out "Na", which is valid according to
	the hirigana table I am looking at, because it is too
	much like "na" for an English-speaker.
	
Notes on the author:

	My first language is US English, with the accent of the
	Delaware/Maryland area (north/east of Washington DC).  I
	know some German and a bit of Spanish.  

*/

#include <stdio.h>

#include "regular.h"

/*
* pseudo-English phoneme set
*/

static char * english_vowels[] =
        {
        "a",   /* at, sat man */
        "y",   /* it, wish, rib */
        "u",   /* up, but */
        "o",   /* so, low, hope */
        "ar",   /* art, farm, mar */
        "aw",   /* awful, bawl, law */
        "ee",   /* eel, mean, city */
        "er",   /* her, fur, skirt */
        "ey",   /* eye, sigh, time */
        "or",   /* or, tore, floor */
        "oy",   /* oil, boy, loyal */
        };


static char *english_consonants[] =
	{
	/*
	* basic phonemes
	*/
	"b",	/* bet, tab, lobby */
	"d",	/* dill, lid, older */
	"f",	/* fun, muff, iffy */
	"g",	/* get, tag, piggy */
	"j",	/* jean, jump, judge */
	"k",	/* can, lick, picker */
	"l",	/* let, tilt, pullet */
	"m",	/* me, rum, timer */
	"n",	/* no, tan, funny */
	"p",	/* pin, nap, happy */
	"r",	/* run, nor, pert */
	"s",	/* sap, pass, fist */
	"t",	/* ton, pit, metal */
	"v",	/* van, love, cover */
	"w",	/* win, swim, wash */
	"z",	/* zoo, ooze, fuzzy */
	"ch",	/* chin, church, itch */
	"sh",	/* she, push, worship */

	/*
	* double consonants - we can get away with using
	* these because they can occur anywhere in a word.
	* There are additional double consonant sounds
	* that can only occur at the beginning or end
	* of a word, but that makes the algorithm much
	* more complex.
	*/
	"bl",	/* blink */
	"fl",	/* fly */
	"gl",	/* glue, glide */
	"kl",	/* clip, klein*/
	"pl",	/* pleasure */
	"sk",	/* skill */
	"sl",	/* slop */
	"sm",	/* smart */
	"sn",	/* snark */
	"sp",	/* spell */
	"st",	/* strange */
	"sw",	/* swear */

	};


#define n_english_vowels	( sizeof(english_vowels) / sizeof(char *))
#define n_english_consonants	( sizeof(english_consonants) / sizeof(char *) )
#define english_number_base 	( n_english_vowels * n_english_consonants )

/*
* make the number base visible to outsiders so they can skip the
* single consonant-vowel words
*/
int pseudoword_english_number_base = english_number_base;

/*
* pseudo-japanese phoneme set, written with English letters
* so you can represent it in ASCII or ISO-8859-1 character
* encodings.
*/

static char *japanese_sounds[] =
	{
	"wa",
	"ra", "ri", "ru", "re", "ro",
	"ya",       "yu",       "yo",
	"ma", "mi", "mu", "me", "mo",
	"ha", "hi", "hu", "he", "ho",
	"na", "ni", "nu", "ne", "no",
	"ta", "ti", "tu", "te", "to",
	"sa", "si", "su", "se", "so",
	"ka", "ki", "ku", "ke", "ko",

	};

#define japanese_number_base	( sizeof(japanese_sounds) / sizeof(char *) )

/*
* make the number base visible to outsiders so they can skip the
* single phoneme words
*/
int pseudoword_japanese_number_base = japanese_number_base;

/*
* Generate a word.
*
* p is pointer to a character buffer to return the word in.  pass in
*	It can take 4 
*	characters per "digit" where there are NCONS * NVOWEL possible
*	digits.  2 000 000 000 is "koyzhuhshuhfar", so 20 characters 
*	should be sufficient, but why not make 30, 50 or 100 available?
*	If you're so low on memory that that would be a problem,
*	you should carefully analyze the algorithm to see how many
*	you will need.
* n is the number to be translated into a word.  n >= 0.
* phoneme_set selects what language you would like your word to
*	resemble.
*
*/

char *
generate_pseudoword(char *p, int n, enum pseudoword_phoneme_set phoneme_set )
{
	char *ret;
	int this_digit, v, c;

	if (!p) messcrash ("null return buffer p passed to generate_pseudoword") ;
	ret = p;

	switch ( phoneme_set )
		{
        case pseudoword_english:

		if (n == 0)
			{
			strcpy(ret,english_consonants[0]);
			strcat(ret,english_vowels[0]);
			return ret;
			}

		p[0] = 0;

		while (n > 0)
			{
			this_digit = n % english_number_base;
			n = n / english_number_base;

			v = this_digit % n_english_vowels;
			c = this_digit / n_english_vowels;

			strcpy(p,english_consonants[c]);
			while (*p) p++;

			strcpy(p,english_vowels[v]);
			while (*p) p++;
			}

		return ret;

	case pseudoword_japanese:

		if (n == 0)
			{
			strcpy(ret,japanese_sounds[0]);
			return ret;
			}

		p[0] = 0;

		while (n > 0)
			{
			this_digit = n % japanese_number_base;
			n = n / japanese_number_base;

			strcpy(p,japanese_sounds[this_digit]);
			while (*p) p++;
			}


		return ret;

	default:
		messcrash("invalid phoneme parameter %d in pseudoword.c\n",phoneme_set);
		return NULL;
		}
}

#if 0

In principle, you can take the pseudo-word and turn it back into a number
if you know which pseudo-language it came from.  There are words that
exist in both pseudo-word sequences, so you need to know which language
it was trying to look like.

I am turning off this code for now because it is a lot of work to keep it
working with all the changes I am making to support pseudo-japanese.

static int consnum(char *p, int max)
{
	int x;
	for (x=max-1; x >= 0; x--)
		if (strncmp(consonants[x],p,strlen(consonants[x])) == 0)
			return x;
	return -1;
}

static int vowelnum(char *p, int max)
{
	int x;
	for (x=max-1; x >= 0; x--)
		if (strncmp(vowels[x],p,strlen(vowels[x])) == 0)
			return x;
	return -1;
}


int decode_pseudoword(char *p, enum pseudoword_phoneme_set phoneme_set )
{
        int multiplier = 1;
        int v, c;
        int result = 0;

	select_phonemes(phoneme_set);

        while (*p)
                {
		c = consnum(p, nconsonants);
		if (c < 0) goto bad;
		while (p[strlen(consonants[c])] && vowelnum(p+strlen(consonants[c]), nvowels) < 0)
			{
			c = consnum(p, c);
			if (c < 0)
				return -1;
			}
		p += strlen(consonants[c]);

		v = vowelnum(p, nvowels);
		if (v < 0) goto bad;
		while (p[strlen(vowels[v])] && consnum(p+strlen(vowels[v]), nconsonants) < 0)
			{
			v = vowelnum(p, v);
			if (v < 0)
				return -1;
			}
		p += strlen(vowels[v]);

                result = result + ( v + ( c * nvowels ) ) * multiplier;
                multiplier = multiplier * nvowels * nconsonants;
                }

        return result;
bad:
	return -1;
}

#else

int decode_pseudoword(char *p, enum pseudoword_phoneme_set phoneme_set )
{
return -1;
}

#endif
