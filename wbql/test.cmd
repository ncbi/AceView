comment "Select from where first examples"
select p,pub  from p in class person, pub in p->papers where count { p->papers } >= 2  
select p,w,h from p in class person, w in p->weight, h in p->height where w > 150 * (h - 1) and w > 80  

comment "Short forms, acting on the active list"
Find person   // list all members of class person
select K*     // limit to names starting with K
follow papers  // list their papers
list
select ?person k* ; >papers    // idem on a single line
select ?Paper
list

comment "Dealing with constants"
select x, y, z from x = 2 , y = "hello", z = 2 * 3 + 4
select L from L in class line , pi = 3.14 where L->length > 2 * pi / 3
select d from  d = 2 where d == 5 - 3   // correct filter
// select d from  d = 2 where d = 5 - 3    // error, please use == 
select x from x = 2 // full form
select x = 2   // short form
select 2       // shortest form

comment "TITLE"
S -title y from x=1,y=2,z=3 TITLE x:uu,y:vv,z:ww
S -title y,z from x=1,y=2,z=3 TITLE x:uu,y:vv,z:ww
S -title z,y from x=1,y=2,z=3 TITLE x:uu,y:vv,z:ww
S -title x from x=1,y=2,z=3 TITLE x:uu,y:vv,z:ww


comment "arithmetics"
select  (1 + 2) * (3 - ( 4 + 5))  // returns  (3) * (-6) = -18
select  8 modulo 3      // returns 2
select  -1 modulo 3     // returns  2 (math convention) rather than  -1
select  2 * 5 modulo 3  // returns  1 
select  9 ^ 2           // returns square(9) = 9 * 9 = 81
select  9 ^ (1/2)       // returns sqrt(9) = 3
select  1, (-1) ^ .5    // NULL sqrt(-1) is NaN
select  3*2 , (5+1)/2   // returns 6,3

comment "DNA and Proteins"
select ?sequence s* ; DNA      // export a set of DNA sequences


comment "Navigating through the acedb schema"
select p in ?Paper, a in p->Author, dr in a#Professor, z in a->Affiliation
select P in ?Paper, p1 in P->pages, p2 in p1[1] where p2
select P, p1, p2 from P in ?Paper, p in P#pages, p1 in p[1], p2 in p[2] where p2
select P, p1, p2 from P in ?Paper, p1 in P->pages[1], p2 in P->Pages[2] where p2
select P  in ?Paper, p in P#pages, p2 in p[2]
select c, g, a1, a2 from c in ?chromosome, g in c->Gene, a1 in g[1], a2 in g[2]
select c, g, a1, a2 from c in ?chromosome, g in c->Gene, a1 in g[1], a2 in a1[1]
select c in ?chromosome, g in c->Gene, a1 in g[1], a2 in a1[1]
select c, g, a1, a2 from c in ?chromosome, g in c->Gene, a1 in c->Gene[2], a2 in c->Gene[3]

comment "Transitive closure and multivalued data cell"
select p in ?Person, c in p->child where c
select ?Person Eva ; @ => child // Eva's children
select c from w in class "Person" where w == "Eva", c in w => child 
select c from w in class "Person" where w == "Eva", c in w>> child 
select ?paper p1 ; @=>Author

comment "Simplification rules"
select x from x in class "sequence"  // full syntax, or equivalently
select x in class "sequence"  // drop from since x is exported, or
select x in class  sequence   // drop the the quotes. or
select class sequence         // drop the name of the variable, or
select ?sequence              // use ? to represent class
find sequence                 // intuitive equivalent syntax
list

comment "Name recognition"
select x from x in class "person"  where x like "king"  // full syntax
select ?person like "king"   // drop x  and the where keyword
select ?person  king         // drop like and the quotes
find Person king             // intuitive equivalent syntax
list
select ?Person "K?ng"        // select King, Kong but not Kuong 
select ?Person "k*ng"        // select King, Kong, Kuong and "Karl Maning" 
select ?Person < kj          // select Karl, King but not kong
select ?Person =~ '^k.*g'    // =~ invokes Unix regular expressions

comment "protecting named variables and reserved words"
select x in ?Person, y in x->Brother where y > x    // x is a variable
select x in ?Person, y in x->Brother where y > "x"  // "x" is letter x

comment "Chaining queries"
select a from a in class "person" ;   // populate the active list @, stay silent
select a from a in @ where a like "king"   // derive a from the list
find person
select king

comment "Active list is populated by the first selected variable"
select p, j from p in class paper, j in p->journal  // case PJ
select j, p from p in class paper, j in p->journal  // case JP
select j    from p in class paper, j in p->journal  // case J
select a in class "Person" ; select b from b in @ where b like "k*" 
select ?person ; "k*"    //  short equivalent form  
select ?person ; k*    //  short equivalent form    

comment "Filtering with a where clause"
select ?Person Jim 
select p in class "Person" where p == "Jim"
 
comment "Filter on Text"
select p in ?Person where p
select p in ?Person where p#affiliation
select p in ?Person where p->affiliation
select p,q from p in ?Person, q in p#affiliation where q
select p,q from p in ?Person, q in p->affiliation where q
select ?person ; w in @->affiliation where w

comment "Filter on ?Text"
select p in ?Person where p#nickname
select p in ?Person where p->nickname
select p,q from p in ?Person, q in p#nickname where q
select p,q from p in ?Person, q in p->nickname where q
select ?person ; w in @->nickname where w

comment "Filter on Int"
select p in ?Person where p#size
select p in ?Person where p->size
select p,q from p in ?Person, q in p#size where q
select p,q from p in ?Person, q in p->size where q
select p,q from p in ?Person, q in p->size where q > 0 
select p,q from p in ?Person, q in p->size where q == 0
select ?person ; w in @->size where w

comment "Filter on Float"
select p in ?Person where p#weight
select p in ?Person where p->weight
select p,q from p in ?Person, q in p#weight where q
select p,q from p in ?Person, q in p->weight where q
select p,q from p in ?Person, q in p->weight where q > 0
select p,q from p in ?Person, q in p->weight where q == 0
select ?person ; w in @->weight where w

comment "Filter on DateType"
select p in ?Person where p#Birth
select p in ?Person where p->Birth
select p,q from p in ?Person, q in p#Birth where q
select p,q from p in ?Person, q in p->Birth where q
select p,q from p in ?Person, q in p->Birth where q > `2000`
select p,q from p in ?Person, q in p->Birth where q > `2001-03`
select p,q from p in ?Person, q in p->Birth where q > `2001-03-05`
select p,q from p in ?Person, q in p->Birth where q >~ `2001-03-05`
select ?person ; w in @->Birth where w

comment "Dumpy persons"
select p in ?Person, h in p->height, w in p->weight where 2 * (h - 100) < 11.3 * w 

comment "Jokers in text patterns"
select p in ?Person where    p ~ "T?m"     // selects Tam, Tim and Tom,  but not Attim which does not start with T
select p in ?Person where     p ~ "T*m"     // also selects Theotym, but not Thomas which does not end with m
select p in ?Person where     p ~ "*T?m*"   // also selects Attim"
select p in ?Person where     p ~ "*T*m*"   // selects Tam,  Tim, Tom, Theotym, Thomas and  Attim

comment "UNIX regular expressions (equal tilde)"
select p in ?Person where       p =~ 't[io]m'    // selects Tim, Tom and Attim but not Tam
select p in ?Person where       p =~ '^t[io]m'   // selects Tim and Tom, but not Attim which does not start with T
select p in ?Person where       p =~ 't.*m'      // selects Tam,  Tim, Tom, Theotym, Thomas and  Attim
select p in ?Person where       p =~ 't.*m$'     // rejects Thomas which does not ends with m

comment "Spotting missing or out of range data"
select p,w from p in ?person, w in p->weight where w < 10 OR ! w

S ?person tom
S p from p in @ where COUNT {select q from q in p->papers } > 0   
S p from p in @ where COUNT {select q from q in p->papers where q->journal } > 0  

S p from p in ?person where COUNT {select q from q in p->papers where q->journal == nature} > 0
S p from p in ?person where COUNT {select q from q in p->papers } > 0       
S p from p in ?person where COUNT {select q from q in p->papers_junk } > 0  
S p,c from p in ?person, c in COUNT {select q from q in p->papers where q->journal == nature}
S p,c from p in ?person, c in COUNT {select q from q in p->papers }
S p,c from p in ?person, c in COUNT {select q from q in p->papers_junk} 
S p,c from p in ?person, c in COUNT p->papers where c          
S p,c from p in ?person, c in COUNT p->papers where c  > 0
S p,c from p in ?person, c in COUNT p->papers_junk
S p,c from p in ?person where p == tom, c in COUNT {select q from q in p->papers}

comment "acedb dynamic subclasses"
select p in class Prolific_author
select p from p in @ where p ISA Prolific_author

comment logical operators
select x from x=3, a=1, b=4 where (a <= x && x <= b) || (b < x && x < a) 

comment "Counting and embedded virtual queries"
select a in ?person where count a->papers > 1
select a in class person where count {select p in a->papers} > 1
select a in class "person" where count {select p in a->papers where p->journal == "nature"} > 1

comment "Dates: The successive numbers represent Y:year, M:month (1 to 12), D:day (1 to 31), h:hour (0 to 24), m:minute (0 to 60) and s:seconds (0 to 60)"
select d = `2016-01-30_22:47:15`
select d = `2016-1-7_3:7`  // No leading zeroes
select x from x="ok" where   `2016-2`  == `2016`       // ---> false
select x from x="ok" where   `2016-2`  >= `2016`       // ---> true
select x from x="ok" where   `2016-2`  <= `2016`       // ---> false

select x from x="ok" where   `2016-2`  =~ `2016`       // ---> true
select x from x="ok" where   `2016-2`  >~ `2016`       // ---> true
select x from x="ok" where   `2016-2`  <~ `2016`       // ---> true

select x from x="ok" where   `2016-2`  == `2016-2-17`  // ---> false
select x from x="ok" where   `2016-2`  == `2016-3`     // ---> false
select x from x="ok" where   `2016-2`  == `2016`       // ---> false

select x from x="ok" where   `2016-2`  =~ `2016-2-17`  // ---> true
select x from x="ok" where   `2016-2`  =~ `2016-3`     // ---> false
select x from x="ok" where   `2016-2`  =~ `2016`       // ---> true

comment "DNA and Proteins"
select s in class sequence, d in DNA(s)
select x = 3, y = 8, s in ?sequence, d in DNA(s,x,y)   
select s,dna, adn, c from s in ?sequence,dna in DNA(s,1,6), adn in DNA(s,6,1), c in s#cds
select PEPTIDE (@,1,2)  // all proteins are shown
select ?sequence CDS ; PEPTIDE // only coding (CDS) sequences are translated

comment "Ordering"
select x, y, z from x in ?Person, y in x->child, z in y->brother where z order_by 1+2+3  // default order
select x, y, z from x in ?Person, y in x->child, z in y->brother where z order_by -y+x+z  // modified order


comment Shortcuts"
select a from a in class "person" where a like "king*"
    select ?person ~ "king*"
    select ?person ~ king*
 select e in ?Person, b in e->Brother where b > e   // var b > var e
 select e in ?Person, b in e->Brother where b > "e" // var b > string e
    select ?person like King  // select person King 
    select p in ?person where p#Professor
    select ?person Professor  // select all persons with the tag
    select ?person ; Professor  // select all persons with the tag
    select ?person King   // which is interpreted as 
    select a in ?person where ( a#King  or a like "King" )
    select ?person t??   // selects tim tam tom
    select ?person < tb  // selects up to Tam, not Tim
    select a from a in class "person"
    select a from a in @ where a ~ "king"
    select ?person
    select king
    select ?paper, ?->journal
    select p, j from p in class paper, j in p->journal
    select j from p in class paper, j in p->journal
    select j, p from p in class paper, j in p->journal
    select j, p from p in class paper, j in p->journal where j
    select ?person ; ~ k* 
    select ?person ; "K*" 
    select a in class person ; select b in @ where b ~ "k*" 
    select a in class person ; b in @ where b ~ "k*" 
    select a in class person ; select b in @ where b == king
    select a in class person ; b in @ where b == king
    select a in class person ; select b in @ where b == "king"
    select a in class person ; b in @ where b == "king" 
    select ?person ; k* ; papers
    select a in class "person"      // select all persons
    select a in @ where a ~ "k*"    // members of previous list called k*
    select a in @ where a#papers    // limit to persons with a paper 
    select a in class "person" ;    // select all persons, output suppressed
    select a in class "person"  order_by -a    // select all persons, report them

