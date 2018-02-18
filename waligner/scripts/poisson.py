#!/usr/bin/env python
###############################################################################################
# poisson.py -- try to graphically fit the distribution of de-uno introns

# meta-x python-shell  gives an interactive python shell inside emacs
###############################################################################################
import getopt
import os
import os.path
import sys
import csv

import matplotlib as mplt
mplt.use('QTAgg')

import string
import time
import numpy as np

import matplotlib.pyplot as plt


if 0 :
    a = np.array([1,2,3,4,5,10,20,50,100,200]);
    b = np.array([1,.619,.442,.356,.302,.179,0.103,0.0449,0.0213,0.00820])
    
    
    plt.plot(np.log(a),1-b)
    plt.plot(np.log(a),1-b,"bo")
    plt.plot(np.log(a),1.3*(1-1/(np.sqrt(a))))
    plt.show()
    exit (0)


#plt.title(r'is there noise $N\;(1 - {1 \over K } \Sum \;N_k\;\exp (-x \. \alpha\.2^{-k}))$')

def p1(N,x):
    return N * (1 - np.exp(-x/N))

def E10(x):
    return 240000.0 + x/150000.0

def E10b(x):
    return 242000  + x/400000.0   - 5000 * np.exp( -(x - 0)/1000000000) -20000 * np.exp( -(x - 0)/300000000)

def E10av1(x):
    return  200000.0 - 1000 + x/30000.0   - 12000 * np.exp( -(x - 0)/80000000) - 20000 * np.exp( -(x - 0)/170000000) - 52000 * np.exp( -(x - 0)/20000000) - 50000 * np.exp( -(x - 0)/5000000)  - 70000 * np.exp( -(x - 0)/15000000) 


def E10a(x):
    return  270000.0 +  0*x/60000.0   - 30000 * np.exp( -(x - 0)/1280000000)  - 30000 * np.exp( -(x - 0)/640000000) - 30000 * np.exp( -(x - 0)/320000000)  - 30000 * np.exp( -(x - 0)/160000000) - 30000 * np.exp( -(x - 0)/80000000) - 30000 * np.exp( -(x - 0)/40000000) - 30000 * np.exp( -(x - 0)/20000000)  - 30000 * np.exp( -(x - 0)/10000000)  - 30000 * np.exp( -(x - 0)/5000000) 



# tt=np.genfromtxt('intronHistoMainNB.txt',delimiter='\t',autostrip=True,skip_header=3,comments="#",names=True,
#                 usecols=range(188)
#                 )

f='intronHistoMainNB.txt'
f='intronHistoMainNB.non_gt_ag.txt'
f='intronHistoMainNB.refseq_encode_magic.txt'
f='intronHistoMainNB.refseq_encode_magic_sorted_new.txt'
#f='intronHistoMainNB.refseq_encode_magic_just_gt_ag_sorted.txt'
#f='intronHistoTitrationACDB.gc_ag.txt'
#f='intronHistoTitrationACDB.acceptor.txt'
#f='intronHistoTitrationACDB.donor.txt'
#f='intronHistoTitrationACDB.fuzzy1000.txt'
#f='intronHistoT.countOnce.txt'

tt=np.genfromtxt(f,delimiter='\t',autostrip=True,skip_header=3,comments="#",names=True,
                 usecols=range(175)
                 )
# ,"Encode-2012 exon junctions dicovered in at least 1 reads"))

xx =100000000.0 * np.arange(60)
xx2=10000000.0 * np.arange(66,306,10)

xx1      =200000000.0 * np.arange(300)
xx10     =20000000.0 * np.arange(300)
xx20     =10000000.0 * np.arange(300)
xx50     =4000000.0 * np.arange(300)
xx100    =2000000.0 * np.arange(300)
xx1000   =400000.0 * np.arange(300)
xx10000  =40000.0 * np.arange(300)

# print tt
print type(tt)
print tt.shape
print tt.dtype.names
print tt[tt.dtype.names[0]]
# print tt['Encode2012_exon_junctions_dicovered_in_at_least_10_reads']


###### comparison of Encode RefSeq AceView genes recovered
if 1 :
    plt.plot ((1000*tt['kb_cumul']),tt['MagicNB2013_exon_junctions_dicovered_in_at_least_50_reads'],label='Magic_NB_2013')
    plt.plot ((1000*tt['kb_cumul']),tt['AceView2010_exon_junctions_dicovered_in_at_least_50_reads'],label='AceView_2010')
    plt.plot ((1000*tt['kb_cumul']),tt['Encode37702013_exon_junctions_dicovered_in_at_least_50_reads'],label='Encode_2013')
    plt.plot ((1000*tt['kb_cumul']),tt['RefSeq2012_exon_junctions_dicovered_in_at_least_50_reads'],label='RefSeq_2013')

    zz = np.ones(tt.size)
    plt.plot((1000*tt['kb_cumul']),469185 * zz,'b--')   # NB magic 100 support gt-ag
    plt.plot((1000*tt['kb_cumul']),369784 * zz,'g--')   # 369784 main chrom ? or ??, 376670 AceView May 31, 2012 SEQC freeze non-fuzzy gt_ag or gc_ag or at_ac
    plt.plot((1000*tt['kb_cumul']),345433 * zz,'r--')   # 345433 encode Homo_sapiens.GRCh37.70 on main chroms (dnd2dna -gtf)
    plt.plot((1000*tt['kb_cumul']),194113 * zz,'c--')   # RefSeq 2013 main chrom

if 0 :
    plt.plot (np.log(1000*tt['kb_cumul']),tt['AceView2010_exon_junctions_dicovered_in_at_least_1_reads'],label='AceView 2010')
    plt.plot (np.log(1000*tt['kb_cumul']),tt['Encode2012_exon_junctions_dicovered_in_at_least_1_reads'],label='Encode 2012')
    plt.plot (np.log(1000*tt['kb_cumul']),tt['RefSeq2012_exon_junctions_dicovered_in_at_least_1_reads'],label='RefSeq 2012')

    zz = np.ones(tt.size)
    plt.plot(np.log(1000*tt['kb_cumul']),383311 * zz,'b--')
    plt.plot(np.log(1000*tt['kb_cumul']),361177 * zz,'g--')
    plt.plot(np.log(1000*tt['kb_cumul']),194113 * zz,'r--')


###### 
# f-scaling P(n.x/n)
if 0 :
    plt.title('Scaled number of new introns discovered at least n times: E(n,x) = P(n,x/n)')
    plt.plot (tt['kb_cumul'],tt['New_exon_junctions_dicovered_in_at_least_1_reads'],label='E1')
    plt.plot (tt['kb_cumul']/2,tt['New_exon_junctions_dicovered_in_at_least_2_reads'],label='E2')
    plt.plot (tt['kb_cumul']/3,tt['New_exon_junctions_dicovered_in_at_least_3_reads'],label='E2')
    plt.plot (tt['kb_cumul']/4,tt['New_exon_junctions_dicovered_in_at_least_4_reads'],label='E2')
    plt.plot (tt['kb_cumul']/5,tt['New_exon_junctions_dicovered_in_at_least_5_reads'],label='E2')
    plt.plot (tt['kb_cumul']/6,tt['New_exon_junctions_dicovered_in_at_least_6_reads'],label='E2')
    plt.plot (tt['kb_cumul']/7,tt['New_exon_junctions_dicovered_in_at_least_7_reads'],label='E2')
    plt.plot (tt['kb_cumul']/8,tt['New_exon_junctions_dicovered_in_at_least_8_reads'],label='E2')
    plt.plot (tt['kb_cumul']/9,tt['New_exon_junctions_dicovered_in_at_least_9_reads'],label='E2')
    plt.plot (tt['kb_cumul']/10,tt['New_exon_junctions_dicovered_in_at_least_10_reads'],label='E10')
    plt.plot (tt['kb_cumul']/100,tt['New_exon_junctions_dicovered_in_at_least_100_reads'],label='E100')
    plt.plot (tt['kb_cumul']/1000,tt['New_exon_junctions_dicovered_in_at_least_1000_reads'],label='E1000')
    plt.plot (tt['kb_cumul']/10000,tt['New_exon_junctions_dicovered_in_at_least_10000_reads'],label='E10000')

# g-scaling P(n.x/n-1)
if 1 :
    a=.30
    plt.title('Scaled number of new introns seen at least n times: E(n,x) = P(n,x/n-1 + ' +str(a) + ')')
    plt.plot (tt['kb_cumul']/a,tt['New_exon_junctions_dicovered_in_at_least_1_reads'],label='E1')
    plt.plot (tt['kb_cumul']/(a+1),tt['New_exon_junctions_dicovered_in_at_least_2_reads'],label='E2')
    plt.plot (tt['kb_cumul']/(a+2),tt['New_exon_junctions_dicovered_in_at_least_3_reads'],label='E2')
    plt.plot (tt['kb_cumul']/(a+3),tt['New_exon_junctions_dicovered_in_at_least_4_reads'],label='E2')
    plt.plot (tt['kb_cumul']/(a+4),tt['New_exon_junctions_dicovered_in_at_least_5_reads'],label='E2')
    plt.plot (tt['kb_cumul']/(a+5),tt['New_exon_junctions_dicovered_in_at_least_6_reads'],label='E2')
    plt.plot (tt['kb_cumul']/(a+6),tt['New_exon_junctions_dicovered_in_at_least_7_reads'],label='E2')
    plt.plot (tt['kb_cumul']/(a+7),tt['New_exon_junctions_dicovered_in_at_least_8_reads'],label='E2')
    plt.plot (tt['kb_cumul']/(a+8),tt['New_exon_junctions_dicovered_in_at_least_9_reads'],label='E2')
    plt.plot (tt['kb_cumul']/(a+9),tt['New_exon_junctions_dicovered_in_at_least_10_reads'],label='E10')
    plt.plot (tt['kb_cumul']/(a+99),tt['New_exon_junctions_dicovered_in_at_least_100_reads'],label='E100')
    plt.plot (tt['kb_cumul']/(a+1000),tt['New_exon_junctions_dicovered_in_at_least_1000_reads'],label='E1000')
    plt.plot (tt['kb_cumul']/(a+10000),tt['New_exon_junctions_dicovered_in_at_least_10000_reads'],label='E10000')


# compute Q(L) by taking the derivative of P(1,x/a)
if 0 :
    a=.45
    a = 1
    xa = tt['kb_cumul']/a
    ya = tt['New_exon_junctions_dicovered_in_at_least_2_reads']
    na = tt.shape[0]
    print 'na=', na
    dya = np.zeros(shape=(na))
    dxa = np.zeros(shape=(na))
    za =  np.zeros(shape=(na))
    for i in np.arange(na-1) :
        dxa[i+1] = xa[i+1] - xa[i]
        dya[i+1] = ya[i+1] - ya[i]
        za[i] = dya[i+1]/dxa[i+1]
        # print i,xa[i],ya[i],dxa[i+1],dya[i+1],za[i]
        # print i,xa[i]
    za[na-1] = 0
    plt.plot (np.log(xa),np.log(za),label='Derivative')
    
    
if 0 :
    plt.plot(xx10+660000000, E10b(xx10),label='b')


    plt.plot(xx1/10,E10a(xx1),'r--',label='1')
    plt.plot(xx10,E10a(xx10),'r--',label='10')
    plt.plot(10*xx100,E10a(xx100),'r--',label='100')
    plt.plot(50*xx1000,E10a(xx1000) ,'r--',label='1000')
    plt.plot(400*xx10000,E10a(xx10000) ,'r--',label='10000')

type = 'Any'
ccumul = tt['kb_cumul']
cnew1 = tt['New_exon_junctions_dicovered_in_at_least_' + str(1) + '_reads']

for i in (2,3,4,5,6,7,8,9,10,20,50,100,200,500,1000) :
    cnew = tt['New_exon_junctions_dicovered_in_at_least_' + str(i) + '_reads']
    # plt.plot(1*ccumul, (i/1.0)*(cnew)/(cnew1+1),label="C" + str(i))    


for i in (1,) :
    cany = tt['Any_exon_junctions_dicovered_in_at_least_' + str(i) + '_reads']
    cnew2 = tt['New_exon_junctions_dicovered_in_at_least_' + str(2) + '_reads']
    cnew3 = tt['New_exon_junctions_dicovered_in_at_least_' + str(3) + '_reads']
    cnew4 = tt['New_exon_junctions_dicovered_in_at_least_' + str(4) + '_reads']
    cnew5 = tt['New_exon_junctions_dicovered_in_at_least_' + str(5) + '_reads']
    cnew6 = tt['New_exon_junctions_dicovered_in_at_least_' + str(6) + '_reads']
    cnew7 = tt['New_exon_junctions_dicovered_in_at_least_' + str(7) + '_reads']
    cnew8 = tt['New_exon_junctions_dicovered_in_at_least_' + str(8) + '_reads']
    cnew9 = tt['New_exon_junctions_dicovered_in_at_least_' + str(9) + '_reads']
    cnew10 = tt['New_exon_junctions_dicovered_in_at_least_' + str(10) + '_reads']

    #    cnew6 = tt['New_exon_junctions_dicovered_in_at_least_' + str(6) + '_reads']
    cnew = tt['New_exon_junctions_dicovered_in_at_least_' + str(i) + '_reads']
    # plt.plot (1*ccumul,1*cany - 0*cnew,label='A'+str(i))
    # plt.plot (1*ccumul,1*cany - 1*cnew,label='K'+str(i))
    # plt.plot (1*ccumul,0*cany + 1*cnew,label='N'+str(i))
    # plt.plot (1*ccumul,cnew/(1+cnew1),label='N'+str(i))


if 0 :
    plt.plot(1*ccumul, (3/1.0)*(cnew2 - cnew3)/(cnew1 - cnew2+1),label="M2")
    plt.plot(1*ccumul, (4/2.0)*(cnew3 - cnew4)/(cnew2 - cnew3+1),label="M3")
    plt.plot(1*ccumul, (5/3.0)*(cnew4 - cnew5)/(cnew3 - cnew4+1),label="M4")
    plt.plot(1*ccumul, (6/4.0)*(cnew5 - cnew6)/(cnew4 - cnew5+1),label="M5")
    plt.plot(1*ccumul, (7/5.0)*(cnew6 - cnew7)/(cnew5 - cnew6+1),label="M6")
    plt.plot(1*ccumul, (8/6.0)*(cnew7 - cnew8)/(cnew6 - cnew7+1),label="M7")
    plt.plot(1*ccumul, (9/7.0)*(cnew8 - cnew9)/(cnew7 - cnew8+1),label="M8")
    plt.plot(1*ccumul, (10/8.0)*(cnew9 - cnew10)/(cnew8 - cnew9+1),label="M9")

if 0 :
    plt.plot(1*ccumul, (3/1.0)*(cnew2 - cnew3)/(cnew1 - cnew2+1),label="M2")
    plt.plot(1*ccumul, (4/2.0)*(cnew3 - cnew4)/(cnew2 - cnew3+1),label="M3")
    plt.plot(1*ccumul, (5/3.0)*(cnew4 - cnew5)/(cnew3 - cnew4+1),label="M4")
    plt.plot(1*ccumul, (6/4.0)*(cnew5 - cnew6)/(cnew4 - cnew5+1),label="M5")
    plt.plot(1*ccumul, (7/5.0)*(cnew6 - cnew7)/(cnew5 - cnew6+1),label="M6")
    plt.plot(1*ccumul, (8/6.0)*(cnew7 - cnew8)/(cnew6 - cnew7+1),label="M7")
    plt.plot(1*ccumul, (9/7.0)*(cnew8 - cnew9)/(cnew7 - cnew8+1),label="M8")
    plt.plot(1*ccumul, (10/8.0)*(cnew9 - cnew10)/(cnew8 - cnew9+1),label="M9")



def multiP(N0,m,n,x,K):
    z = 0
    z1 = np.exp(-m)
    zn = z1
    for k in range(K):
        u = 1.0
        m1 = 0.0 + np.copy(m)
        for i in range(2,n+1,1) :
            u += m1
            if 0 :
                print 'i=',i,'m1=',m1,'u=',u
            m1 *= m/i
        z += u * zn
        if 0 :
            print 'k=',k,'z=',z,'zn=',zn
        zn *= np.sqrt(zn)
        m = m*np.sqrt(2) 
    return N0 * ( 1 - z/K)  

def ENew2(x,n,a) :
    N0 = 10300000
    K = 35
    alpha = a*32*512*13000000000.0
    m = x / alpha
    if n == 2 :
        N0 *= 0.6224
    if n == 3 :
        N0 *= 0.4401
    if n == 4 :
        N0 *= 0.351
    if n == 5 :
        N0 *= 0.296
    # N0 = 9000000
    return multiP(N0,m,n,x,K)       

def ENew1(x,n) :
    return  .6 * ENew2(x,n,2) + .15 * ENew2(x,n,.10)
    return  .2 * ENew2(x,n,1.4) + .17 * ENew2(x,n,.10)

def EKnown1(x,n) :
    N0 = 340000
    K = 33
    K = 1
    alpha = 1700000000.0
    alpha = alpha/13 # gros bug la courbe N = 5 est negative (-300k) pour x = 10^9
    m = x / alpha
    if n == 27 :
        N0 /= 1.015
    if n == 57 :
        N0 /= 1.07
    if n == 107 :
        N0 /= 4.5
    return multiP(N0,m,n,x,K)       

def EAny1(x,n):
    return ENew1(x,n) + EKnown1(x,n)


   

def testP(x):
    m=2*x
    if 1 :
        m1 = 0.0 + np.copy(m)
        for i in range(5) :
            if 1 :
                print 'i=',i,'m=',m,'m1=',m1
            m1 *= m/(i+1)
    return 1

if 0 :
    x=(1.0+np.arange(2))
    y=testP(x)
    print "x=",x,'y=',y
    x=1.0
    y=testP(x)
    print "x=",x,'y=',y
    x = 2
    y=testP(x)
    print "x=",x,'y=',y

    exit(0)
    y=EKnown1(x,5)
    print "x=",x,'y=',y
    x=2000000000
    y=EKnown1(x,5)
    print "x=",x,'y=',y
    
if 0 :
    plt.plot(xx20, ENew1(xx20,1),'r--',label='bN')
    plt.plot(xx20, ENew1(xx20,2),'g--',label='bN2')
    plt.plot(xx20, ENew1(xx20,3),'g--',label='bN3')
    plt.plot(xx20, ENew1(xx20,4),'g--',label='bN4')
    plt.plot(xx20, ENew1(xx20,5),'b--',label='bN5')
    # plt.plot(xx20, ENew1(xx20,10),'g--',label='b10')

if 0 :
    #plt.plot(xx20, EKnown1(xx20,1),'r--',label='bK')
    #plt.plot(xx20, EKnown1(xx20,2),'g--',label='bK2')
    plt.plot(xx20, EKnown1(xx20,5),'b--',label='bK5')
    #plt.plot(xx20, EKnown1(xx20,10),'g--',label='b10')

if 0 :
    plt.plot(xx20, EAny1(xx20,1),'r--',label='bA')
    plt.plot(xx20, EAny1(xx20,2),'g--',label='bA2')
    plt.plot(xx20, EAny1(xx20,3),'b--',label='bA3')
    plt.plot(xx20, EAny1(xx20,4),'b--',label='bA4')
    plt.plot(xx20, EAny1(xx20,5),'b--',label='bA5')
    #plt.plot(xx20, EAny1(xx20,10),'g--',label='b10')

print range(8)
plt.legend(loc='lower center',shadow=True)
plt.ylabel('de novo introns')
plt.xlabel('Cumulated sequencing in Terabases')


# plt.title('is there noise')
# plt.title('Number of junctions supported at least once\n among 194,113 junctions annotated in RefSeq-2012\n394,249 junctions annotated in Encode-2012\n383,111 junctions annotated in Aceview-2010')

# plt.text(1000000000,550000,r'fit: $ N \; (1 \;- \; \sum_{k=1}^{31} \;  e^{- \.2^{k/2}\. \alpha x } )   $')


plt.show()
exit(0)

lines=plt.show()
plt.setp(lines,linewidth=1.0)

exit(0)

