#!/bin/bash
#Use this script with caution, 'right' outcome is not guarranteed!!!
#This script should be used under a Unix-like system(eg, Ubuntu) with both python and perl installed.
#You need run 'chmod +x ./PKSclusterfinder.sh' to get program executeble if this is the first time.
#intermediate--|-script
#              |-file
#Usage ./PKSclusterfinder.sh PfamA.hmm(1) KSN.hmm(2) KSC.hmm(3) input_protein.faa(4) input_feature_table.txt(5) output(6) BGCborderlimit(7)(default:20000)

#This is compiled by Chi Zhang @ Westlake University.

#check if key files exist
if [ ! -e $1 ]
then
    echo "PfamA file is missing."
    exit 1
fi

if [ ! -e $2 ]
then
    echo "First domain file is missing."
    exit 1
fi

if [ ! -e $3 ]
then
    echo "Second domain file is missing."
    exit 1
fi

if [ ! -e $4 ]
then
    echo "Protein file is missing."
    exit 1
fi

if [ ! -e $5 ]
then
    echo "Feature file is missing."
    exit 1
fi

if [ -e $6 ]
then
    echo "Warning output file exists!."
    echo 'Do you want to overwrite it? (Y/N):'
    read choice
    if [ "$choice" = 'Y' -o "$choice" = 'y' ]
    then
      echo "Overwrite and Continue!"
    else
      exit 1
    fi
fi

#check if prameter number is right and set BGC threshold to 20000 if not offered.
if [ "$#" -eq 6 ]
then
    set $1 $2 $3 $4 $5 $6 20000
elif [ "$#" -lt 6 -o "$#" -gt 7 ]
then
    echo "There should be 7 parameters! and the last one is BGC size(default:20000)"
    echo "You just entered $#"
    exit 1
fi

#01_hmmsearch to get proteins with offered domains
esl-sfetch --index $4
mkdir -p intermediate$$/files/01_hmmsearch/
hmmsearch --domE 1e-5 --domtblout intermediate$$/files/01_hmmsearch/"`basename $2`".dtbl $2 $4 >/dev/null
hmmsearch --domE 1e-5 --domtblout intermediate$$/files/01_hmmsearch/"`basename $3`".dtbl $3 $4 >/dev/null
grep -v "^#" intermediate$$/files/01_hmmsearch/"`basename $2`".dtbl | awk '{print $1}' >> intermediate$$/files/01_hmmsearch/firstdata.txt
grep -v "^#" intermediate$$/files/01_hmmsearch/"`basename $3`".dtbl | awk '{print $1}' >> intermediate$$/files/01_hmmsearch/firstdata.txt
sort -u intermediate$$/files/01_hmmsearch/firstdata.txt | esl-sfetch -f $4 - > intermediate$$/files/01_hmmsearch/earlysearch.fa

#02_hmmscan to confirm these proteins are indeed targets
mkdir -p intermediate$$/files/02_hmmscanfirst/
hmmscan --domtblout intermediate$$/files/02_hmmscanfirst/domtscanout.dtbl --noali -E 1e-5 $1 intermediate$$/files/01_hmmsearch/earlysearch.fa >/dev/null
cat intermediate$$/files/02_hmmscanfirst/domtscanout.dtbl | grep -v '^#' | awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' | sed 's/ /\t/g' | \
sort -k 3,3 -k 8n -k 9n | \
perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' | \
perl -e 'while(<>){chomp;@a=split(/\t/,$_);if(($a[-1]-$a[-2])>80){print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<=1e-5;}else{print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<=1e-3;}}' | awk '$NF>0.3' | sort -k 3,3 -k 8,8n -k 9,9n | awk -v OFS='\t' '{print $3,$1}' > intermediate$$/files/02_hmmscanfirst/parsedscanout

#03_check if hmmseach result matches hmmscan and filter.
mkdir intermediate$$/scripts/
echo "
inputfile1 = open('$2')
inputfile2 = open('$3')
inputfile3 = open('intermediate$$/files/02_hmmscanfirst/parsedscanout')
outputfile = open('intermediate$$/files/03_TargetID.txt', 'w')

lines1 = inputfile1.readlines()
lines2 = inputfile2.readlines()
lines3 = inputfile3.readlines()

inputfile1.close()
inputfile2.close()
inputfile3.close()

NameA=lines1[1].strip().split()[1]
NameB=lines2[1].strip().split()[1]

OUTBOX=[]
for record in lines3:
  if record.strip().split()[1] == ( NameA or NameB ) and record.strip().split()[0]+'\n' not in OUTBOX:
    OUTBOX.append(record.strip().split()[0]+'\n')

for STR in OUTBOX:
  outputfile.write(STR)

outputfile.close()
" > intermediate$$/scripts/03_hmmsearchchecker.py
python3 intermediate$$/scripts/03_hmmsearchchecker.py
#get target protein IDs

#04_get all protein information from genome and get target protein information.
grep -v "^#" $5 | grep "^CDS" | awk -F'\t' '{print $3, $7, $8, $9, $10, $11, $19}' | sort -k 1,1 -k 2,2 -k 3,3n > intermediate$$/files/04_AllInform.txt #sort was added 05202020.
#get the chromosome information
grep -Fwf intermediate$$/files/03_TargetID.txt intermediate$$/files/04_AllInform.txt | sort -k 1,1 -k 2,2 -k 3,3n > intermediate$$/files/04_TargetInform.txt
#get the chromsome information of target proteins

#05_decide whether two target proteins are close enough to make a cluster.
echo "
inputfile = open('intermediate$$/files/04_TargetInform.txt')
outputfile = open('intermediate$$/files/05_newdata.txt', 'w')
lines = inputfile.readlines()
lines = lines[:]+lines[-1:]
inputfile.close()
LENGTH=len(lines)
Result=[]
Box=[]
SPACE='\t'
Num = 0
while Num < (LENGTH-1):
  ASSEMB, CHROM, START, STOP, STRAND, ID, PROLENTH = lines[Num].strip().split()
  Num += 1
  ASSEMB1, CHROM1, START1, STOP1, STRAND1, ID1, PROLENTH1 = lines[Num].strip().split()
  if ASSEMB == ASSEMB1 and CHROM == CHROM1 and ( int(START1) - int(STOP) <= $7 ) and Num != (LENGTH-1):
    ADD = ASSEMB+SPACE+CHROM+SPACE+START+SPACE+STOP+SPACE+STRAND+SPACE+ID+SPACE+PROLENTH+SPACE+'B\n'
    Box.append(ADD)
  else:
    ADD = ASSEMB+SPACE+CHROM+SPACE+START+SPACE+STOP+SPACE+STRAND+SPACE+ID+SPACE+PROLENTH+SPACE+'T\n'
    Box.append(ADD)
    if len(Box) == 1:
        Box[0]=Box[0][0:-2]+'M\n'
    else:
        Box[0]=Box[0][0:-2]+'H\n'
    Result.append(Box)
    Box = []

for LST in Result:
    for STR in LST:
        outputfile.write(STR)
outputfile.close()

" > intermediate$$/scripts/05_getclustercorstructure.py
python3 intermediate$$/scripts/05_getclustercorstructure.py

#06_get 20kb(default) region information on both side of core proteins.

awk -F'\t' '{print $2}' 05_newdata.txt | uniq > intermediate$$/files/06_justchrom.txt

split -l 150000 intermediate$$/files/06_justchrom.txt intermediate$$/files/06_justchrom_split_  

for x in intermediate$$/files/06_justchrom_split_*; do grep -Fwf $x intermediate$$/files/04_AllInform.txt >>intermediate$$/files/06_ShortAllInform.txt; done  #make Allinform smaller


echo "

def get_protein_surroundings(proteinID, LIST):
    Truelist = []
    for a in LIST:
        b = a.strip().split()
        if len(b) == 7:
            Truelist.append(b)
    for sublist in Truelist:
        if sublist[5] == proteinID:
            ASSEM = sublist[0]
            CHROM = sublist[1]
            START = sublist[2]
            STOP = sublist[3]
            MINI = int(START) - int(STOP)
    box = []
    for sublist in Truelist:
        a = int(START) - int(sublist[3])
        b = int(sublist[2]) - int(STOP)
        if sublist[0] == ASSEM and sublist[1] == CHROM and (($7 >= a >= MINI) or ($7 >= b >= MINI)):
            sublist = ' '.join(sublist) + '\n'
            box.append(sublist)
    return box


# notice that here 'B' was added.
# because it may happen that H is 20k to first B but far from second B (first B is close two second B though)

inputfile1 = open('intermediate$$/files/06_ShortAllInform.txt')
inputfile2 = open('intermediate$$/files/05_newdata.txt')
outputfile = open('intermediate$$/files/06_20Kblist.txt', 'a')
lines1 = inputfile1.readlines()
lines2 = inputfile2.readlines()
inputfile1.close()
inputfile2.close()
LENGTH1 = len(lines1)
LENGTH2 = len(lines2)

chromlist1 = []
for a in lines1:
    b=a.split()[1]
    chromlist1.append(b)
chromlistR1 = chromlist1[:]
chromlistR1.reverse()

LIST1=[]
Result = []
Box = []
Newbox = []
Finalbox = []
COUNT = 1
Num = 0
while Num < LENGTH2:
    ASSEMB, CHROM, START, STOP, STRAND, ID, PROLENTH, TAG = lines2[Num].strip().split()
    a=chromlist1.index(CHROM)
    b=LENGTH1-chromlistR1.index(CHROM)
    LIST1=lines1[a:b]
    Num += 1
    if TAG == 'M':
        ADD = get_protein_surroundings(ID, LIST1)
        for c, s in enumerate(ADD, 1):
            Box.append(str(COUNT) + ' ' + str(c) + ' ' + s)
        for STR in Box:
                outputfile.write(STR)
        COUNT += 1
        Box = []
    elif TAG in 'HB':
        ADD = get_protein_surroundings(ID, LIST1)
        for x in ADD:
            Box.append(x)
    elif TAG == 'T':
        ADD = get_protein_surroundings(ID, LIST1)
        for y in ADD:
            Box.append(y)
        for z in Box:
            if z not in Newbox:
                Newbox.append(z)
        for c, s in enumerate(Newbox, 1):
            Finalbox.append(str(COUNT) + ' ' + str(c) + ' ' + s)
        for STR in Finalbox:
            outputfile.write(STR)
        COUNT += 1
        Box = []
        Newbox = []
        Finalbox = []

#  elif TAG == 'B':
#      pass


outputfile.close()

" > intermediate$$/scripts/06_get20Kb.py
python3 intermediate$$/scripts/06_get20Kb.py

#07_get the annotation of gene clusters.
awk '{print $8}' intermediate$$/files/06_20Kblist.txt | sort -u > intermediate$$/files/07_uniqclusterproteinID.txt
esl-sfetch -f $4 intermediate$$/files/07_uniqclusterproteinID.txt > intermediate$$/files/07_clusterprotein.fa
mkdir intermediate$$/files/07_hmmscansecond/
hmmscan --domtblout intermediate$$/files/07_hmmscansecond/domtscanout.dtbl --noali -E 1e-5 $1 intermediate$$/files/07_clusterprotein.fa >/dev/null
cat intermediate$$/files/07_hmmscansecond/domtscanout.dtbl | grep -v '^#' | awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' | sed 's/ /\t/g' | \
sort -k 3,3 -k 8n -k 9n | \
perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' | \
perl -e 'while(<>){chomp;@a=split(/\t/,$_);if(($a[-1]-$a[-2])>80){print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<=1e-5;}else{print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<=1e-3;}}' | awk '$NF>0.3' | sort -k 3,3 -k 8,8n -k 9,9n | awk -v OFS='\t' '{print $3,$1,$8,$9,$4}' > intermediate$$/files/07_hmmscansecond/parsedscanout

echo "
inputfile = open('intermediate$$/files/07_hmmscansecond/parsedscanout')
outputfile = open('intermediate$$/files/07_hmmscansecond/parsedscanoutinform', 'w')
lines = inputfile.readlines()
inputfile.close()
LENGTH=len(lines)
Result=[]
TAG=''
SPACE='\t'
KS=['ketoacyl-synt','Ketoacyl-synt_C','KAsynt_C_assoc']
AT='Acyl_transf_1'
ACP='PP-binding'
KR='KR'
ER=['ADH_N','ADH_zinc_N','ADH_zinc_N_2']
DH='PS-DH'
DE='PKS_DE'
TE='Thioesterase'
#DOMAINPOOL=['Docking','ketoacyl-synt','Ketoacyl-synt_C','KAsynt_C_assoc','Acyl_transf_1','PP-binding','KR','ADH_N','ADH_zinc_N','ADH_zinc_N_2','PS-DH','PKS_DE','Thioesterase']
#Notice that Docking and PKSDE were ignored
DICT={'Condensation':'C','AMP-binding_C':'A','AMP-binding':'A','Thioesterase':'TE','Docking':'','PKS_DE':'','PS-DH':'DH','ADH_N':'ER','ADH_zinc_N':'ER','ADH_zinc_N_2':'ER','KR':'KR','PP-binding':'ACP','Acyl_transf_1':'AT','ketoacyl-synt':'KS','Ketoacyl-synt_C':'KS','KAsynt_C_assoc':'KS'}
TAG=''
ID, DOMAIN, START, STOP, LEN = lines[0].strip().split()
if DOMAIN in DICT.keys():
    DOMAINABB=DICT[DOMAIN]
else:
    DOMAINABB=DOMAIN
TAG=ID+' '+DOMAINABB
Num = 0
while Num < (LENGTH-1):
    ID, DOMAIN, START, STOP, LEN = lines[Num].strip().split()
    if DOMAIN in DICT.keys():
        DOMAINABB=DICT[DOMAIN]
    else:
        DOMAINABB=DOMAIN
    Num += 1
    ID1, DOMAIN1, START1, STOP1, LEN1 = lines[Num].strip().split()
    if DOMAIN1 in DICT.keys():
        DOMAINABB1=DICT[DOMAIN1]
    else:
        DOMAINABB1=DOMAIN1
    if ID == ID1:
        if DOMAIN != DOMAIN1 and DOMAINABB == DOMAINABB1:
            pass
        else:
            if DOMAINABB1 != '' and not TAG.endswith('\n'):
                TAG=TAG+'_'+DOMAINABB1
            elif DOMAINABB1 != '' and TAG.endswith('\n'):
                TAG=TAG+DOMAINABB1
            else:
                TAG=TAG
    else:
        TAG=TAG+'\n'
        Result.append(TAG)
        if DOMAINABB1 != '':
            TAG=ID1+' '+DOMAINABB1
        else:
            TAG=TAG
        TAG=ID1+' '+DOMAINABB1
if DOMAIN1 in DICT.keys():
    DOMAINABB1=DICT[DOMAIN1]
else:
    DOMAINABB1=DOMAIN1
for STR in Result:
    outputfile.write(STR)
outputfile.close()
" > intermediate$$/scripts/07_moduleanalyzer.py
python3 intermediate$$/scripts/07_moduleanalyzer.py

#08_write the final result
echo "
inputfile1 = open('intermediate$$/files/06_20Kblist.txt')
inputfile2 = open('intermediate$$/files/07_hmmscansecond/parsedscanoutinform')
inputfile3 = open('intermediate$$/files/03_TargetID.txt')
inputfile4 = open('intermediate$$/files/07_uniqclusterproteinID.txt')
outputfile = open('$6', 'w')
lines1 = inputfile1.readlines()
lines2 = inputfile2.readlines()
lines3 = inputfile3.readlines()
lines4 = inputfile4.readlines()
inputfile1.close()
inputfile2.close()
inputfile3.close()
inputfile4.close()

list3=[]
for a in lines3:
    list3.append(a.strip())

pfamdic={}
for a in lines2:
    name, pfam = a.strip().split()
    pfamdic[name]=pfam

for b in lines4:
    if b.strip() not in pfamdic.keys():
        pfamdic[b.strip()]=''

LENGTH1=len(lines1)
Result=[]
Box=[]
SPACE='\t'

for b in lines1:
  TBN, PRN, ASSEMB, CHROM, START, STOP, STRAND, ID, PROLENTH = b.strip().split()
  if ID in list3:
      ADD=TBN+' PRT'+PRN+' '+ASSEMB+' '+CHROM+' '+START+' '+STOP+' '+STRAND+' '+ID+' '+PROLENTH+' '+pfamdic[ID]+' '+'[---------------PKS----------------]\n'
      Box.append(ADD)
  else:
      ADD=TBN+' PRT'+PRN+' '+ASSEMB+' '+CHROM+' '+START+' '+STOP+' '+STRAND+' '+ID+' '+PROLENTH+' '+pfamdic[ID]+'\n'
      Box.append(ADD)

Record=[]
#genomeBGCnum totalnum proteinnum genome chromosome start end strand protein proteinlength conserved domains
COUNT=1
Num=1
Record.append('BGC'+str(Num)+' '+Box[0])
while COUNT < len(Box):
    if Box[COUNT].strip().split()[0] == Record[-1].strip().split()[1] and Box[COUNT].strip().split()[2] == Record[-1].strip().split()[3]:
        Record.append('BGC'+str(Num)+' '+Box[COUNT])
        COUNT+=1
    elif Box[COUNT].strip().split()[0] != Record[-1].strip().split()[1] and Box[COUNT].strip().split()[2] == Record[-1].strip().split()[3]:
        Record.append('\n')
        Num+=1
        Record.append('BGC'+str(Num)+' '+Box[COUNT])
        COUNT+=1

    elif Box[COUNT].strip().split()[0] != Record[-1].strip().split()[1] and Box[COUNT].strip().split()[2] != Record[-1].strip().split()[3]:
        Record.append('\n\n\n\n\n')
        Num=1
        Record.append('BGC'+str(Num)+' '+Box[COUNT])
        COUNT+=1

for STR in Record:
  outputfile.write(STR)
outputfile.close()

" > intermediate$$/scripts/08_from20Kbtofinalresult.py
python3 intermediate$$/scripts/08_from20Kbtofinalresult.py
