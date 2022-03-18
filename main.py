from Bio.SubsMat import MatrixInfo
import re                                                                                                   #regular expression
#import Bio.Align.substitution_matrices

#__________________FUNCTIONS__________________

#_REMOVE LOW COMPLEXITY REGIONS AND SEQUENCE REPEATS__

def remove_lowComplexity_seqRepeats(protein_seq_query):
    corresponding_blosumValue = ""
    for amino_acid in protein_seq_query:                                                                    #bshof al score ben al a.a w nfso f blosum
        if(blosum.get((amino_acid, amino_acid)) == 11):
            corresponding_blosumValue += '*'                                                                #3shan al (W,W) ht return 11, f htakhod 2 places fl string, (X,X) return -1 , htakhod brdo two places wahda ll -ve sign
        elif(blosum.get((amino_acid, amino_acid)) == -1):
            corresponding_blosumValue += '@'
        else:
            corresponding_blosumValue += str(blosum.get((amino_acid, amino_acid)))

    print("the corresponding blosum value for each a.a with itself :- ", corresponding_blosumValue)
    repeated_seq= []
    seen = {"0": 0}
    n=2                                                                                                     #hat check al awl atnenat then tlatat.... kol al possibilities
    while(n <= len(corresponding_blosumValue)/2):                                                           #awl iteration hakon bdwr 3la repeated sequences of 2 a.a , second iteration 3la seq of three a.a and so on lhd len/2 , why?? lw aktr mn kda htl3 bara al range bt3 al array asln, akbr repeated seq momkn ykon mwgod at least TWICE al length bt3o hwa nos al length bt3 al seq
        i = 0
        while(i+n <= len(corresponding_blosumValue)):                                                       #general idea any bshof kol al possible substrings mn al string da, b store kol substring in a map(as the key) along with number of occurences(value) w al value lw wslt 2 h3mlo store fl repeated_seq
            subseq = corresponding_blosumValue[i: i+n]
            #print(subseq)
            seen[subseq] = seen.get(subseq, 0) + 1                                                         #lw already mwgod fl map increment l value +1 else set it b zero
            if seen.get(subseq) == 2:
                repeated_seq.append(subseq)
            i = i+1
        n = n+1
    print("repeated sequences are (not consecutive) : - ", repeated_seq)                                 #repeated_seq feha kol al zahar aktr mn marten mlesh d3wa lesa kanp wara b3d wla la
    print("***************")

    for x in repeated_seq:
        if corresponding_blosumValue.count(x) > 1:                                                       #lw al substring da mwgoddd aktr mn maraa adkhol shelo, leh aktr mn mara? msh ana kda kda mot2kda an al substrings de zhrt akr mn mara using al loop ale fo2ya?? 3shan kan momkn ykon fe overlapping between 2 sequences zhro mrten, f lama shlt awl wahed al number of it's occurrence n2s wahed(try : ARNDCQEGHILK)
            i = 0
            n = len(x)
            while i < len(corresponding_blosumValue)-n :

                current = corresponding_blosumValue[i: i + n]

                numOf_consecutive_occurences = 1
                if current == x:
                    while True:
                        next = corresponding_blosumValue[i +numOf_consecutive_occurences * n: i + (numOf_consecutive_occurences+1) * n]
                        if(next != x):
                            break
                        numOf_consecutive_occurences +=1
                if numOf_consecutive_occurences >1:
                    tmp = ""
                    for k in range(0, numOf_consecutive_occurences* n):
                        tmp += "-"
                    corresponding_blosumValue = corresponding_blosumValue[:i] + tmp + corresponding_blosumValue[(i + numOf_consecutive_occurences * n):]
                    i = i+ numOf_consecutive_occurences* n
                else:
                    i = i+1

    print("AFTER REMOVING REPEATS:- ", corresponding_blosumValue)
    protein_seq_query_noRepeats=""
    for i in range (0, len(protein_seq_query)):
        if(corresponding_blosumValue[i] != "-"):
            protein_seq_query_noRepeats+=protein_seq_query[i]

    return protein_seq_query_noRepeats

#_GENERATING THE ORIGINAL WORDS FROM THE ORIGINAL SEQUENCE AFTER REMOVING LCRs__

def generate_originalWords(protein_seq_query_noRepeats, word_length):
    i = 0
    original_words = []
    while (i + word_length <= len(protein_seq_query_noRepeats)):
        temp = protein_seq_query_noRepeats[i: i + word_length]
        original_words.append(temp)
        i= i+1
    return original_words
#_FINDING NEIGHBORING WORDS FOR EACH W-LETTER WORD__


def permute(words):
    letters = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    list1 = []
    for i in words:
        for j in letters:
            list1.append(words.replace(i, j))
            # print(dummy.replace(i, j))
    return list1


def get_seed(threshold, words):
    listseeds = []
    listscore = []
    listd = []
    for i in range(0, len(words)):
        dummy = words[i]
        listd.append(permute(dummy))
        #  print( words[i], " Permutations are:")
        # print(listd[i])
        for g in range(0, len(listd[i])):
            # print('g:',listd[i][g])
            sum = 0
            for letter in range(0, len(listd[i][g])):
                # print('LETTER:',listd[i][g][letter])
                # print(dummy[letter])
                # print(blosum.get((dummy[letter],listd[i][g][letter])))
                if (blosum.get((dummy[letter], listd[i][g][letter])) == None):
                    sum += blosum.get((listd[i][g][letter], dummy[letter]))
                else:
                    sum += blosum.get((dummy[letter], listd[i][g][letter]))
                    # print(sum)
            if (sum >= threshold):
                if (listseeds.__contains__(listd[i][g])):
                    continue
                else:
                    listseeds.append(listd[i][g])
                    listscore.append(sum)
    return listseeds

def HSP_PAIR(LeftQ,RightQ,querySEQ,database,Leftdb,Rightdb,seed):
    currentMax =0

    for s in seed:
        currentMax += blosum.get((s, s))
   
    while (LeftQ > 0 and Leftdb > 0) or (RightQ < len(querySEQ)-1 and Rightdb < len(database)-1):
        if (LeftQ > 0 and Leftdb > 0) and (RightQ < len(querySEQ)-1 and Rightdb < len(database)-1):

            if (blosum.get((querySEQ[LeftQ], database[Leftdb])) == None):
                l = blosum.get((database[Leftdb], querySEQ[LeftQ]))
            else:
                l = blosum.get((querySEQ[LeftQ], database[Leftdb]))
            if (blosum.get((querySEQ[RightQ], database[Rightdb])) == None):
                r = blosum.get((database[Rightdb], querySEQ[RightQ]))
            else:
                r = blosum.get((querySEQ[RightQ], database[Rightdb]))
            sum = l + r

        elif (LeftQ < 0) or (Leftdb < 0):
            if(blosum.get((querySEQ[RightQ], database[Rightdb])) ==  None):
                sum = blosum.get((database[Rightdb],querySEQ[RightQ]))
            else:
                sum = blosum.get((querySEQ[RightQ], database[Rightdb]))

        elif (RightQ > len(querySEQ)-1) or (Rightdb > len(database)-1):

            if (blosum.get((querySEQ[LeftQ], database[LeftQ])) == None):
                sum = blosum.get((database[Leftdb], querySEQ[LeftQ]))
            else:
                sum = blosum.get((querySEQ[LeftQ], database[Leftdb]))

        currentMax += sum        #Getting global max
        #previousMax = currentMax - sum   #Getting local max
        #difference=4
        #if(previousMax-difference>currentMax):
         #   return previousMax

        if (LeftQ != 0):
            LeftQ -= 1
        if(Leftdb != 0):
            Leftdb -= 1
        if (RightQ != len(querySEQ) - 1):
            RightQ += 1
        if (Rightdb != len(database)-1):
            Rightdb += 1



    return currentMax


def ExactAlignment(databasesequence,Qseq,listseeds,T):
    for line in databasesequence:
        for seed in listseeds:
            #print(seed)
            for match in re.finditer(seed, line):
                print("DB:\n","DBline:",databasesequence.index(line),"\nStart index:",match.start(),"\nEnd index:", match.end(),"\nSeed:",seed)
                for c in re.finditer(seed, Qseq):
                    #print(seed, listseeds.index(seed), databasesequence.index(line))
                    HSP_Score = HSP_PAIR(c.start() - 1, c.end() + 1, Qseq, line, match.start() - 1, match.end() + 1,seed)
                    #print("Highest Score Computed is ", HSP_Score)

                    if(HSP_Score > T):
                        print("\nSeed:",seed,"\n Seed start index:",c.start(),"\nSeed End index:", c.end(),"\nDB_ID:",databasesequence.index(line),"\nDB start index:",match.start(),"\nDB end index:", match.end(),"\nHSP_Score:",HSP_Score)

                    # print("Q:\n","Start index:",c.start(),"\nEnd index:", c.end(),"\nSeed:", seed)


#__________________MAIN__________________

blosum = MatrixInfo.blosum62;

#___OPENING/READING/DISPLAYING DATABASE FILE_____

f = open(r"C:\Users\bradl\Desktop\Bio.txt")
DB_sequences = f.readlines()
#for line in DB_sequences:
#    print(line)

protein_seq_query = input("ENTER THE PROTEIN SEQUENCE QUERY :- ")
print("*************\n")
word_length = int(input("ENTER WORD LENGTH :- "))
print("********\n")
word_threshold = int(input("ENTER THE WORD THRESHOLD :- "))
print("**********\n")

#HSP_threshold = int(input("ENTER HSP THRESHOLD :- "))
#print("********\n")

protein_seq_query_noRepeats = remove_lowComplexity_seqRepeats(protein_seq_query)
print("SEQUENCE WITHOUT REPEATS: ", protein_seq_query_noRepeats)
print("********\n")
original_words = generate_originalWords(protein_seq_query_noRepeats, word_length)
print("ORIGINAL WORDS: ", original_words)
print("*****\n")

seeds = get_seed(word_threshold, original_words)
print("SEEDS ARE :  ")
print("*****\n")
print(seeds)
SThreshold = int(input("ENTER SCORE THRESHOLD: "))
ExactAlignment(DB_sequences, protein_seq_query, seeds, SThreshold)