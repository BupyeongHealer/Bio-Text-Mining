from Bio import Entrez
import math
from collections import deque

def probability(abstracts, keywords):
    p = 0.0
    count = 0.0
    total = len(abstracts)

    for words in abstracts:
        has_terms = True
        for t in keywords:
            if not t in words:
                has_terms = False
        if has_terms:
            count = count + 1

    return count / total

kegg_names = {}
name_kegg = {}
f = open('C:\\Users\wonde\Desktop\genes.txt', 'r')
for line in f.readlines():

    t1 = line.split(';')[0]
    t2 = t1.split('\t')
    kegg_id = t2[0]
    kegg_names[kegg_id] = []
    for name in t2[1].split(','):
                name = name.strip()
                kegg_names[kegg_id].append(name)
                name_kegg[name] = kegg_id

f.close()

disease = 'Schizophrenia'.upper()

print('Download abstracts...')

Entrez.email = 'tkddudwhswk@naver.com'

handle = Entrez.esearch(db='pubmed', term=disease, retmax=10000)
record = Entrez.read(handle)

downloaded_abstracts = []

cnt = 0
for pubmed_id in record['IdList']:
    cnt = cnt + 1

    print (cnt, '/', len(record['IdList']))
    abstract = Entrez.efetch('pubmed', id=pubmed_id, retmode='text', rettype='abstracts').read()
    downloaded_abstracts.append(abstract)

keywords_in_abstract = []
for ab in downloaded_abstracts:
    keyword_box = []
    words = ab.replace('.',' ').split(' ')
    for w in words:
        if w.upper() == disease:
            keyword_box.append(w.upper())
        else:
            if w in name_kegg:
                keyword_box.append(name_kegg[w])

    keywords_in_abstract.append(keyword_box)

print('Calculating MI....')

scores = {}

p_disease = probability(keywords_in_abstract, [disease])       #dent로 하면 결과 달라짐
for kegg_id in kegg_names:
    p_gene = probability(keywords_in_abstract, [kegg_id])
    p_gene_disease = probability(keywords_in_abstract, [kegg_id, disease])

    if p_gene !=0 and p_disease !=0 and p_gene_disease !=0:
        mi = math.log2(p_gene_disease / (p_gene * p_disease))
        scores[  kegg_names[kegg_id][0]  ] = mi

f2 = open('C:\\Users\wonde\Desktop\\result.txt','w')

for key in sorted(scores, key=scores.__getitem__, reverse=True):
    f2.write(key + '\t' + str(scores[key]) + '\n')
    print(key + '\t' + str(scores[key]))
f2.close()
