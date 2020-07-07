## Text mining from PubMed abstracts

A research purpose is to extract Thalassemia associated genes from abstracts through text mining strategy.  

### Development Environment
- IDE : Pycharm
- Package : Biopython
- Python Version : 3.6 & 3.7

### Installation

1. We need to install `biopython` package.  
```cmd
pip install biopython
```
2. Update old `biopython` package.
```cmd
pip install biopython --upgrade
```

## Overview

[0. Settings](#setting)  
[1. Load the Kegg genes file](#load)  
[2. Download the abstracts from PubMed](#download)  
[3. Text mining the abstracts](#textmining)  
[4. Scoring the selected words](#scoring)  
[5. Discussion](#discussion)  

### 0. Settings <a id="setting"></a>
```python
from Bio import Entrez
import math
import time
import random

time.sleep(random.randint(1, 3))

```

### 1. Load the Kegg genes file <a id="load"></a>
```python
kegg_names = {}
name_kegg = {}

f = open('C:\\Users\PARK\\genes.txt', 'r')
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
```

### 2. Download the abstracts from Pubmed <a id="download"></a>

``` python
disease = 'Thalassemia'.upper()

print('Download abstracts...')

Entrez.email = 'tkddudwhswk@naver.com'

# noinspection PyInterpreter
handle = Entrez.esearch(db='pubmed', term=disease, retmax=10000)
record = Entrez.read(handle)

downloaded_abstracts = []

cnt = 0
for pubmed_id in record['IdList']:
    cnt = cnt + 1

    print(cnt, '/', len(record['IdList']))
    abstract = Entrez.efetch('pubmed', id=pubmed_id, retmode='text', rettype='abstracts').read()
    downloaded_abstracts.append(abstract)
```

### 3. Text mining the abstracts <a id="textmining"></a>
``` python
keywords_in_abstract = []
for ab in downloaded_abstracts:
    keyword_box = []
    words = ab.replace('.', ' ').split(' ')
    for w in words:
        if w.upper() == disease:
            keyword_box.append(w.upper())
        else:
            if w in name_kegg:
                keyword_box.append(name_kegg[w])

    keywords_in_abstract.append(keyword_box)
```

### 4. Scoring the selected words <a id="scoring"></a>

``` python
def probability(abstracts, keywords):
    count = 0.0
    total = len(abstracts)

    for WORDS in abstracts:
        has_terms = True
        for t in keywords:
            if not t in WORDS:
                has_terms = False
        if has_terms:
            count = count + 1

    return count / total

print('Calculating MI....')

scores = {}

p_disease = probability(keywords_in_abstract, [disease])
for kegg_id in kegg_names:
    p_gene = probability(keywords_in_abstract, [kegg_id])
    p_gene_disease = probability(keywords_in_abstract, [kegg_id, disease])

    if p_gene != 0 and p_disease != 0 and p_gene_disease != 0:
        mi = math.log2(p_gene_disease / (p_gene * p_disease))
        scores[kegg_names[kegg_id][0]] = mi

f2 = open('C:\\Users\PARK\\result.txt', 'w')

for key in sorted(scores, key=scores.__getitem__, reverse=True):
    f2.write(key + '\t' + str(scores[key]) + '\n')
    print(key + '\t' + str(scores[key]))
f2.close()

```

### 5. Discussion <a id="discussion"></a>
741 different genes are collected by my text mining codes. 541 genes got a plus score and 200 genes got a minus score.  
A Plus score means that It is more likely to exist genes and disease at same time in abstracts than alone.  
But a Minus score means that the minus scored genes will be likely not to exist together. They will exist alone in abstracts.  
If the genes and disease exist in abstracts simultaneously, we can draw the conclusion that there may be a significant correlation between a disease and genes
