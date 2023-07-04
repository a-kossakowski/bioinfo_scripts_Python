from Bio import Entrez
from bioservices import UniProt

#retrives scientific papers data (title, abstract & auythors) linked to specified gene or genes 
#from UniProt, using a set of user-defined keywords (terms)

#email  for authentification
Entrez.email = 'dummy@mail.com'

#list of gene names, UniProt
gene_names = ['Example_Gene']  #replace with gene names

#initialize UniProt
u = UniProt()

gene_ids = []

for gene_name in gene_names:
    #search for gene in UniProt
    result = u.search(gene_name, frmt="tab", columns="id")
    gene_id_list = result.strip().split('\n')
    
    if len(gene_id_list) > 1:
        print(f"Multiple entries found for gene name: {gene_name}")
    elif len(gene_id_list) == 0:
        print(f"No entry found for gene name: {gene_name}")
    else:
        gene_id = gene_id_list[0]
        gene_ids.append((gene_name, gene_id))

publication_details = []

for gene_name, gene_id in gene_ids:
    term = f'{gene_id} AND (cilia OR ciliopathy OR flagella OR motility)'
    handle = Entrez.esearch(db='pubmed', term=term, retmax=1)
    record = Entrez.read(handle)
    if record['IdList']:
        publication_id = record['IdList'][0]
        handle = Entrez.efetch(db='pubmed', id=publication_id, retmode='xml')
        record = Entrez.read(handle)

        #extract details from record
        title = record['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle']
        abstract = record['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText']
        authors = record['PubmedArticle'][0]['MedlineCitation']['Article']['AuthorList']

        publication_details.append((gene_name, title, abstract, authors))
        handle.close()

#display results
for gene_name, title, abstract, authors in publication_details:
    print(f'Gene Name: {gene_name}')
    print(f'Title: {title}')
    print(f'Abstract: {abstract}')
    print(f'Authors: {", ".join(authors)}')
    print('---')
