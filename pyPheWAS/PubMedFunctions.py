def search(query,rs,f):
    Entrez.email = dev_email
    handle = Entrez.esearch(db='pubmed',
                           sort='relevance',
                           retstart=rs,
                           retmax=10000,
                           retmode='xml',
                           field=f,
                           # usehistory='y',
                           term=query)

    results = Entrez.read(handle)
    handle.close()
    return results