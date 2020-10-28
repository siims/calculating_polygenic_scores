"""
Use case:
A) get results base on my snp values

    1) read in vcf files
    2) create mapping from rsID to alleles
    3) ask for a search word
    4) request data from snpedia
    5) compile a view

"""

import mwclient

agent = 'MySNPBot. Run by User:Xyz. xyz@foo.com Using mwclient/' + mwclient.__version__
# tokens and secrets are only necessary if your bot will write into SNPedia.
# get your own tokens at http://bots.snpedia.com/index.php/Special:OAuthConsumerRegistration
site = mwclient.Site(('https', 'bots.snpedia.com'), path='/',
                    clients_useragent=agent,
                    consumer_token='???',
                    consumer_secret='???',
                    access_token='???',
                    access_secret='???')


site.search("Alzheimer")

# for i, page in enumerate(site.Categories['Is_a_snp']):
#     print(i, page.name)
