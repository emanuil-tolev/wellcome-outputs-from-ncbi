import sys
import logging
from urllib2 import URLError

from Bio import Entrez

def main(argv=None):
    if not argv:
        argv = sys.argv

    LOG_FORMAT = '%(asctime)-15s %(message)s'
    logging.basicConfig(level=logging.DEBUG, format=LOG_FORMAT)
    log = logging.getLogger(__name__)
    
    Entrez.email = 'emanuil@cottagelabs.com'
    
    query = "Wellcome[GRNT]"
    
    log.info('Sending this query to NCBI: {0}'.format(query))
   
    try:
        handle = Entrez.esearch(db="pubmed", term=query)
    except URLError:
        log.critical('''NCBI query failed due to an URL Error. Are you
        connected to the internet? (It could be that the EUtils API is
        down or Biopython is generating the wrong URL as well.)''')

    record = Entrez.read(handle)
    
    log.info('NCBI returned {0} results.'.format(record['Count']))
    
    #Out[69]: {u'Count': '70680', u'RetMax': '20', u'IdList': ['23803849', '23802516', '23802441', '23774321', '23766552', '23766329', '23760639', '23750340', '23750059', '23747310', '23735787', '23729657', '23728646', '23728081', '23719378', '23719345', '23708966', '23700959', '23698362', '23698061'], u'TranslationStack': [{u'Count': '70680', u'Field': 'GRNT', u'Term': 'Wellcome[GRNT]', u'Explode': 'N'}, 'GROUP'], u'TranslationSet': [], u'RetStart': '0', u'QueryTranslation': 'Wellcome[GRNT]'}
    
    resultset = {}
    
    for pmid in record['IdList']:
        resultset[pmid] = {}
    
        individual_handle = Entrez.efetch(db='pubmed', retmode='xml',
                id=pmid)
        individual_record = Entrez.read(individual_handle)
    
        if len(individual_record) > 1:
            resultset[pmid]['warning'] = 'NCBI response contains multiple items in the individual record list'
    
        if 'PubmedData' not in individual_record[0]:
            resultset[pmid]['error'] = 'NCBI response did not contain PubmedData key'
    
        if 'ArticleIdList' not in individual_record[0]['PubmedData']:
            resultset[pmid]['error'] = 'NCBI response did not contain the ArticleIdList key in the PubmedData dict'
    
        for identifier in individual_record[0]['PubmedData']['ArticleIdList']:
            pass

if __name__ == '__main__':
    main()
