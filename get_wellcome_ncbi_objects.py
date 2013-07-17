import sys
import logging
import copy

from urllib2 import URLError
from httplib import HTTPException

from Bio import Entrez

RESULTS_FILE = 'dois-1000-per-line.txt'

LOG_FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(level=logging.DEBUG, format=LOG_FORMAT)
log = logging.getLogger(__name__)

def fail(msg, original_exception):
    global log
    log.critical(msg)
    raise original_exception

def to_file(filename, s):
    try:
        with open(filename, 'wb') as f:
            f.write(str(s))
    except IOError as e:
        fail(e, '''Could not open results file for writing.
IOError {}.
Data to write {}.'''.format(e, s))

def append_file(filename, s):
    try:
        with open(filename, 'ab') as f:
            f.write(str(s))
    except IOError as e:
        fail(e, '''Could not open results file for appending.
IOError {}.
Data to write {}.'''.format(e, s))


def main(argv=None):
    if not argv:
        argv = sys.argv

    resume_at_result_number = 0
    resuming = False
    if len(argv) > 1:
        if argv[1]:
            resume_at_result_number = argv[1]
            resuming = True
    
    Entrez.email = 'emanuil@cottagelabs.com'
    
    query = "Wellcome[GRNT]"
    
    log.info('Sending this query to NCBI: {0}'.format(query))
    log.info('Starting from result number: {}'.format(resume_at_result_number))
   
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=100000,
                retstart=resume_at_result_number)
    except URLError as e:
        fail(e, '''NCBI query failed due to an URL Error. Are you connected
        to the internet? (It could be that the EUtils API is down or
        Biopython is generating the wrong URL.)''')

    record = Entrez.read(handle, validate=False)
    
    log.info('NCBI holds {} records related to this query.'.format(record['Count']))
    
    results = OAGPrep(RESULTS_FILE, resuming=resuming)
    
    for pmid in record['IdList']:
        try:
            individual_handle = Entrez.efetch(db='pubmed', retmode='xml',
                    id=pmid)
            individual_record = Entrez.read(individual_handle,
                    validate=False)
        except ValueError as e:
            log.warn('''ValueError, Biopython probably couldn\'t parse
            something or the returned XML was invalid. Skipping PMID {}.
            Original error {}'''.format(pmid, e))
            continue
        except (URLError, HTTPException) as e:
            log.warn('''Networking error. Skipping PMID {}.
            Original error {}'''.format(pmid, e))
            continue
    
        if len(individual_record) > 1:
            log.warn('PMID {}: NCBI response contains multiple items in the individual record list'.format(pmid))
    
        if 'PubmedData' not in individual_record[0]:
            log.warn('PMID {}: NCBI response did not contain PubmedData key'.format(pmid))
    
        if 'ArticleIdList' not in individual_record[0]['PubmedData']:
            log.warn('PMID {}: NCBI response did not contain the ArticleIdList key in the PubmedData dict'.format(pmid))
    
        for identifier in individual_record[0]['PubmedData']['ArticleIdList']:
            try:
                if identifier.attributes['IdType'] == 'doi':
                    results.add(identifier)
            except AttributeError:
                log.warn('PMID {}: Can\'t add PMID, no associated DOI.'.format(pmid))

class OAGPrep:
    current_row = []
    rows = [current_row]
    count = 0

    def __init__(self, results_filename, resuming=False):
        self.results_filename = results_filename
        self.do_not_overwrite_files = resuming

        if self.do_not_overwrite_files:
            append_file(self.results_filename, "\n")
        else:
            to_file(self.results_filename, '')

    def add(self, identifier):
        self.count = self.count + 1

        # add commas in front of all items, but skip the comma before
        # the first item of every line
        if self.count % 1000 == 1:
            append_file(self.results_filename, identifier)
        else:
            append_file(self.results_filename, ',' + identifier)

        # add a newline after each set of 1000 items
        if self.count % 1000 == 0:
            append_file(self.results_filename, "\n")
        
        #if len(self.current_row) == 1000:
        #    full_row = copy.copy(self.current_row)
        #    self.rows.insert(0, full_row)
        #    self.current_row = []

        #self.current_row.append(identifier)


    def __str__(self):
        # unroll the rows (lists) into a single string for outputting
        return "\n".join([','.join(row) for row in self.rows])
        

if __name__ == '__main__':
    main()

