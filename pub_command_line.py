import sys # needed to accept args from command line
import time
import logging
from datetime import datetime

import config
import pub_comp_lib

########## Initial Design Notes
# two args: database and time 
# database default is 'all', length check to throw error if more than two sent
# time default is 'all', now uses config.start for current grant cycle
# help shows the list of accepted args
# database list: all, pubmed, pmc, nihms, icite, altmetric
# time list: all, now
###########

### Dev - Check the Params that are passed from the command line
# total arguments
#n = len(sys.argv)
#print("Total arguments passed:", n)
# Arguments passed
#print("\nName of Python script:", sys.argv[0])
#print("\nArguments passed:", end = " ")
#for i in range(1, n):
#    print(sys.argv[i], end = " ")
### End Dev check


start_time = time.time()

# set delay values for pausing during web scraping actions
delay = 2
long_delay = 5

# initiatlize logger
logging.basicConfig(
    filename="test.log",
    level=logging.DEBUG,
    format="%(asctime)s:%(levelname)s:%(message)s"
    )

logger = logging.getLogger(__name__)

# loop over all config grants for cleanup
for x in range(len(config.grant_list)):
    # remove all whitespace, leading or trailing hyphenates - clean_grant.py
    config.grant_list[x] = pub_comp_lib.clean(config.grant_list[x])

### Create list for each grant with 34 grant variations - grant_vari.py
variations = []
for grant in config.grant_list:
    variations.extend(pub_comp_lib.variety(grant))

### check and process arguments passed at command line when script was called
params = pub_comp_lib.check_argv(sys.argv, config.start)

### assign checked command line parameters to variables
db = params[0]
timeframe = params[1]

### based on command line arguments, query databases using pub_comp_lib functions
if 'pubmed' in db:
    pubmed_start = time.time()
    print('PubMed started at '+ str(datetime.now()))
    pub_comp_lib.query_pubmed(logger, variations, config.ncbi_api, config.rc_uri, 
    							config.rc_token)
    print('PubMed query complete in {0:0.1f} minutes' .format((time.time()-pubmed_start)/60))

if 'pmc' in db:
    pmc_start = time.time()
    print('PMC started at '+ str(datetime.now()))
    pub_comp_lib.query_pmc(logger, timeframe, variations, delay, long_delay, config.ncbi_login, 
    						config.ncbi_pass, config.rc_uri, config.rc_token)
    print('PMC query complete in {0:0.1f} minutes' .format((time.time()-pmc_start)/60))
    
if 'nihms' in db:
    nihms_start = time.time()
    print('NIHMS started at '+ str(datetime.now()))
    pub_comp_lib.query_nihms(logger, timeframe, delay, long_delay, config.ncbi_login, 
    						config.ncbi_pass, config.rc_uri, config.rc_token)
    print('NIHMS query complete in {0:0.1f} minutes' .format((time.time()-nihms_start)/60))
    
if 'icite' in db:
    icite_start = time.time()
    print('iCite started at '+ str(datetime.now()))
    pub_comp_lib.query_icite(logger, timeframe, config.rc_uri, config.rc_token)
    print('iCite query complete in {0:0.1f} minutes' .format((time.time()-icite_start)/60))

if 'altmetric' in db:
    altmetric_start = time.time()
    print('Altmetric started at '+ str(datetime.now()))
    pub_comp_lib.query_altmetric(logger, timeframe, config.rc_uri, config.rc_token)
    print('Altmetric query complete in {0:0.1f} minutes' .format((time.time()-altmetric_start)/60))

if len(sys.argv) >= 2 and sys.argv[1] == 'compliance':
	load_compliance_start = time.time()
	print('Loading non-compliant publications started at '+ str(datetime.now()))
	pub_comp_lib.pmc_add_non_compliant(logger, timeframe, variations, delay, long_delay, config.ncbi_login, 
    						config.ncbi_pass, config.rc_uri, config.rc_token)
	print('Loading non-compliant publications complete in {0:0.1f} minutes' .format((time.time()-load_compliance_start)/60))
    
print('All queries complete in {0:0.1f} minutes' .format((time.time()-start_time)/60))
