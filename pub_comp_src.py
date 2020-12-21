import time
import logging
from Bio import Entrez
import numpy as np
import pandas as pd
from redcap import Project
from datetime import datetime
import csv

from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys

### Get access keys from the setup file - config.py
import config
import pub_comp_lib

# ## !!** For DEV
from importlib import reload
# reload(name_of_module)
# ## !!** For DEV

start_time = time.time()

logging.basicConfig(
    filename="test.log",
    level=logging.DEBUG,
    format="%(asctime)s:%(levelname)s:%(message)s"
    )

logger = logging.getLogger(__name__)
#logger = logging.basicConfig(filename='app.log', filemode='w',
#                             format='%(name)s - %(levelname)s - %(message)s')

# loop over all config grants for cleanup
for x in range(len(config.grant_list)):
    # remove all whitespace, leading or trailing hyphenates - clean_grant.py
    config.grant_list[x] = pub_comp_lib.clean(config.grant_list[x])

### Create list for each grant with 34 grant variations - grant_vari.py
variations = []
for grant in config.grant_list:
    variations.extend(pub_comp_lib.variety(grant))

### Get pmids from pubmed for all grant variations
# create variables for pubmed queries
Entrez.email = "Your.Name.Here@example.org"
Entrez.api_key = config.ncbi_api

# create set for unique list of all pmids from querying pubmed with each
# grant variation
pmids = set()
# query pubmed for pmids associated with each grant variation
logger.info("Starting pubmed queries...")
pubmed_results = []

for grant in variations:
    attempt = 1
    while attempt <= 3:
        try:
            handle = Entrez.esearch(db='pubmed', term=grant,
                                    field='grant', retmax=5000,
                                    usehistory='y', retmode='xml')
            record = Entrez.read(handle)
            handle.close()
            if int(record['Count']) > 0:
                pubmed_results.append(record)
                pmids.update(record['IdList'])
                logger.info('Entrez ESearch returns %i Ids for %s' % (int(record['Count']), str(grant)))
            attempt = 4
        except Exception as err:
            logger.warning('Received error from server: %s' % str(err))
            logger.warning('Attempt %i of 3 for grant %s.' % (attempt,
                                                              str(grant)))
            attempt += 1
            time.sleep(2)
    logger.debug('Grant %s queried.' % str(grant))

logger.info('All grant queries complete.')

##### To test for PubMed downtime or blocked access...
#handle = Entrez.esearch(db='pubmed', term=grant[0], field='grant')
#record = Entrez.read(handle)
#handle.close()
#print(record)

### Update pmid set if a REDCap project is being used to track publications
if config.rc_token is not None and config.rc_uri is not None:
    old_pmids = []
    # get the full pmid list from the REDCap project
    project = Project(config.rc_uri, config.rc_token)
    rc_pmids = project.export_records(fields=['pmid'], format='json')
    for rc_pmid in rc_pmids:
        old_pmids.append(rc_pmid['pmid'])
    new_pmids = list(pmids.difference(old_pmids))   # newly discovered pmids
    pmids.update(old_pmids)
    # date of first discovery
    if len(new_pmids) > 0:
        first_disc = [datetime.today().strftime("%Y-%m-%d")]*len(new_pmids)
        # create data frame of new_pmids with date of first dicovery and
        # import into REDCap project
        # create data frame using lists and import into redcap
        first_discovered_frame = pd.DataFrame(np.column_stack([new_pmids, first_disc]),
                            columns=['pmid', 'first_discovered'])
        response = project.import_records(first_discovered_frame)

###################### PubMed Summary Section
### Get table of publication details from pubmed for pmids
# make dataframe of publications
pubs_frame = pub_comp_lib.summary(pmids, config.ncbi_api, variations)
# add compliant pmc status for publications with a pmcid
pubs_frame['pmc_status'] = np.where(pubs_frame.pmcid == '', '', '1')
# write table
pubs_frame.to_csv('batch_pubmed_frame.csv', index=False)

# change blank values to nan- makes column merging easier
pubs_frame[pubs_frame == ''] = np.nan

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pubs_frame = pubs_frame.rename(columns={'pmcid':'pmc_id', 'nihmsid':'nihms_id'})
###################### END PubMed Summary Section


###################### PMC Section

#####  For updated PMC interface
# log into era commons
attempt = 1
while attempt <= 3:
    try:
        driver = pub_comp_lib.ncbi_login(config.ncbi_login, config.ncbi_pass)
        attempt = 4
    except Exception as err:
        logger.warning('Unable to log into ERA Commons, attempt %i; error: %s' % (attempt, str(err)))
        attempt += 1
        time.sleep(2)

# get list of publications with during current grant cycle with no pmcid to check on
# nihms status
pubs_frame['pub_date'] = pd.to_datetime(pubs_frame['pub_date'], format='%Y-%m-%d')
#config.start = datetime.strptime(config.start, '%m/%d/%Y')

#!!!!!!! how much of the pubmed results are going to pmc to check for compliance
status_pmc = list(pubs_frame.pmid[(pubs_frame.pub_date > config.start) & (pubs_frame.pmc_id.isnull())])
#status_pmc = list(pubs_frame.pmid[(pubs_frame.pub_date > config.start) | (pubs_frame.pmc_id.isnull())])
#status_pmc = pubs_frame.pmid[pubs_frame.pmc_id.isnull()]
#status_pmc = list(pubs_frame.pmid[(pubs_frame.pub_date > config.start)])
#status_pmc = list(pubs_frame.pmid)

##!!!!!!!!! DEV ONLY csv file since I can't tell if all pmids are being sent to pmc
#status_pmc.to_csv('pmids_to_check_in_pmc.csv', index=False)
with open("pmids_to_check_in_pmc.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerows(status_pmc)

print("length of status_pmc: " + str(len(status_pmc)))

####################### scrape pmc information in batches
pmc_rows = []
batch_size = 250
count = len(status_pmc)
delay = 2
long_delay = 5

for start in range(0, count, batch_size):
    end = min(count, start+batch_size)
    # reload my bib, clear all publications and load pmids in status_pmc
    print('Starting new batch: ' + str(start) + '-' + str(end) + ', Currently have gathered ' + str(len(pmc_rows)) +' rows.')
    driver.get('https://www.ncbi.nlm.nih.gov/myncbi/collections/mybibliography/')
    time.sleep(long_delay)
    pub_comp_lib.clear_my_bib(driver, delay, logger)
    print('*Cleared MyBib')
    time.sleep(long_delay)
    pub_comp_lib.add_to_my_bib(driver, status_pmc[start:end], delay, long_delay, logger)
    print('**Added ' + str(end-start) + ' pubs to MyBib: ' + str(start) + '-' + str(end) + ' pmids ' + str(status_pmc[start]) + '-' + str(status_pmc[end-1]))
    #print("Start: " + str(start))
    #print("End: " + str(end))
    #print("First PMID: " + status_pmc[0])
    #print("Pubs_Frame PMID: " + pubs_frame.pmid[1])

    # reload my bib and begin scraping each page of citations
    time.sleep(long_delay)
    driver.get('https://www.ncbi.nlm.nih.gov/myncbi/collections/mybibliography/')
    time.sleep(long_delay)

    scrape_more = True
    #### loop for each 'next page' click
    while scrape_more == True:
        soup = BeautifulSoup(driver.page_source, 'lxml')
        cites = soup.find_all('div', 'citation-wrap')
        print('**!Gonna scrape ' + str(len(cites)) + ' citations now.')
        for x in range(len(cites)):
            pmc_rows.append(pub_comp_lib.scrape_citations(cites[x], x, variations, driver, delay, long_delay, logger, start))
        print('**!!Finished a page of scraping citations... got ' + str(len(pmc_rows)) + ' rows.')
        ## check if there's another page of citations to scrape
        time.sleep(delay)
        try:
            next_button = driver.find_element_by_xpath('//*[@id="pager1"]/ul/li[4]/a').get_attribute('onclick')
        except Exception as err:
            next_button = 'return false;'
        if next_button == 'return false;' or driver.find_element_by_xpath('//*[@id="pager2"]/ul/li/span').get_attribute('innerText') == '1':
            scrape_more = False
            print('**!!looks like no next button on publication number:' + str(end) + ' with ' + str(len(pmc_rows)) + ' rows and last pmid logged is ' + str(pmc_rows[-1:]))
        else: driver.find_element_by_xpath('//*[@id="pager1"]/ul/li[4]/a').click()
    #if count > 1000:
    time.sleep(15)
    print('***Completed ' + str(start) + ' : ' + str(end) + '    minutes-{0:0.1f}' .format((time.time()-start_time)/60))
    #time.sleep(delay)
driver.close()

print('All done with scraping.  Got ' + str(len(pmc_rows)) + ' rows and expected ' + str(len(status_pmc)) + ' rows.')

## package the pmc_rows into a data frame
pmc_frame = pd.DataFrame(pmc_rows, columns=['pmid', 'pmc_status', 'pmc_tags', 'all_awards', 'pub_num'])
pmc_frame.to_csv('DEV_batch_pmc_status.csv', index=False)
# change blank values to nan- makes column merging easier
pmc_frame[pmc_frame == ''] = np.nan

# drop the pub_num column after the data frame has been written to csv file
pmc_frame = pmc_frame.drop('pub_num', 1)

# get list of publications with non-compliant pmc status to check on
# nihms status
status_nihms = pmc_frame.pmid[pmc_frame['pmc_status'].isin(['2', '3', '4', ''])]
###################### END PMC Section

############# NEW NIHMS Section

nihms_frame = pub_comp_lib.get_nihms(status_nihms, config.ncbi_login, config.ncbi_pass, 1, 5)

nihms_frame.to_csv('DEV_batch_nihms_status.csv', index=False)
# change blank values to nan- makes column merging easier
nihms_frame[pmc_frame == ''] = np.nan

################# END NEW NIHMS Section

###################### PACM Public Access Compliance Monitor for NIHMS status
# establish the root of the pacm publication url
#pacm_root = 'https://www.ncbi.nlm.nih.gov/pmc/utils/pacm/l/'
#
#if config.pacm == 'y':
#    driver = pub_comp_lib.pacm_login(config.era_login, config.era_pass)
#    time.sleep(3)

    # get list of publications with during current grant cycle with no pmcid to check on
    # nihms status
#    pubs_frame['pub_date'] = pd.to_datetime(pubs_frame['pub_date'], format='%Y-%m-%d')
    #config.start = datetime.strptime(config.start, '%m/%d/%Y')
#    check_status = pubs_frame.pmid[(pubs_frame.pub_date > config.start) & (pubs_frame.pmcid.isnull())]

#    pacm_rows = [pub_comp_lib.parse_pacm(driver, pacm_root, x, variations) for x in check_status]
#    pacm_frame = pd.DataFrame(pacm_rows, columns=['pmid', 'nihms_id', 'nihms_status',
#                                                  'journal_method', 'files_deposited',
#                                                  'initial_approval', 'tagging_complete',
#                                                  'final_approval', 'initial_actor',
#                                                  'latest_actor', 'pacm_grants'])
#    driver.quit()
###################### END PACM Section

########## join pmids, pmc, and nihms tables and upload into REDCap
pub_comp = pd.merge(pubs_frame, pmc_frame, on='pmid', how='outer')
pub_comp = pd.merge(pub_comp, nihms_frame, on='pmid', how='outer')


# include nihms ids from all dataframes into a final column
pub_comp['nihms_id'] = pub_comp['nihms_id_x'].combine_first(pub_comp['nihms_id_y'])
pub_comp['nihms_id'] = pub_comp['nihms_id_y'].combine_first(pub_comp['nihms_id'])

# include pmc ids from all dataframes into a final column
pub_comp['pmc_id'] = pub_comp['pmc_id_x'].combine_first(pub_comp['pmc_id_y'])
pub_comp['pmc_id'] = pub_comp['pmc_id_y'].combine_first(pub_comp['pmc_id'])

# include pmc status from all dataframes into a final column
pub_comp['pmc_status'] = pub_comp['pmc_status_x'].combine_first(pub_comp['pmc_status_y'])
pub_comp['pmc_status'] = pub_comp['pmc_status_y'].combine_first(pub_comp['pmc_status'])

# remove columns now that pmc and nihms ids have been merged
pub_comp = pub_comp.drop(['nihms_id_x', 'nihms_id_y'], axis=1)
pub_comp = pub_comp.drop(['pmc_id_x', 'pmc_id_y'], axis=1)
pub_comp = pub_comp.drop(['pmc_status_x', 'pmc_status_y'], axis=1)

pub_comp['nihms_comm'] = ''

### Update REDCap project if one is being used to track publications
if config.rc_token is not None and config.rc_uri is not None and len(pmids) < 5000:
    pub_comp = pub_comp_lib.RC_update_status(pub_comp)
    success = project.import_records(pub_comp)

# write a copy to a .csv file
pub_comp.to_csv('batch_comprehensive_status.csv', index=False)

print('Publication compliance status update process complete in {0:0.1f} minutes' .format((time.time()-start_time)/60))
