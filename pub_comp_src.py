import time
import logging
from Bio import Entrez
import numpy as np
import pandas as pd
from redcap import Project
from datetime import datetime

### Get access keys from the setup file - config.py
import config
import pub_comp_lib

# ## !!** For DEV
from importlib import reload
# reload(name_of_module)
# ## !!** For DEV

start = time.time()

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

### Get table of publication details from pubmed for pmids
# make dataframe of publications
pubs_frame = pub_comp_lib.summary(pmids, config.ncbi_api, variations)
# add compliant pmc status for publications with a pmcid
pubs_frame['pmc_status'] = np.where(pubs_frame.pmcid == '', '', '1')
# write table
pubs_frame.to_csv('batch_pubmed_frame.csv', index=False)

# change blank values to nan- makes column merging easier
pubs_frame[pubs_frame == ''] = np.nan


###################### PMC Section

# loop batches of ~50 pmids for pmc status and tags check
#pmc_rows = pmc_status.pmc_scrape(pubs_frame.pmid[0:50], variations, config.ncbi_login, config.ncbi_pass)
#pmc_rows = []
#batch = 50
#for x in range(0, len(pubs_frame.pmid[0:101]), batch):
#    pmc_rows.extend(pmc_status.pmc_scrape(pubs_frame.pmid[x:x+batch], variations, config.ncbi_login, config.ncbi_pass))

#pmc_frame = pd.DataFrame(pmc_rows, columns=['pmid', 'pmc_status', 'pmc_tags', 'all_awards'])
#pmc_frame.to_csv('DEV_batch_pmc_status.csv', index=False)
# change blank values to nan- makes column merging easier
#pmc_frame[pmc_frame == ''] = np.nan

# get list of publications with red/grey/yellow pmc status to check on
# nihms status
#check_status = pmc_frame.pmid[pmc_frame['pmc_status'].isin(['2', '3', '4'])]

###################### END PMC Section

###################### Start NIHMS Section
# get list of publications with during current grant cycle with no pmcid to check on
# nihms status
#pubs_frame['pub_date'] = pd.to_datetime(pubs_frame['pub_date'], format='%Y-%m-%d')
#config.start = datetime.strptime(config.start, '%m/%d/%Y')
#check_status = pubs_frame.pmid[(pubs_frame.pub_date > config.start) & (pubs_frame.pmcid.isnull())]


### Check NIHMS status
#nihms_frame = pub_comp_lib.get_nihms(check_status, config.ncbi_login, config.ncbi_pass)
#nihms_frame.to_csv('batch_nihms_status.csv', index=False)

# change blank values to nan- makes column merging easier
#nihms_frame[nihms_frame == ''] = np.nan

### Merge the dataframes for final report
#pub_comp = pd.merge(pubs_frame, pmc_frame, on='pmid', how='outer').merge(nihms_frame, on='pmid', how='outer')
#pub_comp = pd.merge(pubs_frame, nihms_frame, on='pmid', how='outer')


# include nihms ids from all dataframes into a final column
#pub_comp['nihms_id'] = pub_comp['nihmsid_x'].combine_first(pub_comp['nihmsid_y'])
#pub_comp['nihms_id'] = pub_comp['nihmsid_y'].combine_first(pub_comp['nihms_id'])

# include pmc ids from all dataframes into a final column
#pub_comp['pmc_id'] = pub_comp['pmcid_x'].combine_first(pub_comp['pmcid_y'])
#pub_comp['pmc_id'] = pub_comp['pmcid_y'].combine_first(pub_comp['pmc_id'])

# remove columns now that pmc and nihms ids have been merged
#pub_comp = pub_comp.drop(['nihmsid_x', 'nihmsid_y','pmcid_x', 'pmcid_y'], axis=1)
###################### END NIHMS Section

###################### PACM Public Access Compliance Monitor for NIHMS status
# establish the root of the pacm publication url
pacm_root = 'https://www.ncbi.nlm.nih.gov/pmc/utils/pacm/l/'

if config.pacm == 'y':
    driver = pub_comp_lib.pacm_login(config.era_login, config.era_pass)
    time.sleep(3)

    # get list of publications with during current grant cycle with no pmcid to check on
    # nihms status
    pubs_frame['pub_date'] = pd.to_datetime(pubs_frame['pub_date'], format='%Y-%m-%d')
    #config.start = datetime.strptime(config.start, '%m/%d/%Y')
    check_status = pubs_frame.pmid[(pubs_frame.pub_date > config.start) & (pubs_frame.pmcid.isnull())]

    pacm_rows = [pub_comp_lib.parse_pacm(driver, pacm_root, x, variations) for x in check_status]
    pacm_frame = pd.DataFrame(pacm_rows, columns=['pmid', 'nihms_id', 'nihms_status',
                                                  'journal_method', 'files_deposited',
                                                  'initial_approval', 'tagging_complete',
                                                  'final_approval', 'initial_actor',
                                                  'latest_actor', 'pacm_grants'])
    driver.quit()
###################### END PACM Section

pub_comp = pubs_frame.rename(columns={'pmcid':'pmc_id', 'nihmsid': 'nihms_id'})

pub_comp = pd.merge(pub_comp, pacm_frame, on='pmid', how='outer')
# include nihms ids from all dataframes into a final column
pub_comp['nihms_id'] = pub_comp['nihms_id_x'].combine_first(pub_comp['nihms_id_y'])
pub_comp['nihms_id'] = pub_comp['nihms_id_y'].combine_first(pub_comp['nihms_id'])
# remove columns now that pmc and nihms ids have been merged
pub_comp = pub_comp.drop(['nihms_id_x', 'nihms_id_y'], axis=1)

pub_comp.to_csv('batch_comprehensive_status.csv', index=False)

### Update REDCap project if one is being used to track publications
if config.rc_token is not None and config.rc_uri is not None and len(pmids) < 5000:
    pub_comp = pub_comp_lib.RC_update_status(pub_comp)
    success = project.import_records(pub_comp)

print('Publication compliance status update process complete in {0:0.1f} minutes' .format((time.time()-start)/60))
