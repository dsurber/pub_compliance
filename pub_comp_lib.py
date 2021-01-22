from Bio import Entrez
from Bio.Entrez import efetch
from Bio.Entrez import read
import regex as re
import datetime
import time
import logging
import pandas as pd
import time
from bs4 import BeautifulSoup
import unicodedata

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys

def clean(grant):
	#	remove all whitespace
	grant = "".join(grant.split())

	#	check for a hyphenated lead character set and remove
	if re.search('.*?-(.*$)', grant) is not None:
		grant = re.search('.*?-(.*$)', grant).group(1)

	#	check for a hyphenated tail character set and remove
	if re.search('(^.*)-.*?$', grant) is not None:
		grant = re.search('(^.*)-.*?$', grant).group(1)

	return grant


# takes a grant number and creates 10 variations
def variety(grant):
	vari = []

	if len(grant) == 11:
		# create sections of grant number for formatting
		first = grant[:3]
		mid = grant[3:5]
		last = grant[5:]

		# prefix and suffix variations
		prefixes = ['', '1', '5', '9']
		suffixes = ['', '-01']

		# outer loop for suffix variations
		for suffix in suffixes:
			# 1LL######
			vari.append('1' + mid + last + suffix)
			# inner loop for prefix variations
			for prefix in prefixes:
				# core variations
				# LL#LL######
				vari.append(prefix + first + mid + last + suffix)
				# LL#_LL######
				vari.append(prefix + first + ' ' + mid + last + suffix)
				# LL#1LL######
				vari.append(prefix + first + '1' + mid + last + suffix)
				# LL#_LL_######
				vari.append(prefix + first + ' ' + mid + ' ' + last + suffix)
	else: vari.append(grant)

	return vari


def details(pub, variations):
    # remove all white space and \n to help regex function
    pub = ''.join(pub.split('\n'))

    pmid = re.search('<PMID.*?>(.*?)</PMID>', pub).group(1)
    # if pmc exists
    if re.search('pmc\">PMC(.*?)</ArticleId>', pub) is not None:
        pmcid = re.search('pmc\">PMC(.*?)</ArticleId>', pub).group(1)
    else:
        pmcid = ''
    # if nihms exists
    if re.search('mid\">NIHMS(.*?)</ArticleId>', pub) is not None:
        nihmsid = re.search('mid\">NIHMS(.*?)</ArticleId>', pub).group(1)
    else:
        nihmsid = ''
    # if nctid exists
    nctid = []
    if re.search('NCT(.*?)</AccessionNumber>', pub) is not None:
        nctid = re.findall('<AccessionNum.*?(NCT[0-9].*?)</AccessionNum', pub)

    nctid = ', '.join(nctid)

    if re.search('<ArticleTitle>(.*?)</ArticleTitle>', pub) is not None:
        pub_title = re.search('<ArticleTitle>(.*?)</ArticleTitle>', pub).group(1)
    else:
        pub_title = ''

    ## loop to get all author info
    # initialize lists
    authors_lnames = []
    authors_initials = []
    authors_fnames = []
    authors_affil = []
    authors_orcid = []

    # split into xml batches of author info
    author_list = re.split('<Author Valid', pub)

    # loop through to get author info
    for x in range(1,len(author_list)):
        if re.search('<LastName>(.*?)</LastName>', author_list[x]) is not None:
            authors_lnames.append(re.search('<LastName>(.*?)</LastName>',
                                            author_list[x]).group(1))
        else:
            authors_lnames.append('Unknown')
        if re.search('<Initials>(.*?)</Initials>', author_list[x]) is not None:
            authors_initials.append(re.search('<Initials>(.*?)</Initials>',
                                              author_list[x]).group(1))
        else:
            authors_initials.append('Unknown')
        if re.search('<ForeName>(.*?)</ForeName>', author_list[x]) is not None:
            authors_fnames.append(re.search('<ForeName>(.*?)</ForeName>',
                                            author_list[x]).group(1))
        else:
            authors_fnames.append('Unknown')
        if re.search('Affiliation>', author_list[x]) is not None:
            authors_affil.append(re.search('Affiliation>(.*?)</Affiliation',
                                           author_list[x]).group(1))
        else:
            authors_affil.append('')
        if re.search('Identifier Source="ORCID">', author_list[x]) is not None:
            authors_orcid.append(re.search('Identifier Source="ORCID">(.*?)</Identifier>',
                                            author_list[x]).group(1))
        else:
            authors_orcid.append('')
    # combine fname and lname to get full list of author names
    authors = [i+' '+j for i, j in zip(authors_fnames, authors_lnames)]

    authors_lnames = ', '.join(authors_lnames)
    authors_fnames = ', '.join(authors_fnames)
    authors_initials = ', '.join(authors_initials)
    authors_affil = ', '.join(authors_affil)
    authors_orcid = ', '.join(authors_orcid)
    authors = ', '.join(authors)


    ## get pub_date from when journal was published
    if re.search('<JournalIssue.*?<PubDate>(.*?)</PubDate>', pub) is not None:
        publish_date = re.search('<JournalIssue.*?<PubDate>(.*?)</PubDate>',
                                 pub).group(1)
    else:
        publish_date = ''

    # clean up pub_date
    if re.search('<Medline', publish_date) is None:
        if re.search('<Year>[0-9]{4}', publish_date) is None:
            year = '2099'
        else:
            year = re.search('<Year>([0-9]{4})</Year>', publish_date).group(1)
        if re.search('<Month>', publish_date) is None:
            month = '01'
        elif re.search('<Month>([A-Za-z].*?)-.*</Month>',
                       publish_date) is not None:
            month = re.search('<Month>([A-Za-z].*?)-.*</Month>',
                              publish_date).group(1)
        else:
            month = re.search('<Month>(.*)</Month>', publish_date).group(1)
        if re.search('<Day>', publish_date) is None:
            day = '01'
        else:
            day = re.search('<Day>(.*)</Day>', publish_date).group(1)
    else:
        medline = re.search('<Medline.*?>(.*?)</Medline.*?>',
                            publish_date).group(1)
        if re.search('[A-Za-z]', medline) is not None:
            year = re.search('^.*?([0-9]{4}).*?$', medline).group(1)
            month = re.search('^.*?([A-Za-z]{3}).*?[-|/].*$', medline).group(1)
            day = '01'
        else:
            year = re.search('^([0-9]{4}).*?$', medline).group(1)
            month = '01'
            day = '01'

    # remove all whitespace from 'month'
    month = ''.join(month.split())

    # combine month day and year
    pub_date = year+ '-' + month + '-' + day
    # check if month is letters and change into numerical date type
    if re.search('[a-zA-Z]', pub_date) is not None:
        pub_date = datetime.datetime.strptime(pub_date,
                                              "%Y-%b-%d").strftime("%Y-%m-%d")

    if re.search('</JournalIssue>.*?<ISOAbbreviation>(.*?)</ISOAbbre',
                              pub) is not None:
        journal_short = re.search('</JournalIssue>.*?<ISOAbbreviation>(.*?)</ISOAbbre',
                                  pub).group(1)
    else:
        journal_short = 'Unknown'

    if re.search('/JournalIssue>.*?<Title>(.*?)</Title>',
                             pub) is not None:
        journal_full = re.search('/JournalIssue>.*?<Title>(.*?)</Title>',
                                 pub).group(1)
    else:
        journal_full = 'Unknown'

    ## get grant list to clean up and compare with variations to get pubmed tags
    pubmed_tags = []
    grant_list = re.split('<Grant>', pub)
    for x in range(1, len(grant_list)):
        if re.search('<GrantID>', grant_list[x]) is not None:
            if re.search('<GrantID>(.*?)</GrantID>',
                         grant_list[x]).group(1) in variations:
                pubmed_tags.append(re.search('<GrantID>(.*?)</GrantID>',
                                             grant_list[x]).group(1))

    pubmed_tags = ', '.join(pubmed_tags)

    ## get publication types to exclude some pubs from NIH PA Policy
    exclude = ''
    pub_types = []
    type_list = re.split('<PublicationType UI', pub)
    for x in range(1, len(type_list)):
        if re.search('\\">(.*?)</PublicationType>', type_list[x]) is not None:
            pub_types.append(re.search('\\">(.*?)</PublicationType>', type_list[x]).group(1))
            if re.search('\\">(.*?)</PublicationType>', type_list[x]).group(1).lower() in ['letter', 'comment', 'editorial']:
                exclude = '1'
    pub_type_list = ', '.join(pub_types)

    ## get mesh heading major and minor topics with qualifiers
    minor_topics = []
    major_topics = []
    key_topics = []

    mesh_list = re.split('<MeshHeading>', pub)
    for x in range(1, len(mesh_list)):
        major = ''
        minor = ''
        # decide if the descriptor is major or minor and extract descriptor text
        if re.search('="N".*?>(.*?)</D', mesh_list[x]) is not None:
            minor = re.search('>(.*?)</D', mesh_list[x]).group(1)
        elif re.search('="Y".*?>(.*?)</D', mesh_list[x]) is not None:
            major = re.search('>(.*?)</D', mesh_list[x]).group(1)
        # check if major or minor descriptor is one of the key topics, if so, append to key list
        if re.search('">(.*?)</D', mesh_list[x]).group(1).lower() in ['pediatrics', 'translational medical research']:
            key_topics.append(re.search('">(.*?)</D', mesh_list[x]).group(1))
        # create list of qualifiers and assemble extracted qualifier text into a list
        qual_list = re.split('<QualifierName', mesh_list[x])
        qualifier = []
        for y in range(1, len(qual_list)):
            if re.search('">(.*?)</Q', qual_list[y]) is not None:
                qualifier.append(re.search('">(.*?)</Q', qual_list[y]).group(1))
        # add qualifiers to the major or minor descriptor and append to the list
        if len(major) > 0:
            major_topics.append(major + ' (' + '; '.join(qualifier) + ')')
        else: minor_topics.append(minor + ' (' + '; '.join(qualifier) + ')')

    mesh_minor = '; '.join(minor_topics)
    mesh_major = '; '.join(major_topics)
    mesh_key = '; '.join(key_topics)

    ## assemble all values for the row of the dataframe to be returned
    row = [pmid, pmcid, nihmsid,  nctid, pub_title, authors,
            authors_lnames, authors_initials, authors_orcid, authors_affil,
            pub_date, journal_short, journal_full, pubmed_tags, pub_type_list,
            exclude, mesh_major, mesh_minor, mesh_key]

    return row


def summary(pmids, ncbi_key, grants):

	#***!!! developing !!!***
	Entrez.email = "Your.Name.Here@example.org"
	Entrez.api_key = ncbi_key
	logger = logging.getLogger(__name__)
	try:
			from urllib.error import HTTPError # for Python 3
	except ImportError:
			from urllib2 import HTTPError # for Python 2

	count = len(pmids)
	records = []
	attempt = 0

	while attempt < 3:
		attempt += 1
		logger.info('Going to Epost pmid list results')
		try:
			# query pubmed with pmids and post results with ePost
			post_xml = Entrez.epost('pubmed', id=','.join(pmids))
			# read results
			search_results = Entrez.read(post_xml)
			# close the link
			post_xml.close()
			attempt = 4
		except HTTPError as err:
			if 500 <= err.code <= 599:
				logger.warning('Received error from server: %s' % err)
				logger.warning('Attempt %i of 3' % attempt)
				time.sleep(10)
			else:
				raise

	# set paramater values from ePost location to get xml with eFetch
	webenv = search_results['WebEnv']
	query_key = search_results['QueryKey']

	batch_size = 500

	for start in range(0, count, batch_size):
		end = min(count, start+batch_size)
		logger.info('Going to fetch record %i to %i' % (start+1, end))
		attempt = 0
		while attempt < 3:
			attempt += 1
			try:
				# use eFetch to get xml information out of ePost results
				fetch_handle = Entrez.efetch(db='pubmed',
				                             retstart=start, retmax=batch_size,
				                             webenv=webenv, query_key=query_key,
				                             retmode='xml')
				records.extend(str(fetch_handle.read()))
				fetch_handle.close
				attempt = 4
			except HTTPError as err:
				if 500 <= err.code <= 599:
					logger.warning('Received error from server: %s' % err)
					logger.warning('Attempt %i of 3' % attempt)
					time.sleep(10)
				else:
					raise

	pub_list = re.split('<PubmedArticle>', ''.join(records))
	rows = []
	for x in range(1, len(pub_list)):
		# assemble list of publication details
		rows.append(details(pub_list[x], grants))


	pubs_frame = pd.DataFrame(rows, columns=[
	                          'pmid', 'pmcid', 'nihmsid',  'nctid', 'pub_title',
	                          'authors', 'authors_lnames', 'authors_initials',
	                          'orcid', 'authors_affil', 'pub_date', 'journal_short',
	                          'journal_full', 'pubmed_tags', 'pub_type_list', 'exclude',
                              'mesh_major', 'mesh_minor', 'mesh_key'])

	return pubs_frame


### clears text from a webpage element
def clear_text(element):
    length = len(element.get_attribute('value'))
    element.send_keys(length * Keys.BACKSPACE)


def ncbi_login(login, password):
    # set chrome driver options to headless
    options = Options()
    options.headless = True
    driver = webdriver.Chrome(options = options)
    driver.set_window_size(1440, 900)
    driver.get('https://www.ncbi.nlm.nih.gov/myncbi/collections/mybibliography/')
    driver.switch_to.frame(driver.find_element_by_id('loginframe'))
    driver.find_element_by_id('nih').click()
    time.sleep(3)
    driver.find_element_by_id('USER').send_keys(login)
    driver.find_element_by_id('PASSWORD').send_keys(password)
    driver.find_element_by_xpath('//*[@id="CredSelectorNotice"]/div/button').click()
    return driver


def open_my_bib(driver, url, my_bib, delay, long_delay, login, password):
    driver.get(my_bib)
    attempt = 0
    while attempt <= 3:
        try:
            driver.find_element_by_xpath('//*[@id="503-content"]').is_displayed() == True
            time.sleep(long_delay)
            attempt += 1
            if attempt == 3:
                print('My Bibliography webpage is having trouble loading. Please wait a few minutes and trying again...')
                driver.quit()
                return 'Fail'
        except Exception as err:
            attempt = 4
            return driver


def clear_my_bib(driver, delay, logger):
    try:
        select_all = driver.find_element_by_xpath('//*[@id="selectbar"]/div[1]/ul/li[1]/a')
        driver.execute_script("arguments[0].click();", select_all)
        delete_all = driver.find_element_by_id('delete-citations')
        driver.execute_script("arguments[0].click();", delete_all)
        time.sleep(delay)
        #switch to pup up and click 'okay'
        driver.switch_to.alert.accept()
        time.sleep(delay)
        #click 'done' on new popup
        driver.find_element_by_xpath('//*[@id="deleteCitations"]/button').click()
    except Exception as err:
        logger.warning('no citations to delete: %s' %err)


def add_to_my_bib(driver, add_pubs, delay, long_delay, logger):
    # ## Open the 'Add Citations' window
    attempt = 0
    while attempt <= 5:
        try:
            add_citation = driver.find_element_by_xpath('//*[@id="add-drop"]/li[1]/a')
            driver.execute_script("arguments[0].click();", add_citation)
            attempt = 6
        except Exception as err:
            #print('unable to add citations: ', str(err))
            if attempt == 4:
                time.sleep(long_delay*4)
                #print('I really tried to add citations...')
            else:
                time.sleep(delay)
            attempt += 1

    # ## Add citations and search
    attempt = 0
    while attempt <= 3:
        try:
            search_field = driver.find_element_by_xpath('//*[@id="term"]')
            time.sleep(delay)
            attempt = 4
        except Exception as err:
            time.sleep(delay)
            attempt += 1

    attempt = 0
    while attempt <= 3:
        try:
            search_field.send_keys(', '.join(add_pubs))
            time.sleep(delay)
            attempt = 4
        except Exception as err:
            time.sleep(delay)
            attempt += 1

    attempt = 0
    while attempt <= 3:
        try:
            driver.find_element_by_xpath('//*[@id="search-but"]').click()
            attempt = 4
            time.sleep(delay)
        except Exception as err:
            time.sleep(delay)
            attempt += 1

    next_button = True
    while next_button == True:
        time.sleep(delay)

        attempt = 0
        while attempt <= 3:
            try:
                next_button = driver.find_element_by_xpath('//*[@id="nextpage"]').is_displayed()
                attempt = 4
            except Exception as err:
                time.sleep(long_delay)
                attempt += 1

        # ## Select all publications, add, close 'Add Citations' window
        checkboxes = driver.find_elements_by_css_selector("#search-results input[type='checkbox']")
        for checkbox in checkboxes:
            driver.execute_script("arguments[0].click();", checkbox)
        time.sleep(delay)
        if next_button == False:
            #time.sleep(long_delay)
            #add_button = WebDriverWait(driver, long_delay).until(EC.presence_of_element_located((By.XPATH, '//*[@id="add"]')))
            #print('i am trying to click ADD')
            #add_button.submit()
            attempt = 0
            while attempt < 4:
                try:
                    driver.find_element_by_xpath('//*[@id="add"]').click()
                    attempt = 4
                except Exception as err:
                    attempt += 1
                    if attempt == 3:
                        time.sleep(long_delay)
                    else: time.sleep(delay)
        else:
            #print('still trying to click NEXT')
            try:
                driver.find_element_by_xpath('//*[@id="nextpage"]').click()
            except Exception as err:
                print('Failed to click next for some reason: %s', str(err))
                logger.warning('Failed to click next for some reason: %s' %err)
                next_button = False

            time.sleep(long_delay)

    try:
        driver.find_element_by_xpath('/html/body/div[4]/div[1]/button/span[1]').click()
    except Exception as err:
        time.sleep(long_delay)
        driver.find_element_by_xpath('/html/body/div[6]/div[1]/button/span[1]').click()

def scrape_citations(cite, count, grants, driver, delay, long_delay, logger, pub_start):
    load_awards_delay = long_delay
    pmid = re.search('pmid=\\"([0-9]*)\\"', str(cite)).group(1)
    stat = re.search('span class\\=\\"status\\">([\S\s]*)</span>', str(cite)).group(1)
    if str('Complete') in str(stat):
        status = '1'
    elif str('In Process') in str(stat):
        status = '2'
    elif str('Non-compliant') in str(stat):
        status = '3'
    elif str('Not defined') in str(stat):
        status = '4'
    elif str('Citation not in NIHMS') in str(stat):
        status = '4'
    elif str('Exempted') in str(stat):
        status = '5'
    else:
        status = '3'
    cite_award = re.search('<a class=\\"([\S\s]*)\\" href', str(cite)).group(1)
    all_awards = []
    pmc_tag = []
        # if len(re.search('(view8-awards open-award-dialog)', str(awards[x])).group(1))==30:
    if 'view8-awards open-award-dialog' in str(cite_award):
        place = count+1
        # update xpath for award list link
        if place == 1:
            div_place = 'div'
        else:
            div_place = str('div['+str(place)+']')
        awards_list = str('//*[@id="main_content"]/section/div[7]/'+div_place+'/div[2]/div[2]/a')
        attempt = 0
        while attempt < 4:
            try:
                # use updated xpath to open the awards dialog box
                driver.find_element_by_xpath(awards_list).click()
                attempt = 4
            except Exception as err:
                if attempt == 2:
                    time.sleep(long_delay)
                    attempt += 1
                    print('Might have failed to open awards list for pmid: ' + pmid)
                else:
                    time.sleep(delay)
                    attempt += 1
        time.sleep(load_awards_delay)
        # extract the text of the award dialog box
        grant_dialog = driver.find_element_by_xpath('//*[@id="grant-dialog"]').get_attribute('outerHTML')
        all_awards = re.findall('<div class="checked read-only.*?<p>(.*?) -', grant_dialog)
        # Development check for loading time of awards info!!!!!!!!!!!!
        if len(all_awards) == 0:
            print('--May need to wait longer for awards to load, got 0 for pmid: ' + pmid)
            time.sleep(load_awards_delay)
            grant_dialog = driver.find_element_by_xpath('//*[@id="grant-dialog"]').get_attribute('outerHTML')
            all_awards = re.findall('<div class="checked read-only.*?<p>(.*?) -', grant_dialog)
            if len(all_awards) == 0:
                print('---Well maybe there are no awards for pmid: ' + pmid)
            else:
                print('---I waited longer and got them when I tried again.')
        # not sure about the 'in grants' portion below to only list pmc tagged
        # grants that are in the list provided by the user
        for award in all_awards:
            if award in grants:
                pmc_tag.append(award)
        pmc_tag = ', '.join(pmc_tag)
        time.sleep(delay)
        # closes the grant dialog box
        last_try = 'no'
        attempt = 0
        while attempt < 4:
            try:
                driver.find_element_by_id('cancel-association').click()
                attempt = 4
            except Exception as err:
                if attempt == 2:
                    time.sleep(long_delay)
                    attempt += 1
                    print('Might have failed on -cancel association- button for pmid: ' + pmid)
                    last_try = 'yes'
                else:
                    time.sleep(delay)
                    attempt += 1
        if last_try == 'yes':
            time.sleep(long_delay)
            try:
                driver.find_element_by_xpath('/html/body/div[4]/div[1]/button/span[1]').click()
            except Exception as err:
                print('Also may have failed to hit the x for pmid: ' + pmid)
    else:
        pmc_tag = ""

    # assemble the pmid, status, and pmc_tag values into rows of a table
    #row = [pmid, status, pmc_tag, all_awards, (int(count) + int(pub_start))]
    row = [pmid, status, pmc_tag, all_awards, pub_start]
    return row


def scrape_nihms_status(driver, nihms, pmid, delay, long_delay):
    # initialize lists for the nihms status and progress details
    pmc = ''
    reviewer = ''
    files_uploaded = ''
    initial_approval = ''
    nihms_conversion = ''
    final_approval = ''
    pmcid_assigned = ''

    nihms_url = 'https://www.nihms.nih.gov/submission/' + nihms + '/'
    driver.get(nihms_url)
    time.sleep(delay)
    soup = BeautifulSoup(driver.page_source, 'lxml')
    script_info = soup.find_all('div', 'usa-grid')[2]
    # check for a pmcid value and append, else append blank value
    if re.search('PMCID:</dt>\n<dd>([0-9].*?)</dd', str(script_info)) is not None:
        pmc = (re.search('PMCID:</dt>\n<dd>([0-9].*?)</dd', str(script_info)).group(1))

    # check for a reviewer and append, else append blank value
    if re.search('Reviewer:</dt>\n<dd>([A-Za-z].*?)</dd', str(script_info)) is not None:
        reviewer = re.search('Reviewer:</dt>\n<dd>([A-Za-z].*?)</dd', str(script_info)).group(1)

    # package of dates for the progress stages
    script_progress = soup.find_all('div', 'progress')[0]
    all_status = re.findall('<span>\\((.*?)\\ .*?\\)</span>', str(script_progress))
    if len(all_status) >= 5:
        files_uploaded = all_status[0]
        initial_approval = all_status[1]
        nihms_conversion = all_status[2]
        final_approval = all_status[3]
        pmcid_assigned = all_status[4]
    elif len(all_status) == 4:
        files_uploaded = all_status[0]
        initial_approval = all_status[1]
        nihms_conversion = all_status[2]
        final_approval = all_status[3]
    elif len(all_status) == 3:
        files_uploaded = all_status[0]
        initial_approval = all_status[1]
        nihms_conversion = all_status[2]
    elif len(all_status) == 2:
        files_uploaded = all_status[0]
        initial_approval = all_status[1]
    elif len(all_status) == 1:
        files_uploaded = all_status[0]

    row = [pmid, nihms, pmc, reviewer, files_uploaded, initial_approval, nihms_conversion, final_approval, pmcid_assigned]

    return row


def get_nihms(pmids, login, password, delay, long_delay):
    rows = []

    # log into ncbi
    driver = ncbi_login(login, password)

    # navigate to nihms since already logged in to ncbi
    driver.get('https://www.nihms.nih.gov/submission/')
    # click xpath for the NCBI access button - in case of errors, verify this xpath
    driver.find_element_by_xpath('//*[@id="react-app"]/div/div/div[3]/div[3]/a').click()

    # loop through pmids and get nihms status and progress details that are available
    for pmid in pmids:
        search_url = 'https://www.nihms.nih.gov/submission/search/?q=' + pmid
        driver.get(search_url)

        # scrape the search results and see if there's a nihmsid for the pmid
        html = driver.find_element_by_class_name('usa-table-borderless').get_attribute('innerText')

        # initialize lists for the nihms status and progress details
        nihms = ''
        pmc = ''
        reviewer = ''
        files_uploaded = ''
        initial_approval = ''
        nihms_conversion = ''
        final_approval = ''
        pmcid_assigned = ''

        # get the nihmsid
        if re.search('No manuscripts found', html) is not None:
            row = [pmid, nihms, pmc, reviewer, files_uploaded, initial_approval, nihms_conversion, final_approval, pmcid_assigned]
        elif re.search('([0-9].*?)\t', html) is None:
            row = [pmid, 'error', pmc, reviewer, files_uploaded, initial_approval, nihms_conversion, final_approval, pmcid_assigned]
        else:
            nihms = re.search('([0-9].*?)\t', html).group(1)
            row = scrape_nihms_status(driver, nihms, pmid, delay, long_delay)

        rows.append(row)

    driver.close()

    ## package the rows into a data frame
    nihms_frame = pd.DataFrame(rows, columns= ['pmid', 'nihms_id', 'pmc_id', 'reviewer', 'files_uploaded', 'initial_approval', 'nihms_conversion', 'final_approval', 'pmcid_assigned'])

    return nihms_frame


def pacm_login(login, password):
    # set chrome driver options to headless
    options = Options()
    #options.headless = True
    driver = webdriver.Chrome(options = options)
    driver.get('https://auth.nih.gov/CertAuthV2/forms/NIHPivOrFormLogin.aspx')
    driver.set_window_size(1440, 900)
    id_box = driver.find_element_by_id('USER').send_keys(login)
    pass_box = driver.find_element_by_id('PASSWORD').send_keys(password)
    time.sleep(2)
    login_button = driver.find_element_by_id('Image2').click()

    driver.get('https://www.ncbi.nlm.nih.gov/pmc/utils/pacm/')

    driver.find_element_by_xpath('//*[@id="content"]/div/ul/li/a').click()
    driver.switch_to.frame('loginframe')
    driver.find_element_by_xpath('//*[@id="era"]/img').click()
    clear_text(driver.find_element_by_name('USER'))
    login_box = driver.find_element_by_name('USER').send_keys(login)
    pass_box = driver.find_element_by_name('PASSWORD').send_keys(password)
    login_button = driver.find_element_by_xpath('//*[@id="Image2"]').click()
    return driver


def parse_pacm(driver, pacm_root, pmid, grant_list):
    driver.get(pacm_root+pmid)
    soup = BeautifulSoup(driver.page_source, 'lxml')

    nihms_id = ''
    nihms_status = ''
    journal_method = ''
    files_deposited = ''
    initial_approval = ''
    tagging_complete = ''
    final_approval = ''
    initial_actor = ''
    latest_actor = ''

    if re.search('An error has occurred.', soup.text) is not None:
        nihms_status = 'error: pmid not found'
        pacm_grants = ''

    else:
        table = soup.find_all('table')[2].text
        clean_table = unicodedata.normalize("NFKD", table)

        if re.search('NIHMS ID: (.*?)\n', clean_table) is not None:
            nihms_id = re.search('NIHMS ID: (.*?)\n', clean_table).group(1)

        if re.search('Status: (.*?)\n', clean_table) is not None:
            nihms_status = re.search('Status: (.*?)\n', clean_table).group(1)

        if re.search('Method A journal: (.*?)\n', clean_table) is not None:
            journal_method = re.search('Method A journal: (.*?)\n', clean_table).group(1).strip()

        if re.search('Files deposited: (.*?)I', clean_table) is not None:
            files_deposited = re.search('Files deposited: (.*?)I', clean_table).group(1)

        if re.search('Initial approval: (.*?)T', clean_table) is not None:
            initial_approval = re.search('Initial approval: (.*?)T', clean_table).group(1)

        if re.search('Tagging complete: (.*?)F', clean_table) is not None:
            tagging_complete = re.search('Tagging complete: (.*?)F', clean_table).group(1)

        if re.search('Final approval: (.*?)\n', clean_table) is not None:
            final_approval = re.search('Final approval: (.*?)\n', clean_table).group(1)

        if re.search('Initial actor: (.*?)\n', clean_table) is not None:
            initial_actor = re.search('Initial actor: (.*?)\n', clean_table).group(1)

        if re.search('Latest: (.*?)\n', clean_table) is not None:
            latest_actor = re.search('Latest: (.*?)\n', clean_table).group(1)

        if re.search('Associated grants(.*?)Article source:', clean_table) is not None:
            pacm_grants = re.search('Associated grants(.*?)Article source:', clean_table).group(1)
        else:
            pacm_grants = 'not parsed'

    row = [pmid, nihms_id, nihms_status, journal_method, files_deposited, initial_approval, tagging_complete, final_approval, initial_actor, latest_actor, pacm_grants]

    return row


def RC_update_status(pub_comp):
    # Assign REDCap field values to nihms_comm based on date of completion for the nihms steps
    pub_comp.loc[pub_comp['pmc_id'] == '', 'nihms_comm'] = '1'
    pub_comp.loc[pub_comp['pmcid_assigned'] == '', 'nihms_comm'] = '5'
    pub_comp.loc[pub_comp['final_approval'] == '', 'nihms_comm'] = '3'
    pub_comp.loc[pub_comp['initial_approval'] == '', 'nihms_comm'] = '2'
    pub_comp.loc[pub_comp['files_uploaded'] == '', 'nihms_comm'] = '1'

    # Assign REDCap field values to nihms_comm, pmc_status, and author_excluded based on
    # 'Compliant' and 'Excluded' status per PACM scrape
    #pub_comp.loc[pub_comp['pmc_id'].isnull() == False, 'nihms_comm'] = '5'
#    pub_comp.loc[pub_comp['pmc_id'].isin(['']) == False, 'nihms_comm'] = '5'
#    pub_comp.loc[pub_comp['nihms_status'] == 'Compliant', 'nihms_comm'] = '5'
#    pub_comp.loc[pub_comp['nihms_status'] == 'Compliant', 'pmc_status'] = '1'
#    pub_comp.loc[pub_comp['nihms_status'] == 'Excluded', 'nihms_comm'] = '6'
#    pub_comp.loc[pub_comp['nihms_status'] == 'Excluded', 'author_excluded'] = '1'
#    pub_comp.loc[pub_comp['nihms_status'] == 'Excluded', 'pmc_status'] = '5'

    # Assign REDCap field values to journal_method based on Method A documentation in PACM
#    pub_comp['journal_method'] = pub_comp['journal_method'].fillna('')
#    pub_comp.journal_method = pub_comp.journal_method.apply(lambda x: '0' if 'No' in x else x)
#    pub_comp.journal_method = pub_comp.journal_method.apply(lambda x: '1' if 'Yes' in x else x)

    # Update nihms completion dates to importabl REDCap format (YYYY-MM-DD)
#    pub_comp['files_deposited'] = pub_comp['files_deposited'].fillna('')
#    pub_comp.files_deposited = pub_comp.files_deposited.apply(lambda x: x if x in '' else datetime.datetime.strptime(x, "%m/%d/%y").strftime("%Y-%m-%d"))

#    pub_comp['initial_approval'] = pub_comp['initial_approval'].fillna('')
#    pub_comp.initial_approval = pub_comp.initial_approval.apply(lambda x: x if x in '' else datetime.datetime.strptime(x, "%m/%d/%y").strftime("%Y-%m-%d"))

#    pub_comp['final_approval'] = pub_comp['final_approval'].fillna('')
#    pub_comp.final_approval = pub_comp.final_approval.apply(lambda x: x if x in '' else datetime.datetime.strptime(x, "%m/%d/%y").strftime("%Y-%m-%d"))

#    pub_comp['tagging_complete'] = pub_comp['tagging_complete'].fillna('')
#    pub_comp.tagging_complete = pub_comp.tagging_complete.apply(lambda x: x if x in '' else datetime.datetime.strptime(x, "%m/%d/%y").strftime("%Y-%m-%d"))

    return pub_comp


# create function to open a new file to write called MyNCBI-[rppr] then
# write pmids and titles in an amended medline file format
def write_myncbi(pmid_list, title_list, rppr):
    with open('MyNCBI-'+rppr+'.txt', 'w') as myncbi:
        for x in range(len(pmid_list)):
            myncbi.write('PMID '+pmid_list[x])
            myncbi.write('\n')
            myncbi.write('TI '+title_list[x])
            myncbi.write('\n\n')
        myncbi.close
    return 'Complete'
