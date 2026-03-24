# Publication Compliance Tracker

The **Publication Compliance Tracker** is a Python-based automation tool designed to retrieve, process, and manage publication compliance data across multiple NIH-related systems. It supports automated scraping, API queries, data extraction, REDCap integration, and compliance assessment for grant-funded publications.

This tool is especially useful for research administration offices, CTSAs, and compliance personnel responsible for monitoring NIH Public Access Policy compliance.

---
## Features
- PubMed publication lookup via grant numbers
- PMC compliance scraping via MyNCBI Bibliography
- NIHMS manuscript status scraping (Selenium)
- iCite bibliometrics API querying
- Altmetric API querying
- PACM scraping
- REDCap integration
- Automated login to NIH systems

---
## Architecture
- `pub_command_line.py`: Command-line interface and workflow manager
- `pub_comp_lib.py`: Core scraping/processing library
- `config.py`: User configuration file

---
## Requirements
See `requirements.txt`.

---
## Installation
```bash
git clone <repository-url>
cd publication-compliance-tracker
pip install -r requirements.txt
```

---
## Configuration
Copy the template configuration file:
```bash
cp config_blank.py config.py
```
Update grant list, REDCap settings, MyNCBI login, and API keys.

---
## Usage
```bash
python pub_command_line.py [database] [timeframe]
```

**Database options:** all, compliance, pubmed, pmc, nihms, icite, altmetric  
**Timeframes:** all, refresh, current

Example:
```bash
python pub_command_line.py all all
```

---
## Outputs
Generated CSV files include:
- `batch_pubmed_query_output.csv`
- `dev-pmc_query_output.csv`
- `dev-nihms_query_output.csv`
- `dev-icite_query_output.csv`
- `dev-altmetric_query_output.csv`

If REDCap integration is enabled, these results are imported automatically.

---
## Troubleshooting
- Ensure ChromeDriver matches your Chrome version
- Verify REDCap permissions
- Increase delay values if MyNCBI loads slowly

---
## File Structure
```
.
├── pub_command_line.py
├── pub_comp_lib.py
├── config_blank.py
├── config.py
└── requirements.txt
```

---
## License
Insert licensing terms here.
