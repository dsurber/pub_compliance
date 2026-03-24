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

This project is licensed under the MIT License.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the “Software”), to deal
in the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

