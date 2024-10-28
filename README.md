# nmdc-client

* GitHub: https://github.com/jeffbaumes/nmdc-client
* Documentation: https://jeffbaumes.github.io/nmdc-client

This module provides functions to interact with the NMDC API,
including retrieving metadata into dataframes and linking across schema objects.

## Getting started

To install:

```bash
pip install git+https://github.com/jeffbaumes/nmdc-client.git
```

To use:

```python
from nmdc_client import NmdcClient
from nmdc_schema import nmdc

client = NmdcClient()
studies = client.find(nmdc.Study)
```
