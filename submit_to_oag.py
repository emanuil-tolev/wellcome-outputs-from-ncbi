import sys
import requests
from datetime import datetime

with open(sys.argv[1], 'rb') as f:
    content = f.read()

content = content.splitlines()

data = []
for row in content:
    r = requests.get('http://oag.cottagelabs.com/lookup/' + row)
    if r.status_code != 200:
        print datetime.now(), 'OAG returned', r.status_code, '; First 100 resp. chars:', r.text[:100]
