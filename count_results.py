import sys

with open(sys.argv[1], 'rb') as f:
    content = f.read()

content = content.splitlines()

data = []
for row in content:
    data.append(row.split(','))

count = 0
for row in data:
    for item in row:
        count = count + 1

sys.stdout.write(str(count))
