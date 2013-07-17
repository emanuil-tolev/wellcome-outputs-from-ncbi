import sys

with open(sys.argv[1], 'rb') as f:
    content = f.read()

content = content.splitlines()

data = []
for row in content:
    data.append(row.split(','))

duplicate_count = 0
processed_items = []
for row in data:
    for item in row:

        if item in processed_items:
            duplicate_count = duplicate_count + 1

        processed_items.append(item)

sys.stdout.write(str(duplicate_count))
