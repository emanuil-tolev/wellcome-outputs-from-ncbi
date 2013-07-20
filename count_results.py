import sys

def count_results(filename):
    with open(filename, 'rb') as f:
        content = f.read()
    
    content = content.splitlines()
    
    data = []
    for row in content:
        data.append(row.split(','))
    
    count = 0
    for row in data:
        for item in row:
            count = count + 1

    return count

if __name__ == '__main__':
    sys.stdout.write(str(count_results(sys.argv[1])))
