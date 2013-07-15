def doi():
    yield '10.1371/journal.pone.0035089'

def wellcome_ncbi_sim_result():
    import get_wellcome_ncbi_objects as t
    data = t.OAGPrep()

    for i in range(0, 75000):
        data.add(doi().next())

    t.to_file('test_data.txt', data)

def main(argv=None):
    if not argv:
        import sys
        argv = sys.argv

    wellcome_ncbi_sim_result()

if __name__ == '__main__':
    main()
