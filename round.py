import sys

with open(sys.argv[1], "r") as stream:
    for line in stream:
        try:
            n = float(line)
            print(round(n * 10) / 10)
        except:
            print()
