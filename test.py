from simulateMass import massSplit, transferSplitDict

compoundL = ['mn', 'R2', 'R3', 'R4']
aDict = {}

massSplit(compoundL,  aDict)

print aDict

compoundLL = []
transferSplitDict(aDict, compoundLL, [])

compoundLL.sort()
print compoundLL
