import sys # needed to accept args from command line

########## Initial Design Notes
# two args: database and time 
# database default is 'all', length check to throw error if more than two sent
# time default is 'all', now uses config.start for current grant cycle
# help shows the list of accepted args
# database list: all, pubmed, pmc, nihms, icite, altmetric
# time list: all, now
###########


# total arguments
n = len(sys.argv)
print("Total arguments passed:", n)
 
# Arguments passed
print("\nName of Python script:", sys.argv[0])
 
print("\nArguments passed:", end = " ")
for i in range(1, n):
    print(sys.argv[i], end = " ")

