import sys
import os

current_query = ""
hits = []

with open(sys.argv[1], 'w') as out_file:
	for line in sys.stdin:
		temp = line.split( )
		if not current_query:
			current_query = temp[0]
			hits.append([temp[2], temp[3], temp[4]])
		else:
			if temp[0] != current_query:
				no_self = hits[1:]
				pass_fail = []
				if not no_self:
					out_file.write(current_query + "\n")
					out_file.flush()
					os.fsync(out_file)
				
				current_query = temp[0]
				hits = []
				hits.append([temp[2], temp[3], temp[4]])
			else:
				hits.append([temp[2], temp[3], temp[4]])