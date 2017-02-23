
# Houses functions shared across scripts for the GenBank processing pipeline.

# Function to print the product name to a file. Accepts the product name
# and the handle for the outfile. Will process multi-line name if the
# characters-per-line exceeds 80.
def printName(name,out): # need some handling in case of multi-line
    base = " " * 21 + "/product=" # first line
    blank_base = " " * 21 # non-first lines

    # Here replace some common errors found in product names
    name = name.replace('sulphide','sulfide')
    name = name.replace('sulphur','sulfur')
    name = name.replace(' fibre ',' fiber ')

    if len(name) < 48: # good to go, can print on one line
        final = '%s"%s"\n' % (base,name)
        out.write(final)

    # If greater than that length, need to do a multi-line print. Note
    # that this approach does not absolutely maximize line length, but 
    # arbitrarily does not permit this region to extend beyond the ~70 
    # character mark. 
    else: 

        max_len = 47
        words = name.split(' ')
        multi_line = [] # store all the lines in the product qualifier
        single_line = [] # store just the current line being built

        for x in words:
            if not (len(' '.join(single_line)) + len(x)) > max_len: # if i can fit, add it
                single_line.append(x)
            else: # too long, add to multi-line and build a new one
                multi_line.append(' '.join(single_line))
                single_line = []
                single_line.append(x)

        multi_line.append(' '.join(single_line))

        for j in range(0,len(multi_line)):
            line = ""
            if j == 0: # first line
                line = '%s"%s\n' % (base,multi_line[j])
            elif j == (len(multi_line)-1): # last line
                line = '%s%s"\n' % (blank_base,multi_line[j])
            else: # middle line
                line = '%s%s\n' % (blank_base,multi_line[j])
            out.write(line)