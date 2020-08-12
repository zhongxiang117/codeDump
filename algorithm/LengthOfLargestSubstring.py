def lengthOfLongestSubstring(s):
    # a hashtable to append substrings
    dicts = {}
    
    maxlen = start = 0
    for i, value in enumerate(s):
        
        # only for non-repeats
        if value in dicts:
        
            # get real-index, python starts at zero, thus add 1
            new = dicts[value] + 1
            
            # avoid many repeats, only the Final char will be matched
            if new > start:
                start = new
        # length
        nm = i - start + 1
        
        # calculate new largest length
        if nm > maxlen:
            maxlen = nm
        
        # update table
        dicts[value] = i
    return maxlen
