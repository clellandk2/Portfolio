# Problem 152

# This is a solution to Project Euler problem #152, which I've included as a screenshot elsewhere in the directory. This solution is based on a
# derivation I did (also in the directory) which shows that prime numbers and their multiples must appear in certain valid combinations in order
# to even possibly be a part of a valid answer. The criteria is that p^2 must divide a, where p is the prime and a is the numerator after summing up
# all the 1/n^2 terms and factoring out p^2 from the denominator.

# The code uses this result (p^2 | a) to find all possible valid combinations of prime factors and their multiples for primes larger than 3. 
# Finding valid combinations for multiples of 3 proves expensive (due to both a larger number of multiples and a smaller p^2), so instead numbers
# that are strictly multiples of 2 and 3 are handled by tabulating all possible partial sums that can come from these terms (2^14 possibilities),
# which are then stored in a dictionary, along with their corresponding numbers (including duplicates for degenerate sums). 

# The valid combinations from larger prime factors are then combined together in all possible overlaps and the partial sums of each are computed. 
# These partials sums are then subtracted from the target value (1/2) and that number is looked up in the dictionary. If it's a valid entry, 
# then one (or more) correct answers have been found. The dictionary value stores the number of ways that particular sum can be achieved from 
# factors of 2 and 3, so that multiplicity is added to the total answer count. 

# A few other notes: 
# 1) It is a correlary of my proof that prime numbers cannot appear by themselves without other multiples, so we don't need to check prime factors
#    larger than 37 (since n = 80). 

# 2) By computing the sum of all the 1/n^2 terms ranging from 3 to 80 it can be easily shown that 2 must be included in any valid answer, so this 
#    code includes 2 automatically by treating the target sum as 0.25 and excluding 2 from the search space entirely. 

# 3 Runtime:
# The runtime is dominated by generating all valid combinations of multiples of 5 and by the precomputation for the multiples of 2 and 3.  
# For maxNum = 80, this is roughly 2^14 + 2^14 = 2^15 iterations. 
# More generally, the complexity grows exponentially with the number of contributing primes.





from itertools import combinations

def getFrac(p, nums, valid): 
    """ Takes a prime (p), a list containing some of its multiples (nums), and a list (valid) which stores the combo (nums) 
    if it satisfies the p^2 | a criteria. 
    The function evaluates the p^2 | a criteria and does any necessary appending, but returns nothing. """
    
    if (len(nums) == 0):
        return                      # do not add an empty combo 
    reduced = []
    for n in nums:
        if (n%p != 0):              # check just to be sure all the multiples are indeed multiples of p
            print("Problem.")
        reduced.append((n/p)**2)    # list of denomiators in the 1/n^2 series once p^2 has been factored out, as seen in the derivation 
    a = 0
    b = 1                           # represents the fraction a/b - we need to work in fractions, not decimals for this 
    for i in reduced:
        c = a*i + b
        d = b*i                     # the new fraction is c/d 
        a = c
        b = d
    if (a%(p*p) == 0):              # checking p^2 | a 
        valid.append(nums)

def powerSet(mults):
    """Generates all possible subsets (the powerset) of the given list 'mults'. Returns empty set if given an empty set. """
    
    all_subsets = []
    # combinations(mults, r) gives all subsets of length r
    for r in range(len(mults) + 1):
        for combo in combinations(mults, r):
            all_subsets.append(list(combo))
    return all_subsets


def findCombos(p, banned, maxNum):   
    """ This function takes a prime number p (from the primes list in main) and ultimately returns all possible combinations of its multiples 
    that satisfy the inclusion criteria. #p is the prime, banned is a list of three numbers which can be excluded for reasons given in main. """ 
    
    mults = []  
    i = 1 
    while (i*p < maxNum + 1):       # iterates over all multiples up to maxNum
        if (not (i*p in banned)):   # Blocks the banned numbers from being added in
            mults.append(i*p)
        i += 1
    lists = powerSet(mults)         # generates the powerset, all possible combinations of these multiples 
    valid = []                      # This will ultimtaely be a list of combos that meet the inclusion criteria 
    for combo in lists:
        getFrac(p, combo, valid)    # takes each combo and evaluates the inclusion criteria (p^2 | a). Valid combinations are appended onto valid
    return valid


def finalSet(masterList):
    """ Recursively combines all the data in masterList (which is a list of lists of lists) into a list of lists collectively containing 
    # all the possible combinations of numbers that have at least one prime factor greater than 3. """
    
    workingSet = []
    options = masterList.pop(-1)         # Retrives all possible combinations of multiples for a specific prime, p (pop from end for O(1) speed)
    options.append([])                   # this represents the option of not including any multiples of p, since the list only contains valid combinations
    if (len(masterList) == 0):
        return options                   # base case, if this is the last remaining entry (the first prime listed in this case), just return those options
    smallers = finalSet(masterList)      # retrives all options for combinations of the primes which appear further down the recursion tree (the larger ones)
    for opt in options:
        for s in smallers:  #s is [a,b,c] 
            workingSet.append(opt + s)   # for each valid combination for the current prime (opt), combine it with every option for the other primes 
    return workingSet


def sortEm(nums):
    """This function takes a list of arrays, and sorts and removes duplicates for each array """
    
    final = []
    for x in range(len(nums)): 
        n = nums[x] 
        count35 = 0
        count70 = 0
        n.sort()                    # sorting
        x = 0 
        removed = 0
        while (x < (len(n)-1)):
            myNum = n[x]
            if (myNum == n[x+1]):
                removed = n.pop(x)  # removes duplicates
            else:
                x += 1
                
            # This next part handles whether sets overlap correctly (and it only actually comes up for 5 and 7).
            # The idea is if 35 is in the list of all multiples of 7, then it must also be in the list of all multiples of 5, by definition. 
            # With the way the code's written, this shows up as a duplicate here. But sometimes, naievely combining valid combinations for 5 and 7
            # will yield one list that includes 35, and one that doesn't, and that's clearly a nonsensical scenario since both are complete lists 
            # of all the prime's multiples in the current combination of numbers being tried.
            
            if myNum == 35: 
                count35 = 35 - removed 
            elif myNum == 70: 
                count70 = 70 - removed 
        if (count35 == 0):
            if (count70 == 0):
                final.append(n)      
    return final



def rigorousCheck(soln, printSolutions):  
    """A rigorous mathematical check of all solutions that come very close to the correct sum."""
    
    soln.append(2)                        # slightly slower computationally here, but needed to print full solution sets if printSolutions == True
    soln.sort()
    a = 0
    b = 1                                 # representing the fraction a/b  
    for i in soln:
        c = a*i*i + b                     # Exact algebraic computation 
        d = b*i*i                         # the new fraction is c/d 
        a = c
        b = d                             # update a/b
    if (a*2 == b):                        # has to sum to 1/2 by exact arithmetic 
        if (printSolutions == True):      # prints answers if you want
            print(soln)
        return soln
    else: 
        print("FAILED: Likely a machine precision issue. A combination was found that is within 10^-14 of the correct sum but not exact.")
    return []


def smallPrimes(maxNum):
    """ Precomputes all 2^14 possible partial sums from the terms that are only multiples of 2 and 3. """
    
    twosAndThreesFull = [4,8,16,32,64,3, 6, 9, 12, 18, 24, 36, 48, 72]
    twosAndThrees = []
    for n in twosAndThreesFull:
        if (n <= maxNum):
            twosAndThrees.append(n)
    possibilities = powerSet(twosAndThrees)
    partials = {}
    for nums in possibilities:
        total = 0
        for n in nums:
            total += 1/(n**2)                 # Computes sum 
        total = round(total, 14)              # Rounding to high precision, but avoiding machine error. 
        if total in partials:
            theValue = partials[total]
            theValue[0] = theValue[0] + 1
            theValue[1].append(nums)
            partials[total] = theValue
        else:
            partials[total] = [1,[nums]]      # The value used in this dictionary is an integer that indicates the number of ways this sum can be obtained, 
                                              # and a list of which combinations of numbers achieve that. Usually just [1,[a,b,c...]], but sometimes more complicated.
    return partials 



def main(maxNum,printSolutions):
    """main function.
    maxNum - max denomiator, highest n allowed in the problem description 
    (note: using values higher than 80 will require correctly adding to "primes", and "twosAndThrees")
    printSolutions - print the valid solutions in list format. Slower runtime for obvious reasons. 
    This code is designed to solve the problem specifcally for maxNum = 80, 
    but I wanted to make it a bit generalizable for good practice. """ 
    
    primes = [37, 31, 29, 23, 19, 17, 13, 11, 7, 5]   # remember that 2 and 3 are handled separately
    masterList = []                                   # For each entry, stores a list of arrays collectively storing all of the valid multiple combinations for a particular prime. 
    
    banned = [55,77,65]
    
    # This is an optimization I put in later. You can observe in the code's output that since there are no valid combinations involving factors of 11,
    # 55 does not need to be included as a possible multiple of 5. The same is true for 77. 
    # Likewise, since there is only one valid combination of factors for 13 (13, 39, and 52), and it does not include 65, we may exclude 65 as well. 

    print("# of valid combinations:")
    for p in primes:
        print("Prime: ", p)
        result = findCombos(p, banned, maxNum)    # Finds all possible combinations of multiples for a given prime, p.  
        print(len(result))                        # I think this is a helpful number to give more intuition about the search space
        masterList.append(result)                 # save combos 
    candidates = finalSet(masterList)             # appends the lists for each prime so far together 
    candidates = sortEm(candidates)               # sorts + removes duplicates 
    valid = smallPrimes(maxNum)                   # tabulates all possible partial sums from the set of numbers that only have 2's and 3's as factors.
    
    count = 0                                     # counts the number of solutions found 
    for i in range(len(candidates)):
        currentArray = candidates[i]              # iterate over each possible combination (candidate) 
        partialSum = 0
        for x in currentArray:
            partialSum += 1/(x**2)                # computes the partial sum from the larger prime contributions 
        difference = 0.25 - partialSum            # calculates the remaining difference 
        difference = round(difference, 14)        # Rounding to high precision, but avoiding machine error, and exact check is done later anyway.
        if difference in valid:                   # searches the remaining difference in the dictionary to find valid answers 
            dictEntry = valid[difference]
            multiplicity = dictEntry[0]  
            for k in range(multiplicity):
                solution = currentArray + dictEntry[1][k]
                if (len(rigorousCheck(solution, printSolutions)) == 0):  # algebraically checks each solution, this is not too expensive since there aren't very many at this point
                    count -= 1                                           # excludes false combinations from the total count
            count += multiplicity

    print("\nDone. \n\nFinal Count =", count) 
    return count

if __name__ == "__main__":
    main(80, False)
    #main(45, True) 
    

