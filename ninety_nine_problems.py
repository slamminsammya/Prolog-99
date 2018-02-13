import random
import math
import time

## Ninety-nine prolog problems. Meant to practice basic problem solving.
## 1) Find the last element of a list

def get_last(ell):
    return ell[-1]

## 2) Find the second to last element of a list

def second_last(ell):
    return ell[-2]

## 3) Find the k'th element of a list, where first element is 1'st.

def kth_element(ell, k):
    return ell[k - 1]

## 4) Find the number of elements in a list.

def num_elements(ell):
    return len(ell)

## 5) Reverse a list (python has built in list reversal but
## we gonna implement it ourselves.

def reverse(ell):
    
    reversed_list = [ell[j] for j in range( -1, -len(ell) - 1)]
    return reversed_list

## 6) Find out if a list is a palindrome.

def palindrome(ell):
    
    reverse_ell = reverse(ell)
    for j in range(len(ell)):
        if ell[j] != reverse_ell[j]:
            return False
    return True

## 7) Flatten a nested list structure.

def flatten(ell):
    
    if type(ell[0]) == list:
        first_part = flatten(ell[0])

    else:
        first_part = ell[0]

    return first_part + flatten(ell[1:])

## 8) Eliminate consecutive duplicates of list elements

def compress(ell):

    for j in range(len(ell) - 1):
        if ell[j] == ell[j + 1]:
            new_list = ell[ :j + 1] + ell[j + 2: ]
            return compress(new_list)
    return ell

## 9) Pack consecutive duplicates of list elements into sublists.

def pack(ell):

    output = []
    counter = 1
    
    for j in range(len(ell) - 1):
        
        if ell[j] != ell[j + 1]:
            block = counter * [ ell[j] ]
            output.append(block)
            counter = 1
            
        else:
            counter += 1
            
    last_block = counter * [ ell[-1] ]
    output.append(last_block)

    return output
    
## 10) Run-length encoding of a list.

def encode(ell):
    packed_list = pack(ell)
    
    for j in range(len(packed_list)):
        multiplicity = len(packed_list[j])
        element = packed_list[j][0]
        packed_list[j] = [multiplicity, element]
        
    return packed_list

## 11) Modify run-length encoding.

def modified_encode(ell):

    encoded = encode(ell)
    
    for j in range(len(encoded)):
        multiplicity = encoded[j][0]
        
        if multiplicity == 1:
            encoded[j] = encoded[j][1]
            
    return encoded

## 12) Decode a run-length enccoded list as in problem 11.

def decode(encoded):

    output = []

    for j in range(len(encoded)):
        if type(encoded[j]) == list:
            multiplicity = encoded[j][0]
            element = encoded[j][1]
            chunk = multiplicity * [element]

        else:
            chunk = [ encoded[j] ]
        output += chunk

    return output

## 13) Run-length encoding of a list (direct solution)

def alternate_encode(ell):

    output = []
    counter = 1
    
    for j in range(len(ell) - 1):
        
        if ell[j] != ell[j + 1]:
            
            if counter == 1:
                block = ell[j]
                
            else:
                block = [counter, ell[j]]
                
            output.append(block)
            counter = 1
            
        else:
            counter += 1
            
    if counter == 1:
        last_block = ell[-1]
        
    else:
        last_block = [counter, ell[-1]]
        
    output.append(last_block)

    return output

## 14) Duplicate the elements of a list.

def duplicate(ell):

    ouptut = []
    
    for x in ell:
        output += [x, x]
        
    return output

## 15) Duplicate the elements of a list a given number of times.

def duplicate(ell, n):

    output = []
    
    for x in ell:
        chunk = n * [x]
        output += chunk

    return output

## 16) Drop every n'th element of a list.

def drop(ell, n):
    
    if n > len(ell):
        return ell

    else:
        return ell[ :n - 1] + drop(ell[n: ], n)

## 17) Split a list into two parts, the length of the first part is given.

def split(ell, k):

    return [ell[:k], ell[k:]]

## 18) Extract a slice from a list. I.e. a list containing i'th element of
## ell up to the j'th element, inclusive. Here we start counting at 1.

def interval(ell, i, j):
    return ell[i - 1: j]

## 19) Rotate a list n places to the left.

def left_rotate(ell, n):
    return ell[n - 1:] + ell[: n - 1]

## 20) Remove the k'th element from a list. (Solve without using built in
## functions where possible.)

def remove(ell, k):
    return ell[: k - 1] + ell[k:]

## 21) Insert an element at a given position in a list.

def insert(ell, element, k):
    return ell[: k - 1] + [element] + ell[k - 1: ]

## 22) Create a list containing all integers within a given range.

def make_range(j,k):
    number = j
    ell = []
    
    while number <= k:
        ell.append(number)
        number += 1

    return ell

## 23) Extract a given number of randomly selected elements from a list.

def extract_randomly(ell, n):

    output = []
    count = n
    current_list = ell
    
    while count > 0:
        random_index = random.randrange(len(ell))
        
        output.append(ell[random_index])
        ell = remove(ell, random_index + 1)
        ## We remove at random_index + 1 since the problem uses a different
        ## counting convention than python, starting at one, and the remove 
        ## function follows this convention for its argument.
        count = count - 1

    return output

## 24) Draw N different random numbers from the interval 1,..., M

def lotto(tickets, m):
    interval = make_range(1, m)
    
    if n <= m:
        return extract_randomly(interval, n)
    
    else:
        print 'Invalid arguments. Number of choices cannot exceed size of interval.'

## 25) Generate a random permutation of the elements of a list.

def random_perm(ell):
    return extract_randomly(ell, len(ell))

## 26) Generate all combinations of k distinct objects from n distinct objects.

def combinations(ell, k):

    if k == 1:
        singletons = [ [x] for x in ell]
        return singletons

    elif k > len(ell):
        return []

    else:
        first_element = ell[0]
        do_contain = [ [first_element] + x for x in combinations(ell[1:], k - 1) ]
        dont_contain = combinations(ell[1:], k)
        return do_contain + dont_contain
    
    pass

## 27) Group the elements of a set into disjoint subsets. (Partition problem)
## Here we do care about the order of the partition. A | B different than
## B | A

def partition(ell, k):

    if k == 1:
        return [ell]

    elif k > len(ell):
        return []

    else:
        partitions = []
        
        for size in range(1, len(ell) - (k - 1) + 1):
            first_subgroups = combinations(ell, size)
            
            for first in first_subgroups:
                remaining_list = [elements for elements in ell if elements not in first]
                remaining_partitions = partition(remaining_list, k - 1)
                resulting_partitions = [ [first] + [x] for x in remaining_partitions]
                partitions += resulting_partitions
                
        return partitions

    pass
            
            
## 28) Sorting a list of lists according to length of sublists, short to long.
## As a second version, sort them according to their length frequency.

def list_sort_1(ell):

    ## Implement a merge sort, but compare lengths of elements.
    
    if len(ell) == 1:
        return ell
    
    elif len(ell) == 2:
        first = ell[0]
        second = ell[1]
        if len(first) > len(second):
            return [ell[1], ell[0]]
        else:
            return ell

    else:
        midpoint = len(ell) / 2
        sublist_1 = ell[:midpoint]
        sublist_2 = ell[midpoint:]
        sorted_sublist_1 = list_sort_1(sublist_1)
        sorted_sublist_2 = list_sort_1(sublist_2)
        merged_list = []
        
        while len(merged_list) < len(ell):

            if len(sorted_sublist_1) == 0:
                merged_list += sorted_sublist_2

            elif len(sorted_sublist_2) == 0:
                merged_list += sorted_sublist_1

            elif len(sorted_sublist_1[0]) < len(sorted_sublist_2[0]):
                merged_list.append(sorted_sublist_1.pop(0))

            else:
                merged_list.append(sorted_sublist_2.pop(0))

    return merged_list

## Now implement a quicksort, using left pivot.

def list_sort_2(ell):

    if len(ell) == 0 or len(ell) == 1:
        return ell

    else:
        pivot = ell[0]
        sublist_lessthan = []
        sublist_greaterthan = []
        
        for element in ell[1:]:
            
            if len(element) < len(pivot):
                sublist_lessthan.append(element)

            else:
                sublist_greaterthan.append(element)

        sorted_list = list_sort_2(sublist_lessthan) + [pivot] +list_sort_2(
            sublist_greaterthan)
        
    return sorted_list

## Now sorted by length frequency

## Subroutine to take a sorted list of lists and group its elements
## by length. Length of an element of output is therefore frequency of
## given length.

## Subroutine to output longest consecutive sublist of elements of same size,
## beginning with the first element.

def head_samesize(ell):
    
    index = 1
    while index < len(ell):
        if len(ell[index - 1]) != len(ell[index]):
            return ell[:index]

        index += 1
        
    return ell

    
def group_by_size(ell):

    ## List ell must already be sorted by size.
    output = []
    
    while len(ell) > 0:
        group = head_samesize(ell)
        output.append(group)
        ell = ell[len(group):]
    return output

def sort_by_frequency(ell):

    sorted_list = list_sort_2(ell)
    list_of_sizegroups = group_by_size(sorted_list)
    sorted_list_of_sizegroups = list_sort_2(list_of_sizegroups)
    sorted_by_frequency = []

    for group in sorted_list_of_sizegroups:
        sorted_by_frequency += group
        
    return sorted_by_frequency


## Determine whether a given integer is prime. Let us sieve like
## Eranthoses.

## First define a function which, given a list and an index as input,
## changes all entries whose index is a multiple (not zero or one) of
## given index to false.

def sieve(ell, index):

    current_index = 2 * index
    
    while current_index < len(ell):
        ell[current_index] = False
        current_index += index

    return ell

## Given a list sieved and a starting point, finds the first non-sieved
## element following starting point.

def next_prime(ell, last_prime):

    index = last_prime + 1
    while ell[index] == False:
        index += 1
        
    return ell[index]

def sieve_integers(bound):

    integers = range((bound + 1) ** 2)
    integers[0], integers[1] = False, False
    current_prime = 2
    
    while current_prime < bound:
        sieve(integers, current_prime)
        current_prime = next_prime(integers, current_prime)
        
    return integers
    

## Uses previous subroutines. Composites get set to 'False', so we start
## setting 0, 1 to False, and success is if our desired input is an integer,
## i.e. was not set to False. 

def isprime(n):
    
    # We sieve by changing composites to the Boolean false because why not!
    bound = int(math.floor(math.sqrt(n)) + 1)
    integers = sieve_integers(bound)

    return type(integers[n]) == int


## 32) Euclid's algorithm.

def gcd(a, b):
    if a > b:
        a, b = b, a

    if a == 0:
        return b

    else:
        qa = a
        remainder = b - a
        while remainder > a:
            qa += a
            remainder = b - qa

    return gcd(a, remainder)

## 33) Determine whether two positive integers are coprime.

def coprime(a, b):
    return gcd(a, b) == 1

## 34) Calculate Euler totient function.

def totient(m):
    
    phi = 0
    for n in range(1, m):
        if coprime(n, m):
            phi += 1
    return phi

## 35) Determine the prime factors of a positive integer with
## multiplicity.

def get_factors(n):

    if n == 1 or n == 0:
        return []
    
    else:
        bound = int(math.floor(math.sqrt(n)) + 1)
        primes_to_bound = [p for p in sieve_integers(bound) if type(p) == int
                           and p < bound]

        for prime in primes_to_bound:
            if gcd(prime, n) > 1:
                return [prime] + get_factors(int(n / prime))
        return [n]
    
def factor(n):
    
    return encode(get_factors(n))

## 36) Determine prime factors with multiplicity (accidentally done as 35.

## 37) Calculate totient function using your factorization.

def phi(prime_multiplicity):

    prime = prime_multiplicity[1]
    multiplicity = prime_multiplicity[0]
    return (prime - 1) * prime ** (multiplicity - 1)

def totient_2(n):

    totient = 1
    factorization = factor(n)
    for prime_multiplicity in factorization:
        totient = totient * phi(prime_multiplicity)
    return totient

## 38) Compare the two totient algorithms.

def advantage(n):

    good_start = time.time()
    totient_2(n)
    good_finish = time.time()
    good_time = good_finish - good_start

    bad_start = time.time()
    totient(n)
    bad_finish = time.time()
    bad_time = bad_finish - bad_start

    return 'The efficient algorithm had an advantage of %s seconds for n =%s.' % (
        bad_time - good_time, n)

## 39) A list of prime numbers.

## No desire to be efficient. This is extraordinarily redundant.

def prime_list(lower_bound, upper_bound):

    return [p for p in range(lower_bound, upper_bound + 1) if isprime(p)]

## 40) Goldbach's conjecture: Find two prime numbers that sum to a given even.

def goldbach_pair(n):

    upper_bound = int(math.floor(n / 2) + 1)
    primes = prime_list(3, upper_bound)
    for prime in primes:
        if isprime(n - prime):
            return [prime, n - prime]
    return None

## 41) List of Goldbach pairs: Given a range of integers, print a list of
## all even numbers and their Goldbach composition.

def list_goldbachs(a, b):
    for n in range(a,b):
        if n % 2 == 0:
            print goldbach_pair(n)
    pass

## Problems skip to 46?
## 46) 
    


